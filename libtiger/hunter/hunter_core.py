# -*- coding: utf-8 -*-
#
# Copyright (C) 2021-2022  LMAI_team @ TU Dresden:
#     LMAI_team: Zhixu Ni, Maria Fedorova
# Copyright (C) 2016-2021  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
#     SysMedOs_team: Zhixu Ni, Georgia Angelidou, Maria Fedorova
#
# LPPtiger2 is Dual-licensed
#     For academic and non-commercial use: `AGPL License V3` Please read more information by the following link:
#         [GNU Affero General Public License] (https://www.gnu.org/licenses/agpl-3.0.en.html)
#     For commercial use:
#         Please contact the LMAI_team by email.
#         https://tu-dresden.de/med/mf/zml/forschungsgruppen/fedorova/
#
# Please cite our publication in an appropriate form:
#     Ni, Zhixu, Georgia Angelidou, Ralf Hoffmann, and Maria Fedorova.
#     LPPtiger software for lipidome-specific prediction
#     and hunter of oxidized phospholipids from LC-MS datasets
#     Scientific Reports 7, Article number: 15138 (2017).
#     DOI: 10.1038/s41598-017-15363-z
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@tu-dresden.de
#     LPPtiger2 repository: https://github.com/LMAI-TUD/lpptiger2
#
from collections import ChainMap
import json
import math
from multiprocessing import Pool
import os
import time
from datetime import datetime

import pandas as pd

from libtiger.hunter.score_generator import (
    hunt_isotope_score,
    hunt_rank_score,
    hunt_msp_score,
    hunt_fingerprint_score,
    hunt_snr_score,
    get_total_score,
)
from libtiger.reader.space_builder import (
    build_search_space,
    # build_residue_ranges,
    build_search_ranges,
)
from libtiger.reader import mzml_spectra_reader, star_spectra_reader
from libtiger.writer.plotter import plot_spectra
from libtiger.writer.reporter import Result, LogPageCreator
from libtiger.writer.exporter import export_to_xlsx


def hunt_mzml_ms2(
    r,
    neg_pr_df,
    pos_pr_df,
    neg_residue_space,
    pos_residue_space,
    selection_window,
    s_info_df,
    sum_spec_dct,
    ms1_ppm,
    ms1_abs_threshold,
    min_isotope_score,
    lipid_space,
    min_rank_score,
    sum_xic_df,
    ms2_ppm,
    output_fig_abs_path,
    dpi=300,
):
    ms2_result_dct = {}

    sum_isotope_score_info_dct = hunt_isotope_score(
        r,
        neg_pr_df,
        pos_pr_df,
        neg_residue_space,
        pos_residue_space,
        selection_window,
        s_info_df,
        sum_spec_dct,
        ms1_ppm,
        ms1_abs_threshold,
        min_isotope_score,
    )

    for chg_formula in sum_isotope_score_info_dct:
        isotope_score_info_dct = sum_isotope_score_info_dct[chg_formula]
        pr_ms2_df = isotope_score_info_dct.get("ms2_df", pd.DataFrame())
        pr_residue_space = isotope_score_info_dct.get("pr_residue_space")
        # print(isotope_score_info_dct)
        if (
            isotope_score_info_dct.get("isotope_score", 0) >= min_isotope_score
            and not pr_ms2_df.empty
        ):
            pr_formula_candidates = lipid_space.get(
                isotope_score_info_dct.get("charged_formula", ""), {}
            )
            sum_rank_score_info_dct = hunt_rank_score(
                pr_candidates=pr_formula_candidates,
                ms2_df=pr_ms2_df,
                res_space_info=pr_residue_space,
                weight_df=pd.DataFrame(),
                rank_score_threshold=min_rank_score,
            )
            if sum_rank_score_info_dct:
                lib_mz = isotope_score_info_dct.get("mz", 0)
                xic_df = sum_xic_df.query(
                    f"{(1 - ms1_ppm * 0.000001) * lib_mz} <= mz <= {(1 + ms1_ppm * 0.000001) * lib_mz}"
                ).copy()
                xic_df.drop(columns="mz", inplace=True)
                xic_data = list(
                    zip(
                        [round(rt, 4) for rt in xic_df["rt"].to_list()],
                        [round(i, 1) for i in xic_df["i"].to_list()],
                    )
                )
                for lipid_label in sum_rank_score_info_dct:
                    rank_score_info_dct = sum_rank_score_info_dct[lipid_label]
                    msp_score_info_dct = hunt_msp_score(
                        pr_formula_candidates.get(lipid_label, {}),
                        pr_ms2_df,
                        {},
                        ms2_ppm,
                    )
                    fingerprint_score_info_dct = hunt_fingerprint_score(
                        pr_formula_candidates.get(lipid_label, {}),
                        pr_ms2_df,
                        {},
                        ms2_ppm,
                    )

                    snr_labels_dct = {
                        "signal_labels": [
                            rank_score_info_dct.get("signal_labels", None),
                            msp_score_info_dct.get("signal_labels", None),
                            fingerprint_score_info_dct.get("signal_labels", None),
                        ],
                        "noise_labels": [
                            rank_score_info_dct.get("noise_labels", None),
                            msp_score_info_dct.get("noise_labels", None),
                            fingerprint_score_info_dct.get("noise_labels", None),
                        ],
                    }
                    snr_score_info_dct = hunt_snr_score(snr_labels_dct)

                    total_score = get_total_score(
                        {
                            "isotope_score": isotope_score_info_dct.get(
                                "isotope_score", 0
                            ),
                            "rank_score": rank_score_info_dct.get("rank_score", 0),
                            "msp_score": msp_score_info_dct.get("msp_score", 0),
                            "fingerprint_score": fingerprint_score_info_dct.get(
                                "fingerprint_score", 0
                            ),
                            "snr_score": snr_score_info_dct.get("snr_score", 0),
                        }
                    )

                    result = Result(lipid_label)
                    result.add_xic_data(xic_data)
                    result.add_isotope_score_data(isotope_score_info_dct)
                    result.add_rank_score_data(rank_score_info_dct)
                    result.add_msp_score_data(msp_score_info_dct)
                    result.add_fingerprint_score_data(fingerprint_score_info_dct)
                    result.add_snr_score_data(snr_score_info_dct)
                    result.add_total_score_data(total_score)
                    result_dct = result.to_dict()
                    # print(result_dct)
                    result_title = result_dct.get("title", "")
                    if result_dct and result_title:
                        has_output_fig = plot_spectra(
                            result_dct, output_path=output_fig_abs_path, dpi=dpi
                        )
                        if has_output_fig:
                            ms2_result_dct[result_title] = result.to_dict()
    return ms2_result_dct


def hunt_star_ms2(
    r,
    neg_pr_df,
    pos_pr_df,
    neg_residue_space,
    pos_residue_space,
    selection_window,
    s_info_df,
    sum_spec_dct,
    ms1_ppm,
    ms1_abs_threshold,
    min_isotope_score,
    lipid_space,
    min_rank_score,
    sum_xic_df,
    ms2_ppm,
    output_fig_abs_path,
):
    ms2_result_dct = {}

    for chg_formula in sum_isotope_score_info_dct:
        isotope_score_info_dct = sum_isotope_score_info_dct[chg_formula]
        pr_ms2_df = isotope_score_info_dct.get("ms2_df", pd.DataFrame())
        pr_residue_space = isotope_score_info_dct.get("pr_residue_space")
        # print(isotope_score_info_dct)
        if (
            isotope_score_info_dct.get("isotope_score", 0) >= min_isotope_score
            and not pr_ms2_df.empty
        ):
            pr_formula_candidates = lipid_space.get(
                isotope_score_info_dct.get("charged_formula", ""), {}
            )
            sum_rank_score_info_dct = hunt_rank_score(
                pr_candidates=pr_formula_candidates,
                ms2_df=pr_ms2_df,
                res_space_info=pr_residue_space,
                weight_df=pd.DataFrame(),
                rank_score_threshold=min_rank_score,
            )
            if sum_rank_score_info_dct:
                lib_mz = isotope_score_info_dct.get("mz", 0)
                xic_df = sum_xic_df.query(
                    f"{(1 - ms1_ppm * 0.000001) * lib_mz} <= mz <= {(1 + ms1_ppm * 0.000001) * lib_mz}"
                ).copy()
                xic_df.drop(columns="mz", inplace=True)
                xic_data = list(
                    zip(
                        [round(rt, 4) for rt in xic_df["rt"].to_list()],
                        [round(i, 1) for i in xic_df["i"].to_list()],
                    )
                )
                for lipid_label in sum_rank_score_info_dct:
                    rank_score_info_dct = sum_rank_score_info_dct[lipid_label]
                    msp_score_info_dct = hunt_msp_score(
                        pr_formula_candidates.get(lipid_label, {}),
                        pr_ms2_df,
                        {},
                        ms2_ppm,
                    )
                    fingerprint_score_info_dct = hunt_fingerprint_score(
                        pr_formula_candidates.get(lipid_label, {}),
                        pr_ms2_df,
                        {},
                        ms2_ppm,
                    )

                    snr_labels_dct = {
                        "signal_labels": [
                            rank_score_info_dct.get("signal_labels", None),
                            msp_score_info_dct.get("signal_labels", None),
                            fingerprint_score_info_dct.get("signal_labels", None),
                        ],
                        "noise_labels": [
                            rank_score_info_dct.get("noise_labels", None),
                            msp_score_info_dct.get("noise_labels", None),
                            fingerprint_score_info_dct.get("noise_labels", None),
                        ],
                    }
                    snr_score_info_dct = hunt_snr_score(snr_labels_dct)

                    total_score = get_total_score(
                        {
                            "isotope_score": isotope_score_info_dct.get(
                                "isotope_score", 0
                            ),
                            "rank_score": rank_score_info_dct.get("rank_score", 0),
                            "msp_score": msp_score_info_dct.get("msp_score", 0),
                            "fingerprint_score": fingerprint_score_info_dct.get(
                                "fingerprint_score", 0
                            ),
                            "snr_score": snr_score_info_dct.get("snr_score", 0),
                        }
                    )

                    result = Result(lipid_label)
                    result.add_xic_data(xic_data)
                    result.add_isotope_score_data(isotope_score_info_dct)
                    result.add_rank_score_data(rank_score_info_dct)
                    result.add_msp_score_data(msp_score_info_dct)
                    result.add_fingerprint_score_data(fingerprint_score_info_dct)
                    result.add_snr_score_data(snr_score_info_dct)
                    result.add_total_score_data(total_score)
                    result_dct = result.to_dict()
                    # print(result_dct)
                    result_title = result_dct.get("title", "")
                    if result_dct and result_title:
                        has_output_fig = plot_spectra(
                            result_dct, output_path=output_fig_abs_path
                        )
                        if has_output_fig:
                            ms2_result_dct[result_title] = result.to_dict()
    return ms2_result_dct


def hunt_jobs(
    jobs,
    neg_pr_df,
    pos_pr_df,
    neg_residue_space,
    pos_residue_space,
    selection_window,
    s_info_df,
    sum_spec_dct,
    ms1_ppm,
    ms1_abs_threshold,
    min_isotope_score,
    lipid_space,
    min_rank_score,
    sum_xic_df,
    ms2_ppm,
    output_fig_abs_path,
    dpi=300,
):
    job_results_lst = []
    for job in jobs:
        job_result = hunt_mzml_ms2(
            job,
            neg_pr_df,
            pos_pr_df,
            neg_residue_space,
            pos_residue_space,
            selection_window,
            s_info_df,
            sum_spec_dct,
            ms1_ppm,
            ms1_abs_threshold,
            min_isotope_score,
            lipid_space,
            min_rank_score,
            sum_xic_df,
            ms2_ppm,
            output_fig_abs_path,
            dpi=dpi,
        )
        job_results_lst.append(job_result)
    job_results_dct = dict(ChainMap(*job_results_lst))
    return job_results_dct


def hunt(
    spectra_file_path: str,
    output_file_path: str,
    output_folder_path: str,
    search_space: dict,
    tiger_folder: str,
    spectra_type: str,
    mz_range: list = None,
    rt_range: list = None,
    selection_window: float = 0.5,
    ms1_ppm: int = 20,
    ms2_ppm: int = 100,
    ms1_abs_threshold: int = 1000,
    ms2_percent_threshold: float = 0.01,
    min_total_score: int = 50,
    min_isotope_score: int = 80,
    min_rank_score: int = 50,
    min_msp_score: int = 0,
    min_fp_score: int = 0,
    min_snr_score: int = 0,
    save_json: bool = True,
    dpi=300,
    worker: int = 3,
):
    t0 = time.time()
    time_tag = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    has_results = False
    if mz_range is None:
        mz_min = 200
        mz_max = 1200
    elif len(mz_range) == 2:
        mz_min = mz_range[0]
        mz_max = mz_range[1]
    else:
        mz_min = 200
        mz_max = 1200
    if rt_range is None:
        rt_min = 1.0
        rt_max = 30.0
    elif len(rt_range) == 2:
        rt_min = rt_range[0]
        rt_max = rt_range[1]
    else:
        rt_min = 1.0
        rt_max = 30.0

    rt_min -= 0.25  # to avoid the 1st peak is MS2
    if rt_min < 0:
        rt_min = 0
    search_space_size = 0
    if isinstance(search_space, str):
        if search_space.lower().endswith(".json"):
            try:
                if os.path.isfile(search_space):
                    with open(search_space, mode="r", encoding="utf-8") as json_fp:
                        search_space = json.load(json_fp)
                        search_space_size = len(list(search_space.keys()))
                        print(f"Search_space_size: {search_space_size}")
                else:
                    return time.time() - t0, has_results
            except Exception as e:
                print(e)
                return time.time() - t0, has_results
        else:
            return time.time() - t0, has_results
    elif isinstance(search_space, dict):
        search_space_size = len(list(search_space.keys()))
        print(f"Search_space_size: {search_space_size}")
    else:
        return time.time() - t0, has_results

    s_info_df, sum_spec_dct, sum_xic_df = None, None, None
    is_valid_spectra = False
    is_star_csv = False
    is_mzml = False

    run_params_dct = {
        "input_epilipidome": search_space,
        "input_spectra": spectra_file_path,
        "report_folder_path": output_folder_path,
        "export_table_path": output_file_path,
        "ms1_ppm": ms1_ppm,
        "ms2_ppm": ms2_ppm,
        "ms1_abs_threshold": ms1_abs_threshold,
        "ms2_relative_threshold": ms2_percent_threshold,
        "rt_min": rt_min,
        "rt_max": rt_max,
        "mz_min": mz_min,
        "mz_max": mz_max,
        "selection_window": selection_window,
        "total_score_min": min_total_score,
        "isotope_score_min": min_isotope_score,
        "rank_score_min": min_rank_score,
        "msp_score_min": min_msp_score,
        "fingerprint_score_min": min_fp_score,
        "snr_score_min": min_snr_score,
        "save_json": save_json,
    }

    results_dct_lst = []
    results_dct = {}
    output_fig_path = os.path.join(
        output_folder_path, f"LPPtiger_Results_Figures_{time_tag}"
    )
    os.mkdir(output_fig_path)
    output_fig_abs_path = os.path.abspath(output_fig_path)

    if os.path.exists(spectra_file_path):
        if (
            spectra_file_path.lower().endswith(".mzml")
            and spectra_type.lower() == "mzml"
        ):
            (
                s_info_df,
                sum_spec_dct,
                sum_xic_df,
            ) = mzml_spectra_reader.extract_mzml_spectra(
                spectra_file_path, rt_range=[rt_min, rt_max]
            )
            is_mzml = True
        elif spectra_file_path.lower().endswith(".csv") and spectra_type.lower() in [
            "star",
            "star2",
            "lipostar",
            "lipostar2",
        ]:
            (
                s_info_df,
                sum_spec_dct,
                sum_xic_df,
            ) = star_spectra_reader.extract_star_spectra(
                spectra_file_path, rt_range=[rt_min, rt_max]
            )
            is_star_csv = True
        else:
            pass
    else:
        pass

    if (
        isinstance(s_info_df, pd.DataFrame)
        and isinstance(sum_spec_dct, dict)
        and isinstance(sum_xic_df, pd.DataFrame)
    ):
        if not s_info_df.empty and not sum_xic_df.empty:
            is_valid_spectra = True

    if is_valid_spectra:
        pre_s_ms2_info_df = s_info_df[s_info_df["DDA_rank"] >= 0].copy()
        mz_s_ms2_info_df = pre_s_ms2_info_df.query(
            f"{mz_min} <= MS2_PR_mz <= {mz_max}"
        ).copy()
        s_ms2_info_df = mz_s_ms2_info_df.query(
            f"{rt_min} <= scan_time <= {rt_max}"
        ).copy()
        spectra_count = s_ms2_info_df.shape[0]
        print(
            f"Spectra loaded. Total features: {spectra_count} in m/z {mz_min}~{mz_max} @ RT {rt_min}~{rt_max} "
        )
        if "lipid_space" in search_space and "residue_space" in search_space:
            pass
        else:
            search_space = build_search_space(
                search_space,
                charge_modes=["neg", "pos"],
                adducts=["[M-H]-", "[M+HCOO]-", "[M+Na]+"],
                mz_range=mz_range,
            )
        lipid_space = search_space.get("lipid_space", {})
        search_ranges = build_search_ranges(
            search_space, ms1_ppm=ms1_ppm, ms2_ppm=ms2_ppm
        )

        pos_pr_df = search_ranges.get("pos", {}).get("pr", {}).get("space")
        pos_class_df = search_ranges.get("pos", {}).get("class", {}).get("space")
        pos_frag_df = search_ranges.get("pos", {}).get("frag", {}).get("space")
        pos_nl_df = search_ranges.get("pos", {}).get("nl", {}).get("space")
        neg_pr_df = search_ranges.get("neg", {}).get("pr", {}).get("space")
        neg_class_df = search_ranges.get("neg", {}).get("class", {}).get("space")
        neg_frag_df = search_ranges.get("neg", {}).get("frag", {}).get("space")
        neg_nl_df = search_ranges.get("neg", {}).get("nl", {}).get("space")

        pos_residue_space = {
            "class_space": pos_class_df,
            "frag_space": pos_frag_df,
            "nl_space": pos_nl_df,
            "frag_mz_range": search_ranges.get("pos", {})
            .get("frag", {})
            .get("mz_range"),
            "nl_mz_range": search_ranges.get("pos", {}).get("nl", {}).get("mz_range"),
        }
        neg_residue_space = {
            "class_space": neg_class_df,
            "frag_space": neg_frag_df,
            "nl_space": neg_nl_df,
            "class_mz_range": search_ranges.get("neg", {})
            .get("class", {})
            .get("mz_range"),
            "frag_mz_range": search_ranges.get("neg", {})
            .get("frag", {})
            .get("mz_range"),
            "nl_mz_range": search_ranges.get("neg", {}).get("nl", {}).get("mz_range"),
        }
        # if isinstance(pos_pr_df, pd.DataFrame):
        #     pos_pr_formula_lst = pos_pr_df["formula"].tolist()
        # else:
        #     pos_pr_formula_lst = []
        # if isinstance(neg_pr_df, pd.DataFrame):
        #     neg_pr_formula_lst = neg_pr_df["formula"].tolist()
        # else:
        #     neg_pr_formula_lst = []
        #
        # matched_idx_lst = []
        # unmatched_idx_lst = []
        # ident_results = {}

        print(f"Hunting for lipids ...")

        if is_mzml:
            print(f"Mode Tiger mzml ...")
            s_ms2_info_lst = s_ms2_info_df.to_dict(orient="records")
            ms2_spectra_count = len(s_ms2_info_lst)
            if ms2_spectra_count > 1:
                pass
            else:
                ms2_spectra_count = 1
            if isinstance(worker, int):
                if worker > ms2_spectra_count:
                    worker = ms2_spectra_count
                else:
                    worker = min(ms2_spectra_count, worker)
                    if worker > 16:
                        worker = 16
            else:
                worker = min(ms2_spectra_count, 16)
            if 2000000 <= search_space_size:
                worker = 1
            elif 1000000 <= search_space_size < 2000000:
                worker = 3
            elif 500000 <= search_space_size < 1000000:
                worker = 3
            elif 100000 <= search_space_size < 500000:
                worker = 4
            elif 50000 <= search_space_size < 100000:
                worker = 6
            elif 10000 <= search_space_size < 50000:
                worker = 8
            if worker > 1:
                pool = Pool(processes=worker)
                job_steps = math.ceil(ms2_spectra_count / worker)
                job_split_idx_lst = list(range(0, ms2_spectra_count, job_steps))[:-1]
                # job_split_idx_lst.append(ms2_spectra_count-1)
                print("split jobs at index", [job_split + job_steps for job_split in job_split_idx_lst])
                job_lst = []
                for idx in range(1, len(job_split_idx_lst)):
                    job_lst.append(
                        s_ms2_info_lst[
                            job_split_idx_lst[idx - 1] : job_split_idx_lst[idx]
                        ]
                    )
                job_lst.append(
                    s_ms2_info_lst[job_split_idx_lst[-1] :]
                )  # add the last job segment
                print(
                    f"Split #{ms2_spectra_count} tasks to {worker} workers each has {job_steps} tasks"
                )
                print("start multiprocessing...")
                result_lst = []
                for job in job_lst:
                    r = pool.apply_async(
                        hunt_jobs,
                        args=(
                            job,
                            neg_pr_df,
                            pos_pr_df,
                            neg_residue_space,
                            pos_residue_space,
                            selection_window,
                            s_info_df,
                            sum_spec_dct,
                            ms1_ppm,
                            ms1_abs_threshold,
                            min_isotope_score,
                            lipid_space,
                            min_rank_score,
                            sum_xic_df,
                            ms2_ppm,
                            output_fig_abs_path,
                        ),
                    )  # dict
                    result_lst.append(r)
                pool.close()
                pool.join()
                print("multiprocessing finished...")
                print("start merging results...")
                for r_dct in result_lst:
                    if isinstance(r_dct, dict):
                        worker_result_dct = r_dct
                    else:
                        # worker_result_dct = r_dct.get()  # debug only
                        try:
                            worker_result_dct = r_dct.get()
                        except (KeyError, SystemError, ValueError, TypeError):
                            print(f"Cannot load results from multiprocessing...")
                            worker_result_dct = {}
                    results_dct_lst.append(worker_result_dct)
                results_dct = dict(ChainMap(*results_dct_lst))
                print("results merged...")
            else:
                print("start processing...")
                results_dct = hunt_jobs(
                    s_ms2_info_lst,
                    neg_pr_df,
                    pos_pr_df,
                    neg_residue_space,
                    pos_residue_space,
                    selection_window,
                    s_info_df,
                    sum_spec_dct,
                    ms1_ppm,
                    ms1_abs_threshold,
                    min_isotope_score,
                    lipid_space,
                    min_rank_score,
                    sum_xic_df,
                    ms2_ppm,
                    output_fig_abs_path,
                    dpi=dpi
                )
                print("processing finished...")
            if is_star_csv:
                print(f"Mode Tiger star csv ...")
                results_df_lst = []
                s_ms2_info_df["Custom Adduct"] = ""
                s_ms2_info_df["Custom Formula"] = ""
                s_ms2_info_df["LPPtiger Score"] = 0
                s_ms2_info_df["rank_score"] = 0
                s_ms2_info_df["msp_score"] = 0
                s_ms2_info_df["fingerprint_score"] = 0
                s_ms2_info_df["snr_score"] = 0
                for idx, r in s_ms2_info_df.iterrows():
                    chg_pr_formula = r["MS2_PR_Formula"]
                    pr_charge_mode = chg_pr_formula[-1]
                    if pr_charge_mode == "-":
                        pr_residue_space = neg_residue_space
                    elif pr_charge_mode == "+":
                        pr_residue_space = pos_residue_space
                    else:
                        pr_residue_space = neg_residue_space
                    pr_mz = r["MS2_PR_mz"]
                    scan_number = r["scan_number"]
                    scan_time = r["scan_time"]
                    star_pr_isotope_score = r["isotope_score"]
                    feature_label = f"{pr_mz}@{scan_time}"
                    r_spec_index = r["spec_index"]

                    pr_formula_candidates = lipid_space.get(chg_pr_formula, {})
                    pr_ms2_df = sum_spec_dct.get(r_spec_index)
                    rank_score_results = get_rank_score(
                        pr_candidates=pr_formula_candidates,
                        ms2_df=pr_ms2_df,
                        res_space_info=pr_residue_space,
                        weight_df=pd.DataFrame(),
                        rank_score_threshold=min_rank_score,
                    )

                    if rank_score_results:
                        print(f"#{scan_number}# {feature_label}:")

                        for rs_key in rank_score_results:
                            rs_score_info = get_isotope_score(
                                rank_score_results[rs_key], star_pr_isotope_score
                            )
                            rs_score_info = get_msp_score(
                                pr_formula_candidates.get(rs_key, {}),
                                pr_ms2_df,
                                rs_score_info,
                                ms2_ppm=ms2_ppm,
                            )
                            rs_score_info = get_fingerprint_score(
                                pr_formula_candidates.get(rs_key, {}),
                                pr_ms2_df,
                                rs_score_info,
                                ms2_ppm=ms2_ppm,
                            )
                            rs_score_info = get_snr_score(
                                pr_formula_candidates.get(rs_key, {}),
                                pr_ms2_df,
                                rs_score_info,
                            )
                            rs_score_info_lst = [
                                rs_score_info.get("rank_score", 0),
                                rs_score_info.get("isotope_score", 0),
                                rs_score_info.get("msp_score", 0),
                                rs_score_info.get("fingerprint_score", 0),
                                rs_score_info.get("snr_score", 0),
                            ]
                            tiger_sum_score = round(mean(rs_score_info_lst), 1)
                            if tiger_sum_score >= min_tiger_sum_score:
                                signal_peaks_dct = rs_score_info.get("signal_peaks", {})
                                print(
                                    f"{rs_key} # sum score {tiger_sum_score} # "
                                    f"rnk:{rs_score_info.get('rank_score')} | "
                                    f"ist:{rs_score_info.get('isotope_score')} | "
                                    f"msp:{rs_score_info.get('msp_score')} | "
                                    f"fgp:{rs_score_info.get('fingerprint_score')} | "
                                    f"snr:{rs_score_info.get('snr_score')}"
                                )
                                r["Custom Name"] = rs_key
                                r["name"] = rank_score_results[rs_key].get("name", 0)
                                r["adduct"] = rank_score_results[rs_key].get(
                                    "adduct", ""
                                )
                                r["neutral_formula"] = rank_score_results[rs_key].get(
                                    "neutral_formula", ""
                                )
                                r["LPPtiger_score"] = tiger_sum_score
                                r["rank_score"] = rs_score_info.get("rank_score", 0)
                                r["isotope_score"] = rs_score_info.get(
                                    "isotope_score", 0
                                )
                                r["msp_score"] = rs_score_info.get("msp_score", 0)
                                r["fingerprint_score"] = rs_score_info.get(
                                    "fingerprint_score", 0
                                )
                                r["snr_score"] = rs_score_info.get("snr_score", 0)
                                r["matched_peaks"] = len(signal_peaks_dct.get("mz", []))
                                # r["matched_residues"] = json.dumps(rs_score_info.get("obs_residues", []))
                                try:
                                    signal_peaks_tp_lst = list(
                                        zip(
                                            [
                                                str(s_mz)
                                                for s_mz in signal_peaks_dct.get(
                                                    "mz", []
                                                )
                                            ],
                                            [
                                                str(s_i)
                                                for s_i in signal_peaks_dct.get("i", [])
                                            ],
                                        )
                                    )
                                    signal_peaks_lst = [
                                        ";".join(p) for p in signal_peaks_tp_lst
                                    ]
                                    signal_peaks_str = (
                                        "[" + ",".join(signal_peaks_lst) + "]"
                                    )
                                except Exception as e:
                                    print(e)
                                    signal_peaks_str = []
                                r["identified_peak_list"] = signal_peaks_str
                                results_df_lst.append(r.rename(rs_key))

                            # if chg_pr_formula in pr_formula_list:
                    #     pr_candidates = lipid_space[pr_charged_formula]
                    #     matched_idx_lst.append(idx)
                    #     pr_ms1_spec_df = sum_spec_dct.get(scan_number - 1)
                    #     pr_ms2_spec_df = sum_spec_dct.get(scan_number)
                    #     if isinstance(pr_ms2_spec_df, pd.DataFrame):
                    #         rank_score_results = get_rank_score(
                    #             pr_candidates=pr_formula_candidates,
                    #             ms2_df=sum_spec_dct.get(r_spec_index),
                    #             res_space_info=pr_residue_space,
                    #             weight_df=pd.DataFrame(),
                    #             rank_score_threshold=min_rank_score
                    #         )
                    #         rank_score_results = get_rank_score(
                    #             pr_candidates,
                    #             pr_ms2_spec_df,
                    #             res_space_info,
                    #             pd.DataFrame(),
                    #         )
                    #
                    #         ident_remarks = f"{feature_label} | {pr_charged_formula}"
                    #         feature_results = {}
                    #         if rank_score_results:
                    #             for r_candidate_label in rank_score_results:
                    #                 r_score_dct = rank_score_results[r_candidate_label]
                    #                 r_score_dct = get_isotope_score(
                    #                     r_score_dct, star_pr_isotope_score
                    #                 )
                    #                 r_score_dct = get_msp_score(
                    #                     pr_candidates[r_candidate_label],
                    #                     pr_ms2_spec_df,
                    #                     r_score_dct,
                    #                     ppm=ppm,
                    #                 )
                    #                 r_score_dct = get_fingerprint_score(
                    #                     pr_candidates[r_candidate_label],
                    #                     pr_ms2_spec_df,
                    #                     r_score_dct,
                    #                     ppm=ppm,
                    #                 )
                    #                 r_score_dct = get_snr_score(
                    #                     pr_candidates[r_candidate_label],
                    #                     pr_ms2_spec_df,
                    #                     r_score_dct,
                    #                 )
                    #                 r_total_score = get_total_score(r_score_dct)
                    #                 r_score_dct["total_score"] = r_total_score
                    #                 # r_score_dct["label"] = r_candidate_label
                    #                 r_score_dct[
                    #                     "title"
                    #                 ] = f"{ident_remarks} | {r_candidate_label} | {r_total_score}"
                    #                 feature_results[
                    #                     f"{r_total_score} | {r_candidate_label}"
                    #                 ] = r_score_dct
                    #         if feature_results:
                    #             feature_result_keys = list(
                    #                 natsorted(list(feature_results.keys()), reverse=True)
                    #             )
                    #             candidates = {}
                    #             info = {}
                    #             for candidate_key in feature_result_keys:
                    #                 candidate_label = feature_results[candidate_key].get(
                    #                     "label", ""
                    #                 )
                    #                 candidates[candidate_label] = feature_results[
                    #                     candidate_key
                    #                 ].get("total_score", 0.0)
                    #                 info[candidate_label] = feature_results[candidate_key]
                    #             ident_results[feature_label] = {
                    #                 "pr_mz": pr_mz,
                    #                 "rt": scan_time,
                    #                 "ms1": pr_ms1_spec_df[["mz", "i"]].to_dict(
                    #                     orient="list"
                    #                 ),
                    #                 "ms2": pr_ms2_spec_df[["mz", "i"]].to_dict(
                    #                     orient="list"
                    #                 ),
                    #                 "candidates": candidates,
                    #                 "info": info,
                    #             }
                    #     else:
                    #         unmatched_idx_lst.append(idx)
                    # else:
                    #     unmatched_idx_lst.append(idx)

                    # progress_counter += 1
                merged_results_df = pd.concat(results_df_lst, axis=1).T
                if output_file_path.lower().endswith(".csv"):
                    pass
                else:
                    output_file_path += ".csv"
                is_exported = export_to_lipostar(
                    filepath=output_file_path, df=merged_results_df
                )
                if is_exported:
                    print(f"output generated: {output_file_path}")

    has_output_file, final_output_dct = export_to_xlsx(results_dct, output_file_path)

    if has_output_file and results_dct and save_json:
        output_json_path = os.path.join(
            output_fig_abs_path, f"LPPtiger_{time_tag}.json"
        )
        with open(output_json_path, "w") as json_fp:
            json.dump(
                {
                    "parameters": run_params_dct,
                    "results": results_dct,
                    "table": final_output_dct,
                },
                json_fp,
            )
            print(f"Dump results as JSON: {output_json_path}")

    if has_output_file and final_output_dct:
        params_dct = {
            "mz_start": mz_min,
            "mz_end": mz_max,
            "rt_start": rt_min,
            "rt_end": rt_max,
            "ms_th": ms1_abs_threshold,
            "ms2_th": ms2_percent_threshold,
            "ms_ppm": ms1_ppm,
            "ms2_ppm": ms2_ppm,
            "score_filter": min_total_score,
            "isotope_score_filter": min_isotope_score,
            "tiger_folder": tiger_folder,
        }
        log_page_creator = LogPageCreator(
            params_dct, final_output_dct, output_folder_path, time_tag
        )
        log_page_creator.add_all_info()
        log_page_creator.close_page()
        has_results = True

    delta_t = time.time() - t0
    return delta_t, has_results


if __name__ == "__main__":

    test_params = {
        # "spectra": r"./test/mzml/QE_21_24_MF61.mzML",
        # "spectra": r"./test/mzml/QE_18_42_Batch1_RS_Leipzig_72h_neg_rep2.mzML",
        # "spectra": r"./test/mzml/oxPC_neg.mzML",
        "spectra": r"test/mzml/DR1 20 min POS.mzML",
        # "spectra": r"test/mzml/4. oxTG 6h.mzML",
        # "spectra": r"test/mzml/DS1.mzML",
        # "spectra": r"./test/mzml/Ex_22_16_MW24_oxCE_6h.mzML",
        # "spectra": r"./test/mzml/070120_CM_neg_70min_SIN_II.mzML",
        # "spectra": r"./test/mzml/IR_18 min pol switch.mzML",
        # "spectra": r"./test/mzml/IR_18 min neg.mzML",
        # "epilipidome": r"./test/output/oxPC_full.json",
        # "epilipidome": r"./test/output/oxPL_TG_CE.json",
        # "epilipidome": r"./test/output/oxPC.json",
        # "epilipidome": r"./test/output/oxPL_TG_CE_mega.json",
        # "epilipidome": r"test/output/oxPC_PE_TG_CE.json",
        "epilipidome": r"test/output/oxTG_lite.json",
        # "epilipidome": r"test/output/oxCE.json",
        # "epilipidome": r"test/output/oxPL_TG_CE.json",
        # "mz_range": [400, 1200],
        # "mz_range": [819, 825],
        "mz_range": [895, 896],
        # "rt_range": [3, 30],
        # "rt_range": [11.1, 11.5],
        # "rt_range": [5, 18],
        "rt_range": [23.5, 24],
        "selection_window": 0.5,
        "ms1_ppm": 10,
        "ms2_ppm": 50,
        "ms1_abs_threshold": 2000,
        "ms2_percent_threshold": 0.001,
        "min_total_score": 30,
        "min_isotope_score": 80,
        "min_rank_score": 30,
        # "output_file": r"./test/output/fig/oxPC_oxPE_ident.xlsx",
        # "output_file": r"test/output/fig/oxTG_ident.xlsx",
        # "output_file": r"test/output/fig/oxPL_ident.xlsx",
        "output_file": r"test/output/GUI/DR1_oxTG.xlsx",
        "output_folder": r"test/output/GUI",
        "worker_count": 1,
        "spectra_type": "mzml",
    }

    tiger_cwd_str = os.getcwd()
    # os.chdir(tiger_cwd_str)
    # os.chdir("../../")
    print(f"Tiger working path: {os.getcwd()}")

    # with open(test_params["epilipidome"], mode="r", encoding="utf-8") as usr_fp:
    #     test_search_space = json.load(usr_fp)

    hunt_time, has_identification = hunt(
        test_params["spectra"],
        test_params["output_file"],
        test_params["output_folder"],
        test_params["epilipidome"],
        # test_search_space,
        tiger_cwd_str,
        mz_range=test_params["mz_range"],
        rt_range=test_params["rt_range"],
        min_total_score=test_params["min_total_score"],
        min_rank_score=test_params["min_rank_score"],
        ms1_ppm=test_params["ms1_ppm"],
        ms2_ppm=test_params["ms2_ppm"],
        worker=test_params["worker_count"],
        spectra_type=test_params["spectra_type"],
    )
    print(f"Runtime: {hunt_time:.2f} seconds")
    print("FIN")
