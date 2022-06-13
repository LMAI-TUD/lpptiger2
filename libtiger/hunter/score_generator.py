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
import math

from natsort import natsorted
import numpy as np
import pandas as pd
from scipy import spatial

from libtiger.hunter.isotope_score import IsotopeHunter
from libtiger.reader import mzml_spectra_reader


def get_star_isotope_score(score_summary: dict, isotope_score: float = 75.0):
    score_summary["isotope_score"] = round(isotope_score, 6)

    return score_summary


def get_rank_obs_peaks(ms2_df: pd.DataFrame, res_space_info: dict):
    rank_obs_dct = {}
    for sub_res in ["class", "frag", "nl"]:
        residue_df = res_space_info.get(f"{sub_res}_space", pd.DataFrame())
        if isinstance(residue_df, pd.DataFrame) and not residue_df.empty:
            residue_range = res_space_info.get(f"{sub_res}_mz_range", {})
            res_mz_min = residue_range.get("min_mz", 0)
            res_mz_max = residue_range.get("max_mz", 0)
            if res_mz_min and res_mz_max:
                res_range_ms2_df = ms2_df[
                    (ms2_df["mz"] >= res_mz_min)
                    & (ms2_df["mz"] <= res_mz_max)
                ].copy()
            else:
                res_range_ms2_df = ms2_df.copy()
            res_range_ms2_df.sort_values(by="i", ascending=False, inplace=True)
            res_range_ms2_df.reset_index(drop=True, inplace=True)
            res_range_ms2_pkl = list(zip(res_range_ms2_df["mz"], res_range_ms2_df["i"]))
            sub_obs_dct = {}
            rank_counter = 0
            for mz_tp in res_range_ms2_pkl:
                peak_mz = mz_tp[0]
                peak_i = mz_tp[1]
                sub_obs_df = residue_df[
                    (residue_df["min_mz"] <= peak_mz) & (residue_df["max_mz"] >= peak_mz)
                ].copy()
                if not sub_obs_df.empty:
                    sub_obs_df["ppm"] = (
                        1000000 * (peak_mz - sub_obs_df["mz"]) / sub_obs_df["mz"]
                    )
                    sub_obs_df["abs_ppm"] = sub_obs_df["ppm"].abs()
                    sub_obs_df[["name", "adduct"]] = sub_obs_df["label"].str.split(
                        "#", 1, expand=True
                    )
                    if rank_counter < 10:
                        sub_obs_dct[rank_counter] = {
                            "rank": rank_counter,
                            "obs_mz": peak_mz,
                            "obs_i": peak_i,
                            "obs_df": sub_obs_df,
                        }
                        # if sub_obs_df.shape[0] == 1:
                        #     theo_mz = sub_obs_df["mz"].values.tolist()[0]
                        #     label = sub_obs_df["label"].values.tolist()[0]
                        #     ppm = 1000000 * (peak_mz - theo_mz) / theo_mz
                        #     abs_ppm = abs(ppm)
                        # else:
                        #     sub_obs_df.at[:, "ppm"] = (
                        #         1000000 * (peak_mz - sub_obs_df["mz"]) / sub_obs_df["mz"]
                        #     )
                        #     sub_obs_df.at[:, "abs_ppm"] = sub_obs_df["ppm"].abs()
                        #     sub_obs_df.sort_values(by="abs_ppm", ascending=True, inplace=True)
                        #     sub_obs_df.reset_index(drop=True, inplace=True)
                        #     theo_mz = sub_obs_df["mz"].values.tolist()[0]
                        #     label = sub_obs_df["label"].values.tolist()[0]
                        #     ppm = sub_obs_df["ppm"].values.tolist()[0]
                        #     abs_ppm = sub_obs_df["abs_ppm"].values.tolist()[0]
                        # if "#" in label:
                        #     lipid_name, adduct = label.split("#")
                        #     raw_sub_obs_dct[rank_counter] = {
                        #         "rank": rank_counter,
                        #         "mz": peak_mz,
                        #         "i": peak_i,
                        #         "theo_mz": theo_mz,
                        #         "ppm": ppm,
                        #         "abs_ppm": abs_ppm,
                        #         "label": label,
                        #         "name": lipid_name,
                        #         "adduct": adduct,
                        #     }
                        rank_counter += 1
                    else:
                        pass
            # if raw_sub_obs_dct:
            #     sub_obs_peaks_df = pd.DataFrame(data=raw_sub_obs_dct).T
            #     sub_obs_peaks_df.sort_values(
            #         by=["abs_ppm", "i"], ascending=[True, False], inplace=True
            #     )
            #     sub_obs_peaks_df.drop_duplicates(["label"], keep="first", inplace=True)
            #     sub_obs_peaks_df["rank"] = sub_obs_peaks_df["rank"].astype(int)
            #     sub_obs_peaks_df.set_index("rank", drop=True, inplace=True)
            #     sub_obs_dct = sub_obs_peaks_df.to_dict(orient="index")
            # else:
            #     sub_obs_dct = {}

            rank_obs_dct[sub_res] = sub_obs_dct

    return rank_obs_dct


def get_msp_vectors(msp_info: dict, mz_power: float = 3, i_power: float = 0.6):
    msp_vectors = []
    peak_pairs = list(zip(msp_info["mz"], msp_info["i"]))
    for pk_pair in peak_pairs:
        msp_vectors.append((pk_pair[0] ** mz_power) * (pk_pair[1] ** i_power))
    return msp_vectors


def hunt_isotope_score(
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
):
    sum_isotope_score_info_dct = {}
    isotope_hunter = IsotopeHunter()
    ms2_dda_event_idx = r["dda_event_idx"]
    ms2_spec_index = r["spec_index"]
    ms2_pr_mz = r["MS2_PR_mz"]
    ms2_polarity = r["polarity"]

    # Due to the selection window for MS2 (usually 1 or 1.5 --> +/- 0.5 or +/- 0.75)
    # Precursor m/z in mzML MS2 can have large differences than actual MS1 precursor m/z
    if ms2_polarity == "neg":
        pr_candidates_df = neg_pr_df[
            (neg_pr_df["min_mz"] >= ms2_pr_mz - selection_window)
            & (neg_pr_df["min_mz"] <= ms2_pr_mz + selection_window)
            & (neg_pr_df["max_mz"] >= ms2_pr_mz - selection_window)
            & (neg_pr_df["max_mz"] <= ms2_pr_mz + selection_window)
        ].copy()
        pr_residue_space = neg_residue_space

    elif ms2_polarity == "pos":
        pr_candidates_df = pos_pr_df[
            (pos_pr_df["min_mz"] <= ms2_pr_mz + selection_window)
            & (pos_pr_df["max_mz"] >= ms2_pr_mz - selection_window)
        ].copy()
        pr_residue_space = pos_residue_space
    else:
        pr_candidates_df = pd.DataFrame()
        pr_residue_space = {}

    if not pr_candidates_df.empty:
        pr_ms1_df = s_info_df[
            (s_info_df["dda_event_idx"] == ms2_dda_event_idx)
            & (s_info_df["DDA_rank"] == 0)
            & (s_info_df["spec_index"] <= ms2_spec_index)
            & (s_info_df["polarity"] == ms2_polarity)
        ].copy()
        if not pr_ms1_df.empty:
            pr_ms1_spec_idx = sorted(pr_ms1_df["spec_index"].tolist())[0]
            pr_ms1_spec_df = sum_spec_dct.get(pr_ms1_spec_idx)
            pr_ms1_spec_df.sort_values(by="mz", ascending=True, inplace=True)
            unique_formula_lst = list(
                set(
                    zip(
                        pr_candidates_df["formula"].tolist(),
                        pr_candidates_df["mz"].tolist(),
                    )
                )
            )
            obs_ms1_pr_lst = []
            for chg_pr_tp in unique_formula_lst:
                chg_pr_formula = chg_pr_tp[0]
                (
                    ms1_obs_pr_mz,
                    ms1_obs_pr_i,
                    ms1_obs_pr_ppm,
                ) = mzml_spectra_reader.find_ms1_pr(
                    ms1_df=pr_ms1_spec_df,
                    mz_lib=chg_pr_tp[1],
                    ms1_ppm=ms1_ppm,
                )
                obs_ms1_pr_lst.append(
                    [
                        chg_pr_formula,
                        chg_pr_tp[1],
                        ms1_obs_pr_mz,
                        ms1_obs_pr_i,
                        ms1_obs_pr_ppm,
                        abs(ms1_obs_pr_ppm),
                    ]
                )
            raw_obs_ms1_pr_df = pd.DataFrame(
                obs_ms1_pr_lst,
                columns=["chg_pr_formula", "lib_mz", "obs_mz", "i", "ppm", "abs_ppm"],
            )
            obs_ms1_pr_df = raw_obs_ms1_pr_df[
                raw_obs_ms1_pr_df["abs_ppm"] <= ms1_ppm
            ].copy()
            obs_ms1_pr_df.sort_values(
                by=["i", "abs_ppm"], ascending=[False, True], inplace=True
            )
            obs_ms1_pr_df.reset_index(drop=True, inplace=True)
            obs_ms1_pr_df.drop_duplicates(keep="first", inplace=True)
            top_obs_ms1_pr_df = obs_ms1_pr_df.head(3)
            if not top_obs_ms1_pr_df.empty:
                top_obs_ms1_pr_dct = top_obs_ms1_pr_df.to_dict(orient="index")
            else:
                top_obs_ms1_pr_dct = {}
            for obs_ms1_pr_idx in top_obs_ms1_pr_dct:
                chg_pr_formula = top_obs_ms1_pr_dct[obs_ms1_pr_idx].get(
                    "chg_pr_formula"
                )
                ms1_obs_pr_mz = top_obs_ms1_pr_dct[obs_ms1_pr_idx].get("obs_mz", 0)
                isotope_score_info_dct = isotope_hunter.get_isotope_score(
                    ms1_obs_pr_mz,
                    top_obs_ms1_pr_dct[obs_ms1_pr_idx].get("i"),
                    chg_pr_formula,
                    pr_ms1_spec_df,
                    isotope_number=2,
                )
                isotope_score = isotope_score_info_dct.get("isotope_score", 0)
                observed_i = top_obs_ms1_pr_dct[obs_ms1_pr_idx].get("i", 0)
                observed_ppm = top_obs_ms1_pr_dct[obs_ms1_pr_idx].get("ppm", 999)
                print(isotope_score, observed_ppm, observed_i)
                if (
                    observed_i >= ms1_abs_threshold
                    and isotope_score >= min_isotope_score
                    and abs(observed_ppm) <= ms1_ppm
                ):

                    isotope_score_info_dct["charged_formula"] = chg_pr_formula
                    isotope_score_info_dct["pr_residue_space"] = pr_residue_space
                    isotope_score_info_dct["mz"] = top_obs_ms1_pr_dct[
                        obs_ms1_pr_idx
                    ].get("lib_mz", 0)
                    isotope_score_info_dct["observed_mz"] = ms1_obs_pr_mz
                    isotope_score_info_dct["observed_i"] = observed_i
                    isotope_score_info_dct["observed_ppm"] = observed_ppm
                    isotope_score_info_dct["ms1_scan_id"] = pr_ms1_spec_idx
                    isotope_score_info_dct["ms1_scan_time"] = pr_ms1_spec_df["rt"].max()
                    isotope_score_info_dct["ms2_polarity"] = pr_ms1_spec_df["rt"].max()
                    isotope_score_info_dct["ms1_df"] = pr_ms1_spec_df.drop(columns="rt")
                    isotope_score_info_dct["dda_event_id"] = ms2_dda_event_idx
                    isotope_score_info_dct["ms2_scan_id"] = ms2_spec_index
                    isotope_score_info_dct["ms2_scan_time"] = r["scan_time"]
                    isotope_score_info_dct["ms2_polarity"] = ms2_polarity
                    isotope_score_info_dct["ms2_precursor_mz"] = ms2_pr_mz
                    isotope_score_info_dct["ms2_df"] = sum_spec_dct.get(ms2_spec_index)

                    print(
                        f"isotope score check passed: {isotope_score} for {chg_pr_formula} @ {ms1_obs_pr_mz}"
                    )
                    sum_isotope_score_info_dct[
                        f"{chg_pr_formula}@{ms2_spec_index}"
                    ] = isotope_score_info_dct

    return sum_isotope_score_info_dct


def hunt_rank_score(
    pr_candidates: dict,
    ms2_df: pd.DataFrame,
    res_space_info: dict,
    weight_df: pd.DataFrame,
    rank_score_threshold: int = 30,
):
    rank_score_results = {}
    rank_obs_dct = get_rank_obs_peaks(
        ms2_df,
        res_space_info,
    )

    for pr_candidate in pr_candidates:
        pr_candidate_info = pr_candidates[pr_candidate]
        pr_name = pr_candidate_info.get("name", "")
        pr_info = pr_candidate_info.get("info", {})
        pr_rank_score = 0

        signal_labels = {}
        noise_labels = {}
        observed_fragments = {}
        observed_neutral_losses = {}
        observed_specific_peaks = {}
        observed_residues = []
        obs_matches_df = pd.DataFrame()
        res_name_lst = []
        for res_idx in ["FA1", "FA2", "FA3"]:
            res_name = pr_info.get(res_idx, {}).get("name", None)
            if res_name:
                res_name_lst.append(res_name)
        unique_res = len(list(set(res_name_lst)))
        frag_label_dct = pr_candidate_info.get("ions", {}).get("frag_info", {})
        nl_label_dct = pr_candidate_info.get("ions", {}).get("nl_info", {})
        class_label_dct = pr_candidate_info.get("ions", {}).get("class_info", {})
        nl_name_lst = list(set([nl_lb.split("#")[0] for nl_lb in nl_label_dct]))

        for sub_res in ["frag", "nl"]:
            if sub_res in rank_obs_dct:
                sub_rank_obs_dct = rank_obs_dct.get(sub_res, {})
                for rank in sub_rank_obs_dct:
                    is_count_for_score = False

                    if rank < 10:
                        obs_mz = sub_rank_obs_dct[rank].get("obs_mz", 0)
                        obs_i = sub_rank_obs_dct[rank].get("obs_i", 0)
                        obs_matches_df = sub_rank_obs_dct[rank].get("obs_df", pd.DataFrame)
                        obs_peak = [obs_mz, obs_i]
                        # obs_label = sub_rank_obs_dct[rank].get("label")
                        # obs_name = sub_rank_obs_dct[rank].get("name")
                        # obs_ppm = round(sub_rank_obs_dct[rank].get("ppm"), 1)
                        # obs_peak = [
                        #     sub_rank_obs_dct[rank].get("mz"),
                        #     sub_rank_obs_dct[rank].get("i"),
                        # ]

                        if obs_matches_df.shape[0] > 1:
                            if sub_res == "frag":
                                fit_obs_res_df = obs_matches_df[
                                    obs_matches_df["name"].isin(res_name_lst)
                                ].copy()
                            elif sub_res == "nl":
                                fit_obs_res_df = obs_matches_df[
                                    obs_matches_df["name"].isin(nl_name_lst)
                                ].copy()
                            else:
                                fit_obs_res_df = pd.DataFrame()
                            if fit_obs_res_df.shape[0] == 1:
                                obs_res_df = fit_obs_res_df
                            elif fit_obs_res_df.shape[0] > 1:
                                fit_obs_res_df.sort_values(
                                    by="abs_ppm", ascending=True, inplace=True
                                )
                                obs_res_df = fit_obs_res_df.head(1)
                            else:
                                obs_res_df = obs_matches_df.head(1)
                        else:
                            obs_res_df = obs_matches_df.copy()
                        obs_res_df.reset_index(drop=True, inplace=True)
                        obs_res_dct = obs_res_df.to_dict(orient="index")
                        obs_label = obs_res_dct[0].get("label")
                        obs_name = obs_res_dct[0].get("name")
                        obs_ppm = round(obs_res_dct[0].get("ppm"), 1)

                        if sub_res == "frag":
                            if obs_name in res_name_lst:
                                if unique_res >= 1 and obs_name not in observed_residues:
                                    is_count_for_score = True
                                elif unique_res == 1 and res_name in observed_residues:
                                    is_count_for_score = True
                                else:
                                    pass
                                if is_count_for_score:
                                    pr_rank_score += 40 * (1 - rank/10)
                                    observed_fragments[obs_label] = {
                                        "ppm": obs_ppm,
                                        "rank": rank,
                                    }
                                    signal_labels[obs_label] = obs_peak
                                    observed_residues.append(obs_name)
                                else:
                                    if obs_label not in frag_label_dct:
                                        noise_labels[obs_label] = obs_peak
                                    elif (
                                        obs_label in frag_label_dct
                                        and obs_label not in signal_labels
                                    ):
                                        signal_labels[obs_label] = obs_peak
                                        observed_neutral_losses[obs_label] = {
                                            "ppm": obs_ppm,
                                            "rank": rank,
                                        }
                                    else:
                                        pass
                            else:
                                if obs_label not in frag_label_dct:
                                    noise_labels[obs_label] = obs_peak

                        elif sub_res == "nl":
                            if (
                                obs_name in nl_name_lst
                                and unique_res > 1
                                and obs_name not in observed_residues
                            ):
                                is_count_for_score = True
                            elif (
                                obs_name in nl_name_lst
                                and unique_res == 1
                                and obs_name in observed_residues
                            ):
                                is_count_for_score = True
                            else:
                                pass
                            if is_count_for_score:
                                pr_rank_score += 10 * (10 - rank) * 0.1
                                observed_neutral_losses[obs_label] = {
                                    "ppm": obs_ppm,
                                    "rank": rank,
                                }
                                signal_labels[obs_label] = obs_peak
                                observed_residues.append(obs_name)
                            else:
                                if obs_label not in nl_label_dct:
                                    noise_labels[obs_label] = obs_peak
                                elif (
                                    obs_label in nl_label_dct
                                    and obs_label not in signal_labels
                                ):
                                    signal_labels[obs_label] = obs_peak
                                    observed_neutral_losses[obs_label] = {
                                        "ppm": obs_ppm,
                                        "rank": rank,
                                    }
                                else:
                                    pass
                        else:
                            if obs_label not in nl_label_dct:
                                noise_labels[obs_label] = obs_peak
                    else:
                        pass  # skip lower signals

        class_rank_obs_dct = rank_obs_dct.get("class", {})
        raw_class_obs_df_lst = []
        if class_rank_obs_dct:
            for c_r in class_rank_obs_dct:
                c_obs_df = class_rank_obs_dct[c_r].get("obs_df", pd.DataFrame)
                if not c_obs_df.empty:
                    c_obs_df["mz"] = class_rank_obs_dct[c_r].get("obs_mz", 0)
                    c_obs_df["i"] = class_rank_obs_dct[c_r].get("obs_i", 0)
                    raw_class_obs_df_lst.append(c_obs_df)

        class_obs_df_lst = [
            c_df for c_df in raw_class_obs_df_lst if isinstance(c_df, pd.DataFrame)
        ]

        if class_obs_df_lst:
            class_obs_df = pd.concat(class_obs_df_lst)
            class_obs_df.sort_values(by="i", ascending=False, inplace=True)
            class_obs_df.drop_duplicates(subset="label", keep="first", inplace=True)
            class_labels = list(class_label_dct.keys())
            matched_class_obs_df = class_obs_df[
                class_obs_df["label"].isin(class_labels)
            ]
            for c_i, c_r in matched_class_obs_df.iterrows():
                signal_labels[c_r["label"]] = [c_r["mz"], c_r["i"]]
                observed_specific_peaks[c_r["label"]] = {"ppm": c_r["ppm"]}
        else:
            pass

        unique_obs_names = list(set(observed_residues))

        has_all_residues = False
        if len(unique_obs_names) >= unique_res:
            has_all_residues = True
        # if pr_rank_score > 0:
        #     print(f"rank score: {pr_rank_score} | {pr_name} # {'; '.join(unique_obs_names)}")
        checked_noise_labels = {}
        for nl_lb in noise_labels:
            if nl_lb in signal_labels:
                pass
            else:
                checked_noise_labels[nl_lb] = noise_labels[nl_lb]

        if pr_rank_score > 100:
            pr_rank_score = 100
        if pr_rank_score >= rank_score_threshold:
            rank_score_results[pr_candidate] = {
                "name": pr_name,
                "label": pr_candidate,
                "neutral_formula": pr_candidate_info.get("neutral_formula", ""),
                "charged_formula": pr_candidate_info.get("formula", ""),
                "adduct": pr_candidate_info.get("adduct", ""),
                "mz": pr_candidate_info.get("mz", 0),
                "rank_score": pr_rank_score,
                "signal_labels": signal_labels,
                "noise_labels": checked_noise_labels,
                "observed_fragments": observed_fragments,
                "observed_neutral_losses": observed_neutral_losses,
                "observed_specific_peaks": observed_specific_peaks,
                "observed_residues": observed_residues,
                "ms2_df": ms2_df,
                "has_all_residues": has_all_residues,
            }
    else:
        pass

    return rank_score_results


def hunt_msp_score(
    pr_candidate: dict, ms2_df: pd.DataFrame, score_summary: dict, ms2_ppm: int = 50
):
    msp_score = 0
    msp_signal_labels = {}
    msp_absent_labels = {}
    pr_msp_dct = pr_candidate.get("msp", {})
    msp_obs_df_lst = []
    base_peak_i = ms2_df["i"].max()
    msp_factor_i = 0.001 * base_peak_i
    if pr_msp_dct:
        msp_df = pd.DataFrame(pr_msp_dct)
        msp_df.sort_values(by="mz", ascending=True, inplace=True)
        msp_chk_lst = list(zip(msp_df["mz"], msp_df["i"], msp_df["label"]))
        for chk_peak in msp_chk_lst:
            is_msp_matched = False
            mz = chk_peak[0]
            mz_l = mz * (1 - ms2_ppm * 0.000001)
            mz_h = mz * (1 + ms2_ppm * 0.000001)

            tmp_df = ms2_df.query(f"{mz_l} <= mz <= {mz_h}").copy()
            if tmp_df.shape[0] == 1:
                msp_obs_df_lst.append(tmp_df)
                is_msp_matched = True
            elif tmp_df.shape[0] > 1:
                tmp_df = tmp_df.sort_values(by="i", ascending=False)
                msp_obs_df_lst.append(tmp_df.head(1))
                is_msp_matched = True
            else:
                tmp_df = pd.DataFrame(data={"mz": [0.0], "i": [0.0]})
                msp_obs_df_lst.append(tmp_df)

            if is_msp_matched:
                msp_signal_labels[chk_peak[2]] = [
                    round(chk_peak[0], 4),
                    chk_peak[1] * msp_factor_i,
                ]
            else:
                msp_absent_labels[chk_peak[2]] = [
                    round(chk_peak[0], 4),
                    chk_peak[1] * msp_factor_i,
                ]

        if msp_obs_df_lst:
            obs_score_df = pd.concat(msp_obs_df_lst)

            # fit theo i to obs i
            obs_i_max = obs_score_df["i"].max()
            pr_msp_dct["abs_i"] = [i * 0.001 * obs_i_max for i in pr_msp_dct["i"]]

            # reduce matrix to 1D array
            obs_flat = get_msp_vectors(
                {
                    "mz": obs_score_df["mz"].values.tolist(),
                    "i": obs_score_df["i"].values.tolist(),
                }
            )
            lib_flat = get_msp_vectors(
                {
                    "mz": msp_df["mz"].values.tolist(),
                    "i": msp_df["i"].values.tolist(),
                }
            )
            # get cosine similarity
            if lib_flat and len(lib_flat) == len(obs_flat) and sum(lib_flat) > 0:
                msp_score = 100 * (
                    (1 - spatial.distance.cosine(obs_flat, lib_flat)) ** 2
                )
                msp_score = round(msp_score, 1)
        else:
            msp_score = 0

        if 0 <= msp_score <= 100:
            pass
        elif msp_score > 100:
            msp_score = 100
        else:
            msp_score = 0

    score_summary["msp_score"] = msp_score

    score_summary["signal_labels"] = msp_signal_labels
    score_summary["absent_labels"] = msp_absent_labels

    return score_summary


def hunt_fingerprint_score(
    pr_candidate: dict,
    ms2_df: pd.DataFrame,
    score_summary: dict,
    ms2_ppm: int = 50,
    fp_threshold: float = 0.01,
    amplify: float = 1.25,
):
    fingerprint_score = 0
    # signal_peaks = score_summary.get("signal_peaks", {})

    fp_signal_labels = {}
    fp_absent_labels = {}

    fp_lib_lst = []
    fp_obs_lst = []

    max_i = ms2_df["i"].max()
    fp_threshold = max_i * fp_threshold
    pr_fingerprint_lst = sorted(pr_candidate.get("fingerprints", []))

    for mz in pr_fingerprint_lst:
        is_fp_matched = False
        fp_lib_lst.append(1)
        mz_l = mz * (1 - ms2_ppm * 0.000001)
        mz_h = mz * (1 + ms2_ppm * 0.000001)

        tmp_df = ms2_df.query(f"{mz_l} <= mz <= {mz_h} and i >= {fp_threshold}").copy()

        if tmp_df.shape[0] == 1:
            fp_obs_lst.append(1)
            is_fp_matched = True
        elif tmp_df.shape[0] > 1:
            fp_obs_lst.append(1)
            is_fp_matched = True
        else:
            fp_obs_lst.append(0)
            is_fp_matched = False

        if is_fp_matched:
            fp_signal_labels[f"fp@{mz:.2f}"] = [
                round(mz, 4),
                1,
            ]
        else:
            fp_absent_labels[f"fp@{mz:.2f}"] = [
                round(mz, 4),
                0,
            ]

    sum_fp_lib_count = sum(fp_lib_lst)
    fp_count_amplify = 1 + 0.1 * (sum_fp_lib_count - 10)
    if sum_fp_lib_count > 0:
        fingerprint_score = (
            100
            * (1 - spatial.distance.cosine(np.array(fp_obs_lst), np.array(fp_lib_lst)))
            * amplify
            * fp_count_amplify
        )
        fingerprint_score = round(fingerprint_score, 1)
    if not fingerprint_score > 0:
        fingerprint_score = 0
    if fingerprint_score > 100:
        fingerprint_score = 100

    score_summary["fingerprint_score"] = fingerprint_score
    score_summary["signal_labels"] = fp_signal_labels
    score_summary["absent_labels"] = fp_absent_labels

    return score_summary


def hunt_snr_score(
    score_summary: dict,
    amplify_factor: float = 10.4795,
):
    snr_score = 0
    # raw_signal_peaks = score_summary.get("signal_labels", [{}])
    # raw_noise_peaks = score_summary.get("noise_labels", [{}])
    # signal_labels = dict(ChainMap(*raw_signal_peaks))
    # noise_labels = dict(ChainMap(*raw_noise_peaks))
    #
    # checked_noise_labels = {}
    # for nl_lb in noise_labels:
    #     if nl_lb in signal_labels:
    #         pass
    #     else:
    #         checked_noise_labels[nl_lb] = noise_labels[nl_lb]

    signal_df = pd.concat(
        [
            pd.DataFrame(s_dct).T
            for s_dct in score_summary.get("signal_labels", None)
            if s_dct is not None
        ]
    )
    pre_noise_df = pd.concat(
        [
            pd.DataFrame(n_dct).T
            for n_dct in score_summary.get("noise_labels", None)
            if n_dct is not None
        ]
    )
    signal_df.reset_index(inplace=True)
    signal_df.columns = ["label", "mz", "i"]
    signal_df.sort_values(by="i", ascending=False, inplace=True)
    signal_df.drop_duplicates(subset="label", keep="first", inplace=True)
    signal_df.drop_duplicates(subset="mz", keep="first", inplace=True)
    signal_df.sort_values(by="mz", ascending=False, inplace=True)
    signal_peaks = list(
        zip(signal_df["mz"].values.tolist(), signal_df["i"].values.tolist())
    )
    score_summary["signal_labels"] = dict(
        zip(signal_df["label"].values.tolist(), signal_peaks)
    )
    signal_sum_i = signal_df["i"].sum()

    if not pre_noise_df.empty:
        pre_noise_df.reset_index(inplace=True)
        pre_noise_df.columns = ["label", "mz", "i"]
        noise_df = pre_noise_df[~pre_noise_df["label"].isin(signal_df["label"])].copy()
        noise_df.sort_values(by="i", ascending=False, inplace=True)
        noise_df.drop_duplicates(subset="label", keep="first", inplace=True)
        noise_df.sort_values(by="mz", ascending=False, inplace=True)
        noise_peaks = list(
            zip(noise_df["mz"].values.tolist(), noise_df["i"].values.tolist())
        )
        score_summary["noise_labels"] = dict(
            zip(noise_df["label"].values.tolist(), noise_peaks)
        )
        noise_sum_i = noise_df["i"].sum()
    else:
        score_summary["noise_labels"] = {}
        noise_sum_i = 1

    if not noise_sum_i > 0:
        noise_sum_i = 1

    # use the SNR equation SNR = 20 * log10(signal/noise)
    # snr_score = 20 * math.log10((signal_sum_i / noise_sum_i))
    # set s/n == 25 --> SNR_SCORE = 100
    # default 3.5767 = 100 / (20 * math.log10(25)) --> 3.5767
    # if set s/n == 3 --> SNR_SCORE = 100
    # default 3.5767 = 100 / (20 * math.log10(3)) --> 10.4795
    sn_ratio = signal_sum_i / noise_sum_i
    snr_score = amplify_factor * 20 * math.log10(sn_ratio)

    if snr_score < 0:
        snr_score = 0.0
    elif snr_score > 100.0:
        snr_score = 100
    else:
        snr_score = round(snr_score, 1)

    score_summary["snr_score"] = snr_score

    return score_summary


def get_total_score(score_info: dict):
    score_lst = []
    score_name_lst = [
        "isotope_score",
        "rank_score",
        "msp_score",
        "fingerprint_score",
        "snr_score",
    ]
    if "msp_score" in score_info:
        msp_score = score_info.get("msp_score", 0)
        if msp_score <= 0:
            score_name_lst.remove("msp_score")

    for score_name in score_name_lst:
        tmp_score = score_info.get(score_name, 0)
        if tmp_score > 0:
            score_lst.append(tmp_score)

    total_score = round(sum(score_lst) / len(score_name_lst), 1)

    return total_score
