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
from typing import List

import pandas as pd

from libtiger.utils.calculator import get_formula
from libtiger.utils.default_values import adduct_elements


def peaks_to_df(peaks: str) -> pd.DataFrame:
    peaks = peaks.strip("[]")
    peaks_lst = peaks.split(";")
    peaks_lst = [p.split(",") for p in peaks_lst]
    peak_float_lst = []
    for pkl in peaks_lst:
        peak_float_lst.append([float(pkl[0]), float(pkl[1])])
    peaks_df = pd.DataFrame(data=peak_float_lst, columns=["mz", "i"])
    peaks_df["mz"] = peaks_df["mz"].astype(float)
    peaks_df["i"] = peaks_df["i"].astype(float)

    return peaks_df


def parse_ms2_cell(ms2: str):
    # e.g. ms2 = "([M-H+2O]-)[52.7992,1519.44;57.7156,1496.90;71.3292,2384.19;71.3324]"
    ms2_pkl_str = ""
    if ")" in ms2:
        ms2_info_lst = ms2.split(")")
        if len(ms2_info_lst) == 2:
            ms2_pkl_str = ms2_info_lst[1]

    return ms2_pkl_str


def extract_star_spectra(
    csv: str, rt_range: List[float] = None
) -> (pd.DataFrame, dict, pd.DataFrame):

    if isinstance(rt_range, list) and len(rt_range) == 2:
        pass
    else:
        rt_range = [-1, 999]

    csv_df = pd.read_csv(csv)
    ms1_xic_df = pd.DataFrame()
    spec_idx = 0
    dda_event_idx = 0
    spec_idx_lst = []
    dda_event_lst = []
    spec_dct = {}
    rt_lst = []
    dda_rank_lst = []
    scan_id_lst = []
    pr_mz_lst = []
    polarity_lst = []
    pr_chg_formula_lst = []
    star_id_lst = []
    isotope_score_lst = []
    csv_df.drop_duplicates(
        subset=["ID", "Iso. Pat.Score", "Mass delta ppm"], keep="first", inplace=True
    )
    for i, r in csv_df.iterrows():
        ms1_str = r["MS1"]
        ms2_str = parse_ms2_cell(r["MS2"])
        scan_rt = r["RT"]
        pr_mz = r["m/z"]
        pr_formula_neutral = r["Formula"]
        pr_adduct = r["Adduct"]
        mass_score = r["Mass Score"]
        isotope_score = r["Iso. Pat.Score"]
        star_id = r["ID"]

        pr_charge_mode = pr_adduct[-1]
        if pr_charge_mode == "-":
            polarity = "neg"
        elif pr_charge_mode == "+":
            polarity = "pos"
        else:
            polarity = "neg"

        if isotope_score > 0:
            pass
        else:
            if mass_score > 0:
                isotope_score = max(mass_score * 0.80, 80)
            else:
                isotope_score = 80
        if pr_adduct in adduct_elements:
            pr_formula_charged = get_formula(pr_formula_neutral, pr_adduct)
            if pr_formula_charged not in pr_chg_formula_lst and star_id not in star_id_lst:
                if "," in ms1_str and "," in ms2_str:

                    ms1_df = peaks_to_df(ms1_str)
                    ms2_df = peaks_to_df(ms2_str)
                    if not ms1_df.empty and not ms2_df.empty:
                        dda_top = 0
                        ms1_df.loc[:, "rt"] = float(scan_rt)
                        ms2_df.loc[:, "rt"] = float(scan_rt)
                        for mz_df in [ms1_df, ms2_df]:

                            spec_dct[spec_idx] = mz_df
                            dda_event_lst.append(int(dda_event_idx))
                            spec_idx_lst.append(int(spec_idx))
                            rt_lst.append(float(scan_rt))
                            dda_rank_lst.append(dda_top)
                            scan_id_lst.append(int(spec_idx))
                            pr_mz_lst.append(float(pr_mz))
                            polarity_lst.append(polarity)
                            star_id_lst.append(star_id)
                            isotope_score_lst.append(isotope_score)
                            pr_chg_formula_lst.append(pr_formula_charged)

                            spec_idx += 1
                            dda_top += 1

                        dda_event_idx += 1

                        ms1_xic_df = ms1_xic_df.append(ms1_df)

    scan_info_dct = {
        "dda_event_idx": dda_event_lst,
        "scan_time": rt_lst,
        "spec_index": spec_idx_lst,
        "DDA_rank": dda_rank_lst,
        "scan_number": scan_id_lst,
        "MS2_PR_mz": pr_mz_lst,
        "polarity": polarity_lst,
        "MS2_PR_Formula": pr_chg_formula_lst,
        "star_id": star_id_lst,
        "isotope_score": isotope_score_lst,
    }
    scan_info_df = pd.DataFrame(
        data=scan_info_dct,
        columns=[
            "dda_event_idx",
            "spec_index",
            "scan_time",
            "DDA_rank",
            "scan_number",
            "MS2_PR_mz",
            "polarity",
            "MS2_PR_Formula",
            "star_id",
            "isotope_score",
        ],
    )

    scan_info_df.sort_values(by="scan_time", inplace=True)
    if len(rt_range) == 2:  # use rt range filter
        scan_info_df = scan_info_df.query(
            f"{rt_range[0]} <= scan_time <= {rt_range[1]}"
        ).copy()
    scan_info_df = scan_info_df.round({"MS2_PR_mz": 6})
    int_col_lst = ["dda_event_idx", "spec_index", "DDA_rank", "scan_number"]
    scan_info_df[int_col_lst] = scan_info_df[int_col_lst].astype(int)

    return scan_info_df, spec_dct, ms1_xic_df


if __name__ == "__main__":
    usr_csv = r"../../temp/star_csv/export4tiger_top5.csv"
    rt_window = [1, 30]
    s_info_df, sum_spec_dct, xic_df = extract_star_spectra(usr_csv, rt_window)
    print("Spectra loaded.")
