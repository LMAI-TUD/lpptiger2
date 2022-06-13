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
import os
import time

import pandas as pd


def export_to_lipostar(filepath: str, df: pd.DataFrame):
    is_exported = False
    df.sort_values(
        by=["star_id", "LPPtiger Score"], ascending=[True, False], inplace=True
    )
    df.drop_duplicates(subset=["star_id", "LPPtiger Score"], keep="first", inplace=True)
    star_transfer_df = pd.DataFrame(
        df,
        columns=[
            "star_id",
            "Custom Name",
            "name",
            "adduct",
            "neutral_formula",
            "LPPtiger_score",
            "rank_score",
            "isotope_score",
            "msp_score",
            "fingerprint_score",
            "snr_score",
            "matched_peaks",
            # "matched_residues",
            "identified_peak_list",
        ],
    )

    star_transfer_df.rename(columns={"star_id": "ID"}, inplace=True)
    try:
        star_transfer_df.to_csv(filepath, index=False)
    except Exception as e:
        print(e)

    if os.path.isfile(filepath):
        is_exported = True

    return is_exported


def export_to_xlsx(results: dict, file_path: str):
    has_output_file = False
    output_dct = {
        "Title": [],
        "Proposed_structures": [],
        "Lipid_class": [],
        "Neutral_formula": [],
        "Charged_formula": [],
        "Adduct": [],
        "Lib_mz": [],
        "MS1_observed_mz": [],
        "ppm": [],
        "MS1_observed_i": [],
        "MS2_precursor_mz": [],
        "MS2_scan_time": [],
        "MS2_scan_id": [],
        "Total_score": [],
        "Isotope_score": [],
        "Rank_score": [],
        "MatchSpectra_score": [],
        "Fingerprint_score": [],
        "SNR_score": [],
        "Figure": [],
    }
    print("=== ==> --> Prepare to generate output table")
    output_round_dct = {
        "Lib_mz": 4,
        "MS1_obs_mz": 4,
        "ppm": 1,
        "MS1_obs_i": 1,
        r"MS2_PR_mz": 4,
        "MS2_scan_time": 2,
        "Total_score": 1,
        "Isotope_score": 1,
        "Rank_score": 1,
        "MatchSpectra_score": 1,
        "Fingerprint_score": 1,
        "SNR_score": 1,
    }
    for title in results:
        candidate = results[title]
        output_dct["Title"].append(title)
        output_dct["Proposed_structures"].append(
            candidate.get("lipid", {}).get("lipid_name", "")
        )
        output_dct["Lipid_class"].append(
            candidate.get("lipid", {}).get("lipid_class", "")
        )
        output_dct["Neutral_formula"].append(
            candidate.get("lipid", {}).get("neutral_formula", "")
        )
        output_dct["Charged_formula"].append(
            candidate.get("lipid", {}).get("charged_formula", "")
        )
        output_dct["Adduct"].append(candidate.get("lipid", {}).get("adduct", ""))
        output_dct["Lib_mz"].append(candidate.get("lipid", {}).get("mz", ""))
        output_dct["MS1_observed_mz"].append(
            candidate.get("lipid", {}).get("ms1_observed_mz", 0.0)
        )
        output_dct["ppm"].append(
            candidate.get("lipid", {}).get("ms1_observed_ppm", 0.0)
        )
        output_dct["MS1_observed_i"].append(
            candidate.get("lipid", {}).get("ms1_observed_i", 0.0)
        )
        output_dct["MS2_precursor_mz"].append(
            candidate.get("lipid", {}).get("ms2_precursor_mz", 0.0)
        )
        output_dct["MS2_scan_time"].append(
            candidate.get("spectra", {}).get("ms2", {}).get("scan_time", 0)
        )
        output_dct["MS2_scan_id"].append(
            candidate.get("spectra", {}).get("ms2", {}).get("scan_id", 0)
        )
        output_dct["Total_score"].append(
            candidate.get("scores", {}).get("total_score", 0.0)
        )
        output_dct["Isotope_score"].append(
            candidate.get("scores", {}).get("isotope_score", 0.0)
        )
        output_dct["Rank_score"].append(
            candidate.get("scores", {}).get("rank_score", 0.0)
        )
        output_dct["MatchSpectra_score"].append(
            candidate.get("scores", {}).get("msp_score", 0.0)
        )
        output_dct["Fingerprint_score"].append(
            candidate.get("scores", {}).get("fingerprint_score", 0.0)
        )
        output_dct["SNR_score"].append(
            candidate.get("scores", {}).get("snr_score", 0.0)
        )
        output_dct["Figure"].append(candidate.get("figure", ""))

    output_df = pd.DataFrame(output_dct)

    output_df.sort_values(
        by=["MS1_observed_mz", "MS2_scan_time", "Total_score", "Isotope_score"],
        ascending=[True, True, False, False],
        inplace=True,
    )
    output_df.reset_index(drop=True, inplace=True)
    output_df.index += 1
    final_output_df = output_df.round(output_round_dct).copy()
    try:
        final_output_df.to_excel(file_path, index=False)
        print(file_path)
    except IOError:
        final_output_df.to_excel(
            "%s-%i%s" % (file_path[:-5], int(time.time()), ".xlsx"), index=False
        )
        print(file_path)

    if os.path.isfile(file_path):
        has_output_file = True
        print(f">>> Output table saved as: {file_path}")
    final_output_dct = final_output_df.to_dict(orient="index")

    return has_output_file, final_output_dct
