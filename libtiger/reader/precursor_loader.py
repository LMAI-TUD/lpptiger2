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
import re

from natsort import natsorted
import pandas as pd


def parse_precursors(df: pd.DataFrame) -> (list, dict):
    lipid_class_use_lst = ["PA", "PC", "PE", "PG", "PI", "PIP", "PS"]
    lipid_df = df[df["Converted_Name"].str.contains("\(", na=False)].copy()
    lipid_df["lipid_class"] = lipid_df["Converted_Name"].str.split("(").str[0]
    lipid_df = lipid_df[lipid_df["lipid_class"].isin(lipid_class_use_lst)]
    precursor_df = lipid_df.copy()
    precursor_dct = lipid_df.to_dict(orient="index")
    fa_lst = []
    for pr in precursor_dct:
        pr_info_sct = precursor_dct[pr]
        pr_name_str = pr_info_sct["Converted_Name"]
        pr_name_fas_str = re.sub(r"P\w\(", "", pr_name_str)
        pr_name_fas_str = re.sub(r"\)", "", pr_name_fas_str)
        precursor_dct[pr]["FA_residues"] = natsorted(pr_name_fas_str.split("_"))
        fa_lst.extend(precursor_dct[pr]["FA_residues"])
    lipid_df["FA"] = lipid_df["Converted_Name"].str.split("(").str[0]
    unique_fa_lst = natsorted(set(fa_lst))
    lipid_class_lst = natsorted(lipid_df["lipid_class"].unique().tolist())

    return lipid_class_lst, unique_fa_lst, precursor_dct, precursor_df


if __name__ == "__main__":

    usr_dct = {
        "Converted_Name": ["PC(16:0/18:1)", "PE(18:0/18:2)", "PS(18:1/18:2)", "nan"],
        "Formula": ["C42H82NO8P", "C41H78NO8P", "C42H76NO10P", "BAD"],
    }

    usr_df = pd.DataFrame(data=usr_dct)
    x = parse_precursors(df=usr_df)
    print(x)
