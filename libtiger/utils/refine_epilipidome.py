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
import json
import re

import pandas as pd


def refine_epilipidome(space, origins: str):

    if isinstance(space, str):
        with open(space, "r") as f:
            space = json.load(f)
    else:
        pass
    print("Space loaded...")

    origin_df = pd.read_excel(origins, engine="openpyxl")
    lipid_col_name = "lipid"
    headers_lst = origin_df.columns.tolist()
    for col in headers_lst:
        if re.match(r'.*lipid.*', col, re.IGNORECASE):
            lipid_col_name = col
    origin_lst = origin_df[lipid_col_name].tolist()
    print(f"Refine by lipids: {origin_lst}")

    refined_space_dct = {}

    for k in space:
        k_space = {}
        for k_chg in space[k]:
            k_info = space[k].get(k_chg, {})
            mod_type = k_info.get("mod_type", "")
            is_candidate = 0
            if mod_type and mod_type != "unmod":
                k_origin_lst = k_info.get("origins", [])
                for o in origin_lst:
                    if o in k_origin_lst:
                        is_candidate = 1
                    else:
                        pass
            if is_candidate:
                # print(k_info.get("name", ""))
                k_space[k_chg] = k_info
            else:
                pass
        if k_space:
            refined_space_dct[k] = k_space

    return refined_space_dct


if __name__ == "__main__":
    usr_space = r"../../test/input/oxPC.json"
    usr_refine_list = r"../../test/input/unmodified_lipids.xlsx"

    usr_refined_space_dct = refine_epilipidome(usr_space, usr_refine_list)
    print(list(usr_refined_space_dct.keys()))

    print('FIN')
