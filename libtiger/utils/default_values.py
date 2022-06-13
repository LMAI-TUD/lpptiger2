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
import regex as re

# iupac '97
periodic_table_dct = {
    "H": [(1.0078250321, 0.999885), (2.0141017780, 0.0001157)],
    "D": [(2.0141017780, 0.0001157)],
    "[2H]": [(2.0141017780, 0.0001157)],
    "C": [(12.0, 0.9893), (13.0033548378, 0.0107)],
    "[13C]": [(13.0033548378, 0.0107)],
    "N": [(14.0030740052, 0.99632), (15.0001088984, 0.00368)],
    "O": [
        (15.9949146221, 0.99757),
        (16.99913150, 0.00038),
        (17.9991604, 0.00205),
    ],
    "Na": [(22.98976967, 1.0)],
    "P": [(30.97376151, 1.0)],
    "S": [
        (31.97207069, 0.9493),
        (32.97145850, 0.0076),
        (33.96786683, 0.0429),
        (35.96708088, 0.0002),
    ],
    "K": [
        (38.9637069, 0.932581),
        (39.96399867, 0.000117),
        (40.96182597, 0.067302),
    ],
}

nominal_mass_dct = {
    "H": 1,
    "C": 12,
    "N": 14,
    "O": 16,
    "Na": 23,
    "P": 31,
    "S": 32,
    "K": 39,
}

adduct_elements = {
    "[M-H]-": [{"H": -1}, "-"],
    "[M+HCOO]-": [{"H": 1, "C": 1, "O": 2}, "-"],
    "[M+CH3COO]-": [{"H": 3, "C": 2, "O": 2}, "-"],
    "[M+H]+": [{"H": 1}, "+"],
    "[M+NH4]+": [{"N": 1, "H": 4}, "+"],
    "[M+Na]+": [{"Na": 1}, "+"],
    "[M+K]+": [{"K": 1}, "+"],
    "[M-H2O]": [{"H": -2, "O": -1}, ""],
    "[M-H2O-H]-": [{"H": -3, "O": -1}, "-"],
    "[M-H2O+H]+": [{"H": -1, "O": -1}, "+"],
    "[M-H2O+Na]+": [{"H": -2, "O": -1, "Na": 1}, "+"],
    "[M-CO2-H]-": [{"C": -1, "H": -1, "O": -2}, "-"],
    "[M-CO2-H2O-H]-": [{"C": -1, "H": -3, "O": -3}, "-"],
    "[M-CH3]-": [{"C": -1, "H": -3}, "-"],
    "[M-CH3-H2O]-": [{"C": -1, "H": -5, "O": -1}, "-"],
    "[M+CH3-H]-": [{"C": 1, "H": 2}, "-"],
    # ox adducts from Lipostar v1
    "[M+O-H]-": [{"H": -1, "O": 1}, "-"],
    "[M+O-2H-H]-": [{"H": -3, "O": 1}, "-"],
    "[M+2O-H]-": [{"H": -1, "O": 2}, "-"],
    "[M+2O-2H-H]-": [{"H": -3, "O": 2}, "-"],
    "[M+2O-4H-H]-": [{"H": -5, "O": 2}, "-"],
    "[M+3O-H]-": [{"H": -1, "O": 3}, "-"],
    "[M+3O-2H-H]-": [{"H": -3, "O": 3}, "-"],
    "[M+3O-4H-H]-": [{"H": -5, "O": 3}, "-"],
    "[M+4O-H]-": [{"H": -1, "O": 4}, "-"],
    "[M+4O-2H-H]-": [{"H": -3, "O": 4}, "-"],
    "[M+4O-4H-H]-": [{"H": -5, "O": 4}, "-"],
    "[M+O+HCOO]-": [{"H": 1, "C": 1, "O": 3}, "-"],
    "[M+O-2H+HCOO]-": [{"H": -1, "C": 1, "O": 3}, "-"],
    "[M+2O+HCOO]-": [{"H": 1, "C": 1, "O": 4}, "-"],
    "[M+2O-2H+HCOO]-": [{"H": -1, "C": 1, "O": 4}, "-"],
    "[M+2O-4H+HCOO]-": [{"H": -3, "C": 1, "O": 4}, "-"],
    "[M+3O+HCOO]-": [{"H": 1, "C": 1, "O": 5}, "-"],
    "[M+3O-2H+HCOO]-": [{"H": -1, "C": 1, "O": 5}, "-"],
    "[M+3O-4H+HCOO]-": [{"H": -3, "C": 1, "O": 5}, "-"],
    "[M+4O+HCOO]-": [{"H": 1, "C": 1, "O": 6}, "-"],
    "[M+4O-2H+HCOO]-": [{"H": -1, "C": 1, "O": 6}, "-"],
    "[M+4O-4H+HCOO]-": [{"H": -3, "C": 1, "O": 6}, "-"],
    "[M+O+CH3COO]-": [{"H": 3, "C": 2, "O": 3}, "-"],
    "[M+O-2H+CH3COO]-": [{"H": 1, "C": 2, "O": 3}, "-"],
    "[M+2O+CH3COO]-": [{"H": 3, "C": 2, "O": 4}, "-"],
    "[M+2O-2H+CH3COO]-": [{"H": 1, "C": 2, "O": 4}, "-"],
    "[M+2O-4H+CH3COO]-": [{"H": -1, "C": 2, "O": 4}, "-"],
    "[M+3O+CH3COO]-": [{"H": 3, "C": 2, "O": 5}, "-"],
    "[M+3O-2H+CH3COO]-": [{"H": 1, "C": 2, "O": 5}, "-"],
    "[M+3O-4H+CH3COO]-": [{"H": -1, "C": 2, "O": 5}, "-"],
    "[M+4O+CH3COO]-": [{"H": 3, "C": 2, "O": 6}, "-"],
    "[M+4O-2H+CH3COO]-": [{"H": 1, "C": 2, "O": 6}, "-"],
    "[M+4O-4H+CH3COO]-": [{"H": -1, "C": 2, "O": 6}, "-"],
    # ox adducts from Lipostar v2 beta
    "[M-H+O]-": [{"H": -1, "O": 1}, "-"],
    "[M-H+O-2H]-": [{"H": -3, "O": 1}, "-"],
    "[M-H+2O]-": [{"H": -1, "O": 2}, "-"],
    "[M-H+2O-2H]-": [{"H": -3, "O": 2}, "-"],
    "[M-H+2O-4H]-": [{"H": -5, "O": 2}, "-"],
    "[M-H+3O]-": [{"H": -1, "O": 3}, "-"],
    "[M-H+3O-2H]-": [{"H": -3, "O": 3}, "-"],
    "[M-H+3O-4H]-": [{"H": -5, "O": 3}, "-"],
    "[M-H+4O]-": [{"H": -1, "O": 4}, "-"],
    "[M-H+4O-2H]-": [{"H": -3, "O": 4}, "-"],
    "[M-H+4O-4H]-": [{"H": -5, "O": 4}, "-"],
    "[M+HCOO+O]-": [{"H": 1, "C": 1, "O": 3}, "-"],
    "[M+HCOO+O-2H]-": [{"H": -1, "C": 1, "O": 3}, "-"],
    "[M+HCOO+2O]-": [{"H": 1, "C": 1, "O": 4}, "-"],
    "[M+HCOO+2O-2H]-": [{"H": -1, "C": 1, "O": 4}, "-"],
    "[M+HCOO+2O-4H]-": [{"H": -3, "C": 1, "O": 4}, "-"],
    "[M+HCOO+3O]-": [{"H": 1, "C": 1, "O": 5}, "-"],
    "[M+HCOO+3O-2H]-": [{"H": -1, "C": 1, "O": 5}, "-"],
    "[M+HCOO+3O-4H]-": [{"H": -3, "C": 1, "O": 5}, "-"],
    "[M+HCOO+4O]-": [{"H": 1, "C": 1, "O": 6}, "-"],
    "[M+HCOO+4O-2H]-": [{"H": -1, "C": 1, "O": 6}, "-"],
    "[M+HCOO+4O-4H]-": [{"H": -3, "C": 1, "O": 6}, "-"],
    "[M+CH3COO+O]-": [{"H": 3, "C": 2, "O": 3}, "-"],
    "[M+CH3COO+O-2H]-": [{"H": 1, "C": 2, "O": 3}, "-"],
    "[M+CH3COO+2O]-": [{"H": 3, "C": 2, "O": 4}, "-"],
    "[M+CH3COO+2O-2H]-": [{"H": 1, "C": 2, "O": 4}, "-"],
    "[M+CH3COO+2O-4H]-": [{"H": -1, "C": 2, "O": 4}, "-"],
    "[M+CH3COO+3O]-": [{"H": 3, "C": 2, "O": 5}, "-"],
    "[M+CH3COO+3O-2H]-": [{"H": 1, "C": 2, "O": 5}, "-"],
    "[M+CH3COO+3O-4H]-": [{"H": -1, "C": 2, "O": 5}, "-"],
    "[M+CH3COO+4O]-": [{"H": 3, "C": 2, "O": 6}, "-"],
    "[M+CH3COO+4O-2H]-": [{"H": 1, "C": 2, "O": 6}, "-"],
    "[M+CH3COO+4O-4H]-": [{"H": -1, "C": 2, "O": 6}, "-"],
}

pl_hg_dct = {
    "LPA": r"OP(O)(OCC(",
    "LPC": r"[O-]P(OCC[N+](C)(C)C)(OCC(",
    "LPC-CH3": r"[O]P(OCC[N](C)C)(OCC(",
    "LPE": r"OP(OCCN)(OCC(",
    "LPG": r"OP(OCC(O)CO)(OCC(",
    "LPS": r"OP(OCC(C(O)=O)N)(OCC(",
    "LPI": r"OP(O[C@H]1[C@H](O)([C@@H](O)[C@H](O)[C@@H](O)[C@H]1O))(OCC(",
    "LPIP": r"OP(O[C@H]1[C@H](O)([C@@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H]1O))(OCC(",
    "LPI4P": r"OP(O[C@H]1[C@H](O)([C@@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H]1O))(OCC(",
    "PA": r"OP(O)(OCC(",
    "PC": r"[O-]P(OCC[N+](C)(C)C)(OCC(",
    "PC-CH3": r"[O]P(OCC[N](C)C)(OCC(",
    "PE": r"OP(OCCN)(OCC(",
    "PG": r"OP(OCC(O)CO)(OCC(",
    "PS": r"OP(OCC(C(O)=O)N)(OCC(",
    "PI": r"OP(O[C@H]1[C@H](O)([C@@H](O)[C@H](O)[C@@H](O)[C@H]1O))(OCC(",
    "PIP": r"OP(O[C@H]1[C@H](O)([C@@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H]1O))(OCC(",
    "PI4P": r"OP(O[C@H]1[C@H](O)([C@@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H]1O))(OCC(",
}

st_hg_dct = {
    "CE_left": r"[C@]12(CC=C3C[C@@H](",
    "CE_right": r")CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@@]([H])([C@@](C)([H])CCCC(C)C)CC[C@@]21[H])[H]",
}

res_rgx = re.compile(
    r"^(?P<head>OC[(])(?P<left_c>C*)(?P<units>(C\\C=C/)+)(?P<right_c>C*)(?P<tail>[)]=O)$"
)

unit_rgx = re.compile(r"C\\C=C/")

charge_cfg_dct = {
    "LPA": {"[M-H]-": ["[M-H]-", "[M-H2O-H]-"]},
    "LPC": {
        "[M-H]-": ["[M-H]-", "[M-H2O-H]-"],
        "[M+HCOO]-": ["[M-H]-", "[M-H2O-H]-"],
        # "[M+CH3COO]-": ["[M-H]-", "[M-H2O-H]-"],
    },
    "LPE": {"[M-H]-": ["[M-H]-", "[M-H2O-H]-"]},
    "LPG": {"[M-H]-": ["[M-H]-", "[M-H2O-H]-"]},
    "LPS": {"[M-H]-": ["[M-H]-", "[M-H2O-H]-"]},
    "LPI": {"[M-H]-": ["[M-H]-", "[M-H2O-H]-"]},
    "PA": {"[M-H]-": ["[M-H]-", "[M-H2O-H]-"]},
    "PC": {
        "[M-H]-": ["[M-H]-", "[M-H2O-H]-"],
        "[M+HCOO]-": ["[M-H]-", "[M-H2O-H]-"],
        # "[M+CH3COO]-": ["[M-H]-", "[M-H2O-H]-"],
    },
    "PE": {"[M-H]-": ["[M-H]-", "[M-H2O-H]-"]},
    "PG": {"[M-H]-": ["[M-H]-", "[M-H2O-H]-"]},
    "PS": {"[M-H]-": ["[M-H]-", "[M-H2O-H]-"]},
    "PI": {"[M-H]-": ["[M-H]-", "[M-H2O-H]-"]},
    "DG": {
        # "[M+H]+": ["[M+H]+", "[M-H2O+H]+"],
        # "[M+NH4]+": ["[M+H]+", "[M-H2O+H]+", "[M+NH4]+"],
        "[M+Na]+": ["[M+Na]+", "[M-H2O+Na]+", "[M+H]+", "[M-H2O+H]+"],
    },
    "TG": {
        # "[M+H]+": ["[M+H]+", "[M-H2O+H]+"],
        # "[M+NH4]+": ["[M+H]+", "[M-H2O+H]+", "[M+NH4]+"],
        "[M+Na]+": ["[M+Na]+", "[M-H2O+Na]+", "[M+H]+", "[M-H2O+H]+"],
    },
    "CE": {
        "[M+Na]+": ["[M+Na]+", "[M-H2O+Na]+", "[M+H]+", "[M-H2O+H]+"],
    },
}

specific_ions_dct = {
    "CE": {
        "frag_info": {
            "Chol:369#[M+H]+": {
                "name": "Chol:369",
                "adduct": "[M+H]+",
                "label": "Chol:369#[M+H]+",
                "remark": "Cholestene cation",
                "mod_type": "unmod",
                "smiles": "",
                "formula": "C27H45+",
                "mz": 369.352125,
            },
        },
    },
    "PA": {
        "frag_info": {},
        "nl_info": {
            "PA:-98#NL": {
                "name": "PA:-98",
                "adduct": "NL",
                "label": "PA:-98#NL",
                "remark": "-PA head group",
                "mod_type": "unmod",
                "smiles": "",
                "formula": "H3O4P",
                "mz": 97.976898,
            },
        },
    },
    "PC": {
        "frag_info": {
            "PC:168#[M-H]-": {
                "name": "PC:168",
                "adduct": "[M-H]-",
                "label": "PC:168#[M-H]-",
                "remark": "demethylated PC head group [M-H]-",
                "mod_type": "unmod",
                "smiles": "",
                "formula": "C4H11O4NP-",
                "mz": 168.042572,
            },
            "PC:224#[M-H]-": {
                "name": "PC:224",
                "adduct": "[M-H]-",
                "label": "PC:224#[M-H]-",
                "remark": "demethylated PC head group dehydrated glycerol ester [M-H]-",
                "mod_type": "unmod",
                "smiles": "",
                "formula": "C7H15O5NP-",
                "mz": 224.068787,
            },
            "PC:242#[M-H]-": {
                "name": "PC:242",
                "adduct": "[M-H]-",
                "label": "PC:242#[M-H]-",
                "remark": "demethylated PC head group glycerol ester [M-H]-",
                "mod_type": "unmod",
                "smiles": "",
                "formula": "C7H17O6NP-",
                "mz": 242.079352,
            },
        },
        "nl_info": {
            "PC:-60#NL": {
                "name": "PC:-60",
                "adduct": "NL",
                "label": "PC:-60#NL",
                "remark": "-methyl (-CH3) [M-H]-",
                "mod_type": "unmod",
                "smiles": "",
                "formula": "CH3",
                "mz": 0.0,
            },
            "PC:-183#NL": {
                "name": "PC:-183",
                "adduct": "NL",
                "label": "PC:-183#NL",
                "remark": "-PC head group",
                "mod_type": "unmod",
                "smiles": "",
                "formula": "C5H15NO4P",
                "mz": 0.0,
            },
        },
    },
    "PE": {
        "frag_info": {
            "PE:140#[M-H]-": {
                "name": "PE:140",
                "adduct": "[M-H]-",
                "label": "PE:140#[M-H]-",
                "remark": "PE head group [M-H]-",
                "mod_type": "unmod",
                "smiles": "",
                "formula": "C2H7O4NP-",
                "mz": 140.011272,
            },
            "PE:196#[M-H]-": {
                "name": "PE:196",
                "adduct": "[M-H]-",
                "label": "PE:196#[M-H]-",
                "remark": "dehydrated PE HG glycerol ester [M-H]-",
                "mod_type": "unmod",
                "smiles": "",
                "formula": "C5H11O5NP-",
                "mz": 196.037487,
            },
        },
        "nl_info": {
            "PE:-141#NL": {
                "name": "PE:-141",
                "adduct": "NL",
                "label": "PE:-141#NL",
                "remark": "-PE head group",
                "mod_type": "unmod",
                "smiles": "",
                "formula": "C2H8O4NP",
                "mz": 141.019097,
            },
            "PE:-43#NL": {
                "name": "PE:-43",
                "adduct": "NL",
                "label": "PE:-43#NL",
                "remark": "-PE Head Group part (ethenamine)",
                "mod_type": "unmod",
                "smiles": "",
                "formula": "C2H5N",
                "mz": 43.042199,
            },
        },
    },
    "PG": {
        "frag_info": {
            "PG:171#[M-H]-": {
                "name": "PG:171",
                "adduct": "[M-H]-",
                "label": "PG:171#[M-H]-",
                "remark": "PG head group [M-H]-",
                "mod_type": "unmod",
                "smiles": "",
                "formula": "C3H8O6P-",
                "mz": 171.005853,
            },
            "PG:153#[M-H]-": {
                "name": "PG:153",
                "adduct": "[M-H]-",
                "label": "PG:153#[M-H]-",
                "remark": "PG head group - H2O [M-H]-",
                "mod_type": "unmod",
                "smiles": "",
                "formula": "C3H6O5P-",
                "mz": 152.995288,
            },
        },
        "nl_info": {
            "PG:-172#NL": {
                "name": "PG:-172",
                "adduct": "NL",
                "label": "PG:-172#NL",
                "remark": "-PG head group",
                "mod_type": "unmod",
                "smiles": "",
                "formula": "C3H9O6P",
                "mz": 172.013678,
            },
        },
    },
    "PI": {
        "frag_info": {
            "PI:241#[M-H]-": {
                "name": "PI:241",
                "adduct": "[M-H]-",
                "label": "PI:241#[M-H]-",
                "remark": "PI head group [M-H]-",
                "mod_type": "unmod",
                "smiles": "",
                "formula": "C6H10O8P-",
                "mz": 241.011333,
            },
        },
        "nl_info": {
            "PI:-162#NL": {
                "name": "PI:-162",
                "adduct": "NL",
                "label": "PI:-162#NL",
                "remark": "-inositol",
                "mod_type": "unmod",
                "smiles": "",
                "formula": "C6H10O5",
                "mz": 162.052824,
            },
        },
    },
    "PS": {
        "frag_info": {
            "PS:184#[M-H]-": {
                "name": "PS:184",
                "adduct": "[M-H]-",
                "label": "PS:184#[M-H]-",
                "remark": "PS head group [M-H]-",
                "mod_type": "unmod",
                "smiles": "",
                "formula": "C3H7NO6P-",
                "mz": 184.001102,
            },
        },
        "nl_info": {
            "PS:-87#NL": {
                "name": "PS:-87",
                "adduct": "NL",
                "label": "PS:-87#NL",
                "remark": "-PS head group part (serine-H2O)",
                "mod_type": "unmod",
                "smiles": "",
                "formula": "C3H5NO2",
                "mz": 87.032029,
            },
        },
    },
}

#
# specific_ions_dct = {
#     "neg_mode": {
#         "PA": {
#             "frag_info": {},
#             "nl_info": {
#                 "PA:-98": {
#                     "name": "PA:-98",
#                     "label": "PA:-98",
#                     "remark": "-PA head group",
#                     "mod_type": "unmode",
#                     "smiles": "",
#                     "formula": "H3O4P",
#                     "exact_mass": 97.976898,
#                     "mz_info": {
#                         "M": {"formula": "H3O4P", "mz": 97.976898},
#                     },
#                 },
#             },
#         },
#         "PC": {
#             "frag_info": {
#                 "PC:168#[M-H]-": {
#                     "name": "PC:168#[M-H]-",
#                     "label": "PC:168",
#                     "remark": "demethylated PC head group [M-H]-",
#                     "mod_type": "unmode",
#                     "smiles": "",
#                     "formula": "C4H12O4NP",
#                     "exact_mass": 169.050397,
#                     "mz_info": {
#                         "[M-H]-": {"formula": "C4H11O4NP-", "mz": 168.042572},
#                     },
#                 },
#                 "PC:224#[M-H]-": {
#                     "name": "PC:224#[M-H]-",
#                     "label": "PC:224",
#                     "remark": "demethylated PC head group dehydrated glycerol ester [M-H]-",
#                     "mod_type": "unmode",
#                     "smiles": "",
#                     "formula": "C7H16O5NP",
#                     "exact_mass": 225.076612,
#                     "mz_info": {
#                         "[M-H]-": {"formula": "C7H15O5NP-", "mz": 224.068787},
#                     },
#                 },
#                 "PC:242#[M-H]-": {
#                     "name": "PC:242#[M-H]-",
#                     "label": "PC:242",
#                     "remark": "demethylated PC head group glycerol ester [M-H]-",
#                     "mod_type": "unmode",
#                     "smiles": "",
#                     "formula": "C7H18O6NP",
#                     "exact_mass": 243.087177,
#                     "mz_info": {
#                         "[M-H]-": {"formula": "C7H17O6NP-", "mz": 242.079352},
#                     },
#                 },
#             },
#             "nl_info": {
#                 "PC:-60": {
#                     "name": "PC:-60",
#                     "label": "PC:-60",
#                     "remark": "-methyl formate (-CH3COOH)",
#                     "mod_type": "unmode",
#                     "smiles": "",
#                     "formula": "C2H4O2",
#                     "exact_mass": 60.021130,
#                     "mz_info": {
#                         "M": {"formula": "C2H4O2", "mz": 60.021130},
#                     },
#                 },
#                 "PC:-183": {
#                     "name": "PC:-183",
#                     "label": "PC:-183",
#                     "remark": "-PC head group",
#                     "mod_type": "unmode",
#                     "smiles": "",
#                     "formula": "C5H14NO4P",
#                     "exact_mass": 183.066047,
#                     "mz_info": {
#                         "M": {"formula": "C5H14NO4P", "mz": 183.066047},
#                     },
#                 },
#             },
#         },
#         "PE": {
#             "frag_info": {
#                 "PE:140#[M-H]-": {
#                     "name": "PE:140#[M-H]-",
#                     "label": "PE:140",
#                     "remark": "PE head group [M-H]-",
#                     "mod_type": "unmode",
#                     "smiles": "",
#                     "formula": "C2H8O4NP",
#                     "exact_mass": 141.019097,
#                     "mz_info": {
#                         "[M-H]-": {"formula": "C2H7O4NP-", "mz": 140.011272},
#                     },
#                 },
#                 "PE:196#[M-H]-": {
#                     "name": "PE:196#[M-H]-",
#                     "label": "PE:196",
#                     "remark": "dehydrated PE HG glycerol ester [M-H]-",
#                     "mod_type": "unmode",
#                     "smiles": "",
#                     "formula": "C5H12O5NP",
#                     "exact_mass": 197.045312,
#                     "mz_info": {
#                         "[M-H]-": {"formula": "C5H11O5NP-", "mz": 196.037487},
#                     },
#                 },
#             },
#             "nl_info": {
#                 "PE:-141": {
#                     "name": "PE:-141",
#                     "label": "PE:-141",
#                     "remark": "-PE head group",
#                     "mod_type": "unmode",
#                     "smiles": "",
#                     "formula": "C2H8O4NP",
#                     "exact_mass": 141.019097,
#                     "mz_info": {
#                         "M": {"formula": "C2H8O4NP", "mz": 141.019097},
#                     },
#                 },
#                 "PE:-43": {
#                     "name": "PE:-43",
#                     "label": "PE:-43",
#                     "remark": "-PE Head Group part (ethenamine)",
#                     "mod_type": "unmode",
#                     "smiles": "",
#                     "formula": "C2H5N",
#                     "exact_mass": 43.042199,
#                     "mz_info": {
#                         "M": {"formula": "C2H5N", "mz": 43.042199},
#                     },
#                 },
#             },
#         },
#         "PG": {
#             "frag_info": {
#                 "PG:171#[M-H]-": {
#                     "name": "PG:171#[M-H]-",
#                     "label": "PG:171",
#                     "remark": "PG head group [M-H]-",
#                     "mod_type": "unmode",
#                     "smiles": "",
#                     "formula": "C3H9O6P",
#                     "exact_mass": 172.013678,
#                     "mz_info": {
#                         "[M-H]-": {"formula": "C3H8O6P-", "mz": 171.005853},
#                     },
#                 },
#                 "PG:153#[M-H]-": {
#                     "name": "PG:153#[M-H]-",
#                     "label": "PG:153",
#                     "remark": "PG head group - H2O [M-H]-",
#                     "mod_type": "unmode",
#                     "smiles": "",
#                     "formula": "C3H7O5P",
#                     "exact_mass": 154.003113,
#                     "mz_info": {
#                         "[M-H]-": {"formula": "C3H6O5P-", "mz": 152.995288},
#                     },
#                 },
#             },
#             "nl_info": {
#                 "PG:-172": {
#                     "name": "PG:-172",
#                     "label": "PG:-172",
#                     "remark": "-PG head group",
#                     "mod_type": "unmode",
#                     "smiles": "",
#                     "formula": "C3H9O6P",
#                     "exact_mass": 172.013678,
#                     "mz_info": {
#                         "[M-H]-": {"formula": "C3H9O6P", "mz": 172.013678},
#                     },
#                 },
#             },
#         },
#         "PS": {
#             "frag_info": {
#                 "PS:184#[M-H]-": {
#                     "name": "PS:184#[M-H]-",
#                     "label": "PS:184",
#                     "remark": "PG head group [M-H]-",
#                     "mod_type": "unmode",
#                     "smiles": "",
#                     "formula": "C3H8NO6P",
#                     "exact_mass": 185.008927,
#                     "mz_info": {
#                         "[M-H]-": {"formula": "C3H7NO6P-", "mz": 184.001102},
#                     },
#                 },
#             },
#             "nl_info": {
#                 "PS:-87": {
#                     "name": "PS:-87",
#                     "label": "PS:-87",
#                     "remark": "-PS head group",
#                     "mod_type": "unmode",
#                     "smiles": "",
#                     "formula": "C3H5NO2",
#                     "exact_mass": 87.032029,
#                     "mz_info": {
#                         "[M-H]-": {"formula": "C3H5NO2", "mz": 87.032029},
#                     },
#                 },
#             },
#         },
#     },
# }
