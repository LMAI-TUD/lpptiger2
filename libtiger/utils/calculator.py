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

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from libtiger.utils.default_values import (
    periodic_table_dct,
    nominal_mass_dct,
    adduct_elements,
)


def get_flat_formula(formula: str):

    while len(re.findall("\(\w*\)", str(formula))) > 0:
        parenthetical = re.findall("\(\w*\)[0-9]+", formula)
        for i in parenthetical:
            p = re.findall("[0-9]+", str(re.findall("\)[0-9]+", i)))
            j = re.findall("\[2H\][0-9]*|\[13C\][0-9]*|[A-Z][a-z]*[0-9]*", i)
            for n in range(0, len(j)):
                numero = re.findall("[0-9]+", j[n])
                if len(numero) != 0:
                    for k in numero:
                        nu = re.sub(k, str(int(int(k) * int(p[0]))), j[n])
                else:
                    nu = re.sub(j[n], j[n] + p[0], j[n])
                j[n] = nu
            newphrase = ""
            for m in j:
                newphrase += str(m)
            formula = formula.replace(i, newphrase)
        if (len((re.findall("\(\w*\)[0-9]+", formula))) == 0) and (
            len(re.findall("\(\w*\)", formula)) != 0
        ):
            formula = formula.replace(r"(", "")
            formula = formula.replace(r")", "")

    return formula


def get_neutral_elem(formula: str):

    elem_dct = {}
    elem_lst = re.findall(
        r"(?P<elem>\[2H\]|\[13C\]|[A-Z][a-z]*)(?P<count>[0-9]*)?", formula
    )
    for _e in elem_lst:
        _elem = _e[0]
        _elem_count = _e[1]
        if len(_elem_count) == 0:
            _elem_count = 1
        else:
            _elem_count = int(_elem_count)
        if _elem in elem_dct:
            elem_dct[_elem] += _elem_count
        else:
            elem_dct[_elem] = _elem_count

    return elem_dct


def get_charged_elem(formula: str, charge="[M-H]-"):

    elem_dct = get_neutral_elem(formula)
    if charge in adduct_elements.keys():
        for k in adduct_elements[charge][0].keys():
            if k in elem_dct:
                elem_dct[k] += adduct_elements[charge][0][k]
            else:
                elem_dct[k] = adduct_elements[charge][0][k]
    else:
        print(f"Not supported charge name: {charge}")
        print(f"Supported types: {list(adduct_elements.keys())}")

    return elem_dct


def get_elements(formula: str, charge: str = None):

    if charge in ["neutral", "Neutral", "", None]:
        elem_dct = get_neutral_elem(formula)
    else:
        elem_dct = get_charged_elem(formula, charge=charge)

    return elem_dct


def get_formula(formula: str, charge=""):

    formula = get_flat_formula(formula)
    elem_dct = get_elements(formula, charge=charge)

    formula_str = ""
    for elem in ["C", "[13C]", "H", "D", "[2H]", "O", "N", "P", "S", "Na", "K"]:
        if elem in elem_dct:
            if elem_dct[elem] == 1:
                formula_str += elem
            elif elem_dct[elem] > 1:
                formula_str += f"{elem}{elem_dct[elem]}"
    if charge not in ["Neutral", "neutral", "", None]:
        formula_str += adduct_elements.get(charge, ["", ""])[1]
    return formula_str


def get_nl_formula(pr_formula: str, nl_formula: str, pr_adduct: str, nl_ion_adduct=""):

    pr_adduct_elements = adduct_elements.get(pr_adduct, [{}, ""])[0]
    nl_ion_adduct_elem, nl_ion_adduct_chg = adduct_elements.get(nl_ion_adduct, [{}, ""])

    pr_formula = get_flat_formula(pr_formula)
    nl_formula = get_flat_formula(nl_formula)
    pr_elem_dct = get_elements(pr_formula, charge="")
    nl_elem_dct = get_elements(nl_formula, charge="")

    nl_ion_elem_dct = {}
    all_elem_lst = list(
        set(
            list(pr_elem_dct.keys())
            + list(pr_adduct_elements.keys())
            + list(nl_elem_dct.keys())
            + list(nl_ion_adduct_elem.keys())
        )
    )
    for elem in all_elem_lst:
        nl_ion_elem_dct[elem] = (
            pr_elem_dct.get(elem, 0)
            - pr_adduct_elements.get(elem, 0)
            - nl_elem_dct.get(elem, 0)
            + nl_ion_adduct_elem.get(elem, 0)
        )

    formula_str = ""
    for elem in ["C", "[13C]", "H", "D", "[2H]", "O", "N", "P", "S", "Na", "K"]:
        if elem in nl_ion_elem_dct:
            if nl_ion_elem_dct[elem] == 1:
                formula_str += elem
            elif nl_ion_elem_dct[elem] > 1:
                formula_str += f"{elem}{nl_ion_elem_dct[elem]}"
            else:
                pass

    formula_str += nl_ion_adduct_chg

    return formula_str


def calc_exact_mass(elem_dct, decimals: int = 6):
    exact_mass = 0.0
    for _elem in elem_dct:
        exact_mass += elem_dct[_elem] * periodic_table_dct[_elem][0][0]

    return round(exact_mass, decimals)


def get_exact_mass(formula: str, charge: str = "neutral", decimals: int = 6):

    elem_dct = get_elements(formula, charge=charge)
    exact_mass = calc_exact_mass(elem_dct, decimals)

    return exact_mass


def calc_nominal_mass(elem_dct):
    nominal_mass = 0
    for _elem in elem_dct:
        nominal_mass += elem_dct[_elem] * nominal_mass_dct.get(_elem, 0)

    return nominal_mass


def get_nominal_mass(formula: str, charge: str = "neutral"):

    elem_dct = get_elements(formula, charge=charge)
    exact_mass = calc_nominal_mass(elem_dct)

    return exact_mass


def get_charged_formula_mz(formula: str, charge: str = "neutral"):

    formula_str = get_formula(formula, charge=charge)
    elem_dct = get_elements(formula, charge=charge)
    mz = calc_exact_mass(elem_dct)

    return {"formula": formula_str, "mz": mz}


def get_ppm_range(exact_mass: float, ppm: float = 5, decimals: int = 6):

    factor_lower = 1 - 0.000001 * ppm
    factor_higher = 1 + 0.000001 * ppm
    ppm_range = f"{exact_mass * factor_lower:.{decimals}f}<= +/- {ppm} ppm <={exact_mass * factor_higher:.{decimals}f}"

    return ppm_range


def get_delta_range(exact_mass: float, delta: float = 0.5, decimals: int = 6):

    delta_range = f"{exact_mass - delta:.{decimals}f}<= +/- {delta} Delta <={exact_mass + delta:.{decimals}f}"

    return delta_range


def calc_smi(
    lipid_class: str, smiles: str = r"OC(CCCCCCC\C=C/C\C=C/CCCCC)=O", decimals: int = 6
):

    mol = Chem.MolFromSmiles(smiles)
    mol_formula = rdMolDescriptors.CalcMolFormula(mol)
    mol_exact_mass = round(rdMolDescriptors.CalcExactMolWt(mol), decimals)
    calc_all_charge = False
    mz_dct = {}
    if isinstance(lipid_class, str):
        if lipid_class.upper() in ["FA", "SN", "RES", "RESIDUE"]:
            mz_dct = {
                "M": get_charged_formula_mz(mol_formula, charge="neutral"),
                "[M-H2O]": get_charged_formula_mz(mol_formula, charge="[M-H2O]"),
                "[M-H]-": get_charged_formula_mz(mol_formula, charge="[M-H]-"),
                "[M-H2O-H]-": get_charged_formula_mz(mol_formula, charge="[M-H2O-H]-"),
                "[M+H]+": get_charged_formula_mz(mol_formula, charge="[M+H]+"),
                "[M-H2O+H]+": get_charged_formula_mz(mol_formula, charge="[M-H2O+H]+"),
                "[M+NH4]+": get_charged_formula_mz(mol_formula, charge="[M+NH4]+"),
                "[M+Na]+": get_charged_formula_mz(mol_formula, charge="[M+Na]+"),
                "[M-H2O+Na]+": get_charged_formula_mz(mol_formula, charge="[M-H2O+Na]+"),
            }
        elif lipid_class.upper() in [
            "PA",
            "PE",
            "PG",
            "PI",
            "PS",
            "LPA",
            "LPE",
            "LPG",
            "LPI",
            "LPS",
        ]:
            mz_dct = {
                "M": get_charged_formula_mz(mol_formula, charge="neutral"),
                "[M-H]-": get_charged_formula_mz(mol_formula, charge="[M-H]-"),
                "[M-H2O-H]-": get_charged_formula_mz(mol_formula, charge="[M-H2O-H]-"),
            }
        elif lipid_class.upper() in ["PC", "LPC"]:
            mz_dct = {
                "M": get_charged_formula_mz(mol_formula, charge="neutral"),
                "[M-H]-": get_charged_formula_mz(mol_formula, charge="[M-H]-"),
                "[M+HCOO]-": get_charged_formula_mz(mol_formula, charge="[M+HCOO]-"),
                "[M+CH3COO]-": get_charged_formula_mz(
                    mol_formula, charge="[M+CH3COO]-"
                ),
            }
        elif lipid_class.upper() in ["DG", "TG"]:
            mz_dct = {
                "M": get_charged_formula_mz(mol_formula, charge="neutral"),
                "[M+H]+": get_charged_formula_mz(mol_formula, charge="[M+H]+"),
                "[M+NH4]+": get_charged_formula_mz(mol_formula, charge="[M+NH4]+"),
                "[M+Na]+": get_charged_formula_mz(mol_formula, charge="[M+Na]+"),
            }
        elif lipid_class.upper() in ["CE"]:
            mz_dct = {
                "M": get_charged_formula_mz(mol_formula, charge="neutral"),
                "[M+Na]+": get_charged_formula_mz(mol_formula, charge="[M+Na]+"),
            }
        else:
            calc_all_charge = True
    else:
        calc_all_charge = True
    if calc_all_charge:
        mz_dct = {
            "M": get_charged_formula_mz(mol_formula, charge="neutral"),
            "[M-H2O]": get_charged_formula_mz(mol_formula, charge="[M-H2O]"),
            "[M+H]+": get_charged_formula_mz(mol_formula, charge="[M+H]+"),
            "[M-H2O+H]+": get_charged_formula_mz(mol_formula, charge="[M-H2O+H]+"),
            "[M+Na]+": get_charged_formula_mz(mol_formula, charge="[M+Na]+"),
            "[M-H2O+Na]+": get_charged_formula_mz(mol_formula, charge="[M-H2O+Na]+"),
            "[M+NH4]+": get_charged_formula_mz(mol_formula, charge="[M+NH4]+"),
            "[M-H]-": get_charged_formula_mz(mol_formula, charge="[M-H]-"),
            "[M-H2O-H]-": get_charged_formula_mz(mol_formula, charge="[M-H2O-H]-"),
            "[M+HCOO]-": get_charged_formula_mz(mol_formula, charge="[M+HCOO]-"),
            "[M+CH3COO]-": get_charged_formula_mz(mol_formula, charge="[M+CH3COO]-"),
        }
    else:
        pass
    return {"formula": mol_formula, "exact_mass": mol_exact_mass, "mz_info": mz_dct}


def get_kendrick_mass(
    formula: str, exact_mass: float = 0.0, base: str = "H", decimals: int = 4
):
    """
    ! kendrick_mass is NOT kendrick_mass_defect (KMD)
    kendrick_mass = exact_mass * (nominal_mass of base / exact_mass of base)
    base is the base unit: e.g. H, O, CH2
    e.g. C28H48O7NP with base of H
    kendrick_mass = 541.316839 * (1 / 1.0078250321) = 537.1139
    """
    if exact_mass > 0:
        pass
    else:
        if isinstance(formula, str):
            exact_mass = get_exact_mass(formula)
            print(exact_mass)
        else:
            exact_mass = 0.0
    kendrick_mass = 0.0
    if base in nominal_mass_dct and base in periodic_table_dct:
        base_nominal_mass = nominal_mass_dct.get(base)
        base_exact_mass = periodic_table_dct.get(base)[0][0]
        if base_nominal_mass and base_exact_mass:
            kendrick_mass = exact_mass * base_nominal_mass / base_exact_mass
    kendrick_mass = round(kendrick_mass, decimals)
    print(kendrick_mass)
    return kendrick_mass


def get_kmd(formula: str, base: str = "H", decimals: int = 4, absolute: bool = False):
    """
    kendrick_mass_defect (KMD) = nominal_mass - kendrick_mass
    base is the base unit: e.g. H, O, CH2
    e.g. C28H48O7NP with base of H
    nominal_mass = 541
    kendrick_mass = 541.316839 * (1 / 1.0078250321) = 537.1139
    KMD(H) = 541 - 537.1139 = 3.8861
    """

    # print(formula)
    nominal_mass = get_nominal_mass(formula=formula)
    print(nominal_mass)
    kendrick_mass = get_kendrick_mass(formula=formula, base=base)
    kmd = nominal_mass - kendrick_mass

    if absolute:
        kmd = abs(kmd)

    kmd = round(kmd, decimals)
    return kmd


if __name__ == "__main__":

    usr_bulk_abbr_lst = ["C28H48NO8P", "C28H48NO7P", "C28H50NO7P", "C28H52NO7P"]

    for usr_formula in usr_bulk_abbr_lst:
        usr_exact_mass = get_charged_formula_mz(usr_formula)
        usr_elements = get_neutral_elem(usr_formula)
        print(usr_exact_mass)
        print(usr_elements)
        print(f'KMD(H): {get_kmd(usr_formula, "H")}')
        print(f'KMD(O): {get_kmd(usr_formula, "O")}')
