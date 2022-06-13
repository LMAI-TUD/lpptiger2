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
from copy import deepcopy
import json
import os
import shutil

import pandas as pd
import regex as re


class Result(object):
    def __init__(
        self, lipid_label, ms1_tolerance_ppm: int = 10, ms2_tolerance_ppm: int = 20
    ):
        self.lipid = {
            "lipid_label": lipid_label,
            "ms1_tolerance_ppm": ms1_tolerance_ppm,
            "ms2_tolerance_ppm": ms2_tolerance_ppm,
        }
        if "#" in lipid_label:
            self.lipid["lipid_name"], self.lipid["adduct"] = lipid_label.split("#")

        lipid_class_rgx = re.compile(
            r"(?P<lipid_class>((L|Lyso)?P[ACEGIS])|[TD]A?G|CE)(.*)"
        )
        lipid_class_match = lipid_class_rgx.match(lipid_label)
        if lipid_class_match:
            self.lipid["lipid_class"] = lipid_class_match.group("lipid_class")

        self.scores = {}
        self.spectra = {"ms1": {}, "ms2": {}, "xic": {}}
        self.labels = {}
        self.data = {
            "title": "",
            "figure": "",
            "lipid": self.lipid,
            "scores": self.scores,
            "spectra": self.spectra,
            "labels": self.labels,
        }

    def to_dict(self):
        return deepcopy(self.data)

    def to_json_str(self):
        return json.dumps(self.data)

    def add_figure_name(self):
        pr_mz = self.spectra.get("ms2", {}).get("precursor_mz", 0)
        ms2_scan_time = self.spectra.get("ms2", {}).get("scan_time", 0)
        lipid_name = self.lipid.get("lipid_name", 0)
        adduct = self.lipid.get("adduct", 0)

        lipid_name_for_figure = re.sub(r"[+_,/\\]", r"_", lipid_name)
        lipid_name_for_figure = re.sub(r"[\(<]", r"[", lipid_name_for_figure)
        lipid_name_for_figure = re.sub(r"[\)>]", r"]", lipid_name_for_figure)
        lipid_name_for_figure = re.sub(r"\:", r"-", lipid_name_for_figure)
        adduct_for_figure = re.sub(r"[+-]", r"_", adduct)

        title = f"{pr_mz:.4f} @ {ms2_scan_time:.2f}min | {lipid_name} # {adduct}"
        figure = f"{pr_mz:.4f}_{ms2_scan_time:.2f}min_{lipid_name_for_figure}_{adduct_for_figure}.png"

        self.data["title"] = title
        self.data["figure"] = figure

    def add_xic_data(self, xic_data):
        self.spectra["xic"]["peaks"] = xic_data

    def add_isotope_score_data(self, isotope_score_data):
        self.lipid["ms1_charged_formula"] = isotope_score_data.get(
            "charged_formula", ""
        )
        self.lipid["ms1_observed_mz"] = isotope_score_data.get("observed_mz", 0)
        self.lipid["ms1_observed_i"] = isotope_score_data.get("observed_i", 0)
        self.lipid["ms1_observed_ppm"] = isotope_score_data.get("observed_ppm", 999)
        self.lipid["ms2_precursor_mz"] = isotope_score_data.get("precursor_mz", 0)

        self.scores["isotope_score"] = isotope_score_data.get("isotope_score", 0)

        self.spectra["ms1"]["scan_id"] = isotope_score_data.get("ms1_scan_id", 0)
        self.spectra["ms1"]["scan_time"] = isotope_score_data.get("ms1_scan_time", 0)
        ms1_df = isotope_score_data.get("ms1_df", pd.DataFrame({"mz": [0], "i": [0]}))
        self.spectra["ms1"]["peaks"] = list(
            zip([round(mz, 4) for mz in ms1_df["mz"].to_list()], ms1_df["i"])
        )

        self.spectra["ms2"]["dda_event_id"] = isotope_score_data.get("dda_event_id", 0)
        self.spectra["ms2"]["scan_id"] = isotope_score_data.get("ms2_scan_id", 0)
        self.spectra["ms2"]["scan_time"] = isotope_score_data.get("ms2_scan_time", 0)
        self.spectra["ms2"]["polarity"] = isotope_score_data.get("ms2_polarity", 0)
        self.spectra["ms2"]["precursor_mz"] = isotope_score_data.get("ms2_precursor_mz", 0)
        ms2_df = isotope_score_data.get("ms2_df", pd.DataFrame({"mz": [0], "i": [0]}))
        self.spectra["ms2"]["peaks"] = list(
            zip([round(mz, 4) for mz in ms2_df["mz"].to_list()], ms2_df["i"])
        )

    def add_rank_score_data(self, rank_score_data):
        self.lipid["neutral_formula"] = rank_score_data.get("neutral_formula", "")
        self.lipid["adduct"] = rank_score_data.get("adduct", "")
        self.lipid["residues"] = rank_score_data.get("obs_residues", 0)
        self.lipid["mz"] = rank_score_data.get("mz", 0)

        self.scores["rank_score"] = rank_score_data.get("rank_score", 0)
        self.scores["has_all_residues"] = rank_score_data.get("has_all_residues", False)
        self.scores["observed_residues"] = rank_score_data.get("observed_residues", {})
        self.scores["observed_fragments"] = rank_score_data.get(
            "observed_fragments", {}
        )
        self.scores["observed_neutral_losses"] = rank_score_data.get(
            "observed_neutral_losses", {}
        )
        self.scores["observed_specific_peaks"] = rank_score_data.get(
            "observed_specific_peaks", {}
        )

        self.labels["rank_labels"] = {
            "signal_labels": rank_score_data.get("signal_labels", {}),
            "noise_labels": rank_score_data.get("noise_labels", {}),
        }
        self.labels["snr_labels"] = {
            "signal_labels": rank_score_data.get("signal_labels", {}),
            "noise_labels": rank_score_data.get("noise_labels", {}),
        }

    def add_msp_score_data(self, msp_score_data):
        self.scores["msp_score"] = msp_score_data.get("msp_score", 0)
        self.labels["msp_labels"] = {
            "signal_labels": msp_score_data.get("signal_labels", {}),
            "absent_labels": msp_score_data.get("absent_labels", {}),
        }

    def add_fingerprint_score_data(self, fingerprint_score_data):
        self.scores["fingerprint_score"] = fingerprint_score_data.get(
            "fingerprint_score", 0
        )
        self.labels["fingerprint_labels"] = {
            "signal_labels": fingerprint_score_data.get("signal_labels", {}),
            "absent_labels": fingerprint_score_data.get("absent_labels", {}),
        }

    def add_snr_score_data(self, snr_score_data):
        self.scores["snr_score"] = snr_score_data.get("snr_score", 0)
        self.labels["snr_labels"] = {
            "signal_labels": snr_score_data.get("signal_labels", {}),
            "noise_labels": snr_score_data.get("noise_labels", {}),
        }

    def add_total_score_data(self, total_score):
        self.scores["total_score"] = total_score
        self.add_figure_name()


class LogPageCreator(object):
    def __init__(self, params, results, output_folder, start_time):
        print(os.getcwd())
        self.results = results
        self.output_folder = output_folder
        self.output_img_folder = (
            output_folder + r"/LPPtiger_Results_Figures_%s" % start_time
        )
        self.main_page = output_folder + r"/LPPtiger_Results_%s.html" % start_time
        self.logo = r"LPPtiger_Results_Figures_%s/tiger_logo.ico" % start_time
        _image_lst_page = (
            r"LPPtiger_Results_Figures_%s/LPPtiger_Results_Figures_list.html"
            % start_time
        )
        _params_lst_page = (
            r"LPPtiger_Results_Figures_%s/LPPtiger_Params_list.html" % start_time
        )
        _idx_lst_page = (
            r"LPPtiger_Results_Figures_%s/LPPtiger_Identification_list.html"
            % start_time
        )
        self.image_lst_page = (
            self.output_img_folder + r"/LPPtiger_Results_Figures_list.html"
        )
        self.params_lst_page = self.output_img_folder + r"/LPPtiger_Params_list.html"
        self.idx_lst_page = (
            self.output_img_folder + r"/LPPtiger_Identification_list.html"
        )

        tiger_folder = params["tiger_folder"]

        with open(self.main_page, "w") as _m_page:
            m_info_lst = [
                "<html>\n",
                '<link rel="icon" href="',
                self.logo,
                '" type="image/x-icon"/>\n' "<title>LPPtiger_Results ",
                start_time,
                '</title>\n<frameset cols="390,*">\n<frameset rows="430,*">\n',
                '<frame src="',
                _params_lst_page,
                '" frameborder="0" >\n',
                '<frame src="',
                _idx_lst_page,
                '" frameborder="0" >\n</frameset>\n',
                '<frame src="',
                _image_lst_page,
                '"name ="results_frame">\n</frameset>\n</html>\n',
            ]
            _m_page.write("".join(m_info_lst))

        with open(self.image_lst_page, "w") as _img_page:
            _img_page.write(
                """
                            <html>\n<body>\n<style type="text/css">\n
                            p {margin-left: 16px; text-decoration: none; font-family: sans-serif;}\n
                            h3 {font-size:20px; margin-left: 16px; text-decoration: none; font-family: sans-serif;}\n
                            body {font-family: sans-serif;}\n
                            table, th, td {font-size:14px;text-align: center; font-family: sans-serif;}\n
                            th{background-color:#0066B2;color:white; margin:center;}\n
                            tr:nth-child(odd){background-color: #B1D3EC;}\n
                            tr:nth-child(even){background-color: #7C94A5;}\n
                            a:link {text-decoration:none; color:black} a:hover{text-decoration:underline; color:black;}
                            a:visited {text-decoration:none; color:black;}\n</style>\n"""
            )

        with open(self.params_lst_page, "w") as _params_page:
            _params_page.write(
                """
                                <html>\n<body>\n<style type="text/css">\n
                                p {margin-left: 20px; text-decoration: none; font-family: sans-serif;}\n
                                body {font-family: sans-serif;}\n table{width:100%s;}\n
                                table, th, td {font-size:14px;text-align: center; font-family: sans-serif;}\n
                                th{background-color:#0066B2;color:white; margin:center;}\n
                                a:link {text-decoration:none} a:hover{text-decoration:underline }\n
                                ul {font-size:14px; width: 260px;}\n </style>\n
                                <h3><img src="tiger_logo.ico" height=36/>  LPPtiger</h3><font size="1">\n
                                <hr> <h3>Parameters:</h3>\n<ul>\n
                                <li>Start time: %s</li>\n
                                <li><i>m/z</i> range: %.1f - %.1f <i>m/z</i></li>\n<li>RT range: %.1f - %.1f min</li>\n
                                <li>MS1 Threshold: %i</li>\n<li>MS2 Threshold: %i</li>\n
                                <li>MS1 ppm: %i</li>\n<li>MS2 ppm: %i</li>\n
                                <li>LPPtiger score > %.1f </li>\n<li>Isotope score > %.1f</li>\n
                                </ul>\n<hr>\n<h3>Lipid hunter list:</h3><font size="1">\n<table>\n<thead>\n
                                <tr style="text-align: center;">\n
                                <th>ID#</th>\n<th> MS1_obs_mz </th>\n<th>RT(min)</th>\n<th>Discrete</th>\n
                                <th>Score</th>\n
                                </tr>\n</thead>\n</table>\n</body>\n</html>\n
                                """
                % (
                    "%",
                    start_time,
                    params["mz_start"],
                    params["mz_end"],
                    params["rt_start"],
                    params["rt_end"],
                    params["ms_th"],
                    params["ms2_th"],
                    params["ms_ppm"],
                    params["ms2_ppm"],
                    params["score_filter"],
                    params["isotope_score_filter"]
                )
            )

        with open(self.idx_lst_page, "w") as _idx_page:
            _idx_page.write(
                """
                            <html>\n<body>\n<style type="text/css">\n
                            p {margin-left: 20px; text-decoration: none; font-family: sans-serif;}\n
                            body {background-color: #B1D3EC;font-family: sans-serif;}\n table{width:100%s;}\n
                            table, th, td {font-size:14px;text-align: center; font-family: sans-serif;}\n
                            th{background-color:#0066B2;color:white;margin:center;}\n
                            tr:nth-child(even){background-color: #7C94A5;}\n
                            a:link {text-decoration:none; color:black}a:hover{text-decoration:underline; color:black;}
                            a:visited {text-decoration:none; color:black;}\n
                            </style>\n<table>\n<tbody>
                            """
                % "%"
            )
        logo_path = os.path.join(tiger_folder, "tiger_logo.ico")
        try:
            shutil.copy(logo_path, self.output_img_folder)
        except IOError:
            print(f" can not find logo at {logo_path}")

    def add_all_info(self):

        with open(self.image_lst_page, "a") as img_page:
            with open(self.idx_lst_page, "a") as idx_page:
                print("try to add hunter to report html")
                idx = 1
                for title in self.results:
                    candidate = self.results.get(title, {})
                    ms1_pr_mz = candidate.get("MS1_observed_mz", 0)
                    ms2_rt = candidate.get("MS2_scan_time", 0)
                    ms2_scan_id = candidate.get("Scan_id", 0)
                    ident_abbr = candidate.get("Proposed_structures", "")
                    score = candidate.get("Total_score", 0)
                    formula_ion = candidate.get("Charged_formula", "")
                    ident_abbr = ident_abbr.replace("<", "&lt;")
                    ident_abbr = ident_abbr.replace(">", "&gt;")
                    charge = candidate.get("Adduct")
                    img_title_str = candidate.get("Title")
                    img_path = candidate.get("Figure")
                    ident_idx = str(idx)
                    img_info_lst = [
                        '<a name="',
                        ident_idx,
                        '"><h3>',
                        '<a href="',
                        img_path,
                        '" target="blank">',
                        img_title_str,
                        "</a></h3></a>",
                        '<a href="',
                        img_path,
                        '" target="blank">',
                        '<img src="',
                        img_path,
                        '" height="800" /></a>',
                        "\n<hr>\n",
                    ]
                    img_page.write("".join(img_info_lst))

                    idx_str = """
                            <tr>\n
                            <td>
                            <a href ="LPPtiger_Results_Figures_list.html#{id}" target ="results_frame">{id}
                            </td>\n<td>
                            <a href ="LPPtiger_Results_Figures_list.html#{id}" target ="results_frame">{mz}
                            </td>\n<td>
                            <a href ="LPPtiger_Results_Figures_list.html#{id}" target ="results_frame">{rt}
                            </td>\n<td>
                            <a href ="LPPtiger_Results_Figures_list.html#{id}" target ="results_frame">{ident}
                            </td>\n<td>
                            <a href ="LPPtiger_Results_Figures_list.html#{id}" target ="results_frame">{score}
                            </td>\n</tr>\n
                            """.format(
                        id=ident_idx,
                        mz="%.4f" % ms1_pr_mz,
                        rt="%.1f" % ms2_rt,
                        ident=ident_abbr,
                        score=score,
                    )
                    idx_page.write(idx_str)
                    idx += 1

            print("==> info added to report html -->")

    def close_page(self):
        with open(self.main_page, "a") as _m_page:
            _m_page.write("\n</body></html>\n")

        with open(self.image_lst_page, "a") as _img_page:
            _img_page.write("\n</body></html>\n")

        with open(self.idx_lst_page, "a") as _idx_page:
            _idx_page.write("\n</tbody>\n</table>\n</body></html>\n")

