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
#         https://tu-dresden.de/med/mf/zml/forschungsgruppen/fedorova/mitarbeiter-innen-der-fedorova-gruppe
#
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
import re
import time
import multiprocessing
from PySide6 import QtCore, QtWidgets

from libtiger.tiger.tiger_gui import Ui_MainWindow
from libtiger.predictor.prediction import get_predictions
from libtiger.utils.inclusion_list_generator import gen_incl
from libtiger.hunter.hunter_core import hunt


class TigerMainWindow(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self, tiger_cwd: str = None, parent=None):
        # QtWidgets.QMainWindow.__init__(self, parent)
        # self.tiger = Ui_MainWindow()
        # super(Tiger2Star, self).__init__(parent)
        super(TigerMainWindow, self).__init__(parent)
        self.setupUi(self)
        core_count = multiprocessing.cpu_count()
        if core_count > 4:
            self.worker_count = core_count - 2
        else:
            self.worker_count = 2
        # self.hide_unused_widgets()
        # self.load_cfg()
        in_fa_path = r"config/1_FA_list_default.csv"

        if tiger_cwd and os.path.isdir(tiger_cwd):
            self.tiger_folder_path = os.path.abspath(tiger_cwd)
        else:
            tmp_cwd = os.path.abspath(os.getcwd())
            if os.path.isfile(os.path.join(tmp_cwd, in_fa_path)):
                self.tiger_folder_path = os.path.abspath(tmp_cwd)
            elif os.path.isfile(
                os.path.join(tmp_cwd, "tiger_core.py")
            ) or os.path.isfile(os.path.join(tmp_cwd, "tiger_main.exe")):
                os.chdir("../..")
                self.tiger_folder_path = os.path.abspath(os.getcwd())
            else:
                self.tiger_folder_path = os.path.abspath(os.getcwd())

        print(f"Tiger folder: {self.tiger_folder_path}")
        if os.path.isfile(in_fa_path):
            abs_in_fa_path = os.path.abspath(in_fa_path)
            self.input_fa_le.setText(abs_in_fa_path)
            print(f"Default FA list loaded: {abs_in_fa_path}")
        else:
            print(f"Can not find default FA list: {in_fa_path}")

        self.toggle_refine()

        # # define single worker
        self.single_worker = SingleWorker()
        self.single_thread = QtCore.QThread()
        self.single_worker.moveToThread(self.single_thread)
        self.single_worker.workRequested.connect(self.single_thread.start)
        self.single_thread.started.connect(self.single_worker.run_prediction_task)
        self.single_worker.finished.connect(self.single_worker_on_finish)
        self.single_worker.info_update.connect(self.single_worker_info_update)

        # init default values

        # slots for tab a
        self.input_fa_pb.clicked.connect(self.open_input_fa)
        self.refine_epilipidome_chkb.clicked.connect(self.toggle_refine)
        self.refine_epilipidome_pb.clicked.connect(self.open_refine_lipid)
        self.output_epilipidome_pb.clicked.connect(self.save_output_epilipidome)
        self.run_prediction_pb.clicked.connect(self.run_single_prediciton)
        # slots for tab b
        self.input_lipid_pb.clicked.connect(self.open_input_lipid)
        self.input_epilipidome_pb.clicked.connect(self.open_input_epilipidome)
        self.output_incl_list_pb.clicked.connect(self.save_output_incl_list)
        self.run_incl_list_pb.clicked.connect(self.run_incl_list)
        # slots for tab c
        self.input_ident_epilipidome_pb.clicked.connect(self.open_input_epilipidome_2)
        self.input_spectra_mzml_pb.clicked.connect(self.open_input_spectra_mzml)
        self.input_spectra_star_pb.clicked.connect(self.open_input_spectra_mzml)
        self.report_folder_pb.clicked.connect(self.save_report_folder)
        self.export_table_pb.clicked.connect(self.save_export_table)
        self.run_identification_pb.clicked.connect(self.run_identification)

    def toggle_refine(self):
        if self.refine_epilipidome_chkb.isChecked():
            self.refine_epilipidome_le.show()
            self.refine_epilipidome_pb.show()
            self.refine_epilipidome_le.setEnabled(True)
            self.refine_epilipidome_pb.setEnabled(True)
        elif not self.refine_epilipidome_chkb.isChecked():
            self.refine_epilipidome_le.hide()
            self.refine_epilipidome_pb.hide()
            self.refine_epilipidome_le.setDisabled(True)
            self.refine_epilipidome_pb.setDisabled(True)
        else:
            self.refine_epilipidome_le.hide()
            self.refine_epilipidome_pb.hide()
            self.refine_epilipidome_le.setDisabled(True)
            self.refine_epilipidome_pb.setDisabled(True)

    def open_file(self, info_str, lb_obj):
        open_file_dialog = QtWidgets.QFileDialog(self)
        open_file_dialog.setNameFilters([info_str])
        open_file_dialog.selectNameFilter(info_str)
        if open_file_dialog.exec():
            lb_obj.clear()
            file_str = open_file_dialog.selectedFiles()[0]
            file_str = os.path.abspath(file_str)
            lb_obj.setText(file_str)

    @staticmethod
    def save_file(typ_filter, le_obj):
        save_path_tp = QtWidgets.QFileDialog.getSaveFileName(
            caption="Save file", filter=typ_filter
        )
        le_obj.clear()
        save_path_str = os.path.abspath("".join(save_path_tp))
        le_obj.setText(save_path_str)

    @staticmethod
    def open_folder(le_obj):
        le_obj.clear()
        folder_str = QtWidgets.QFileDialog.getExistingDirectory()
        if os.path.isdir(folder_str):
            le_obj.setText(folder_str)

    def open_input_fa(self):
        file_info_str = "FA allowed list (*.csv)"
        self.open_file(file_info_str, self.input_fa_le)

    def open_refine_lipid(self):
        file_info_str = "Input lipids (*.xlsx)"
        self.open_file(file_info_str, self.refine_epilipidome_le)

    def open_input_lipid(self):
        file_info_str = "Input lipids (*.xlsx)"
        self.open_file(file_info_str, self.input_lipid_le)

    def open_input_epilipidome(self):
        file_info_str = "Input predicted epilipidome (*.json)"
        self.open_file(file_info_str, self.input_epilipidome_le)

    def open_input_epilipidome_2(self):
        file_info_str = "Input predicted epilipidome (*.json)"
        self.open_file(file_info_str, self.input_ident_epilipidome_le)

    def open_input_spectra_mzml(self):
        file_info_str = "LC-MS spectra (*.mzml)"
        self.open_file(file_info_str, self.input_spectra_mzml_le)

    def open_input_spectra_star(self):
        file_info_str = "Lipostar 2 export (*.csv)"
        self.open_file(file_info_str, self.input_spectra_star_le)

    def save_output_epilipidome(self):
        file_type_str = ".json"
        self.save_file(file_type_str, self.output_epilipidome_le)

    def save_output_incl_list(self):
        file_type_str = ".csv"
        self.save_file(file_type_str, self.output_incl_list_le)

    def save_report_folder(self):
        report_folder_str = QtWidgets.QFileDialog.getExistingDirectory()
        self.report_folder_le.clear()
        self.report_folder_le.setText(report_folder_str)

    def save_export_table(self):
        file_type_str = ".xlsx"
        self.save_file(file_type_str, self.export_table_le)

    def get_prediction_params(self):
        input_fa_path = self.input_fa_le.text()
        refine_epilipidome_path = self.refine_epilipidome_le.text()
        output_epilipidome_path = self.output_epilipidome_le.text()
        max_site = int(self.max_site_spb.value())
        max_o = int(self.max_o_spb.value())
        max_oh = int(self.max_oh_spb.value())
        max_keto = int(self.max_keto_spb.value())
        max_ooh = int(self.max_ooh_spb.value())
        max_ep = int(self.max_ep_spb.value())

        if (
            self.refine_epilipidome_chkb.isChecked()
            and isinstance(refine_epilipidome_path, str)
            and os.path.isfile(refine_epilipidome_path)
        ):
            pass
        else:
            refine_epilipidome_path = None

        lipid_class_lst = []
        lipid_class_chkb_dct = {
            "LPA": self.lpa_chkb,
            "LPC": self.lpc_chkb,
            "LPE": self.lpe_chkb,
            "LPG": self.lpg_chkb,
            "LPI": self.lpi_chkb,
            "LPS": self.lps_chkb,
            "PA": self.pa_chkb,
            "PC": self.pc_chkb,
            "PE": self.pe_chkb,
            "PG": self.pg_chkb,
            "PI": self.pi_chkb,
            "PS": self.ps_chkb,
            "CE": self.ce_chkb,
            "DG": self.dg_chkb,
            "TG": self.tg_chkb,
        }
        for lipid_class in lipid_class_chkb_dct:
            if isinstance(lipid_class_chkb_dct[lipid_class], QtWidgets.QCheckBox):
                if lipid_class_chkb_dct[lipid_class].isChecked():
                    lipid_class_lst.append(lipid_class)

        usr_params = {
            "lipid_class_lst": lipid_class_lst,
            "input_fa_path": input_fa_path,
            "refine_epilipidome_path": refine_epilipidome_path,
            "output_epilipidome_path": output_epilipidome_path,
            "max_site": max_site,
            "max_o": max_o,
            "max_oh": max_oh,
            "max_keto": max_keto,
            "max_ooh": max_ooh,
            "max_ep": max_ep,
            "mod_cfg": os.path.abspath(
                os.path.join(
                    self.tiger_folder_path, r"config/2_DB_Mod_cfg_lv1.csv"
                )
            ),
            "frag_cfg": os.path.abspath(
                os.path.join(
                    self.tiger_folder_path, r"config/4_FragmentationPatterns.xlsx"
                )
            ),
            "worker_count": self.worker_count
        }
        return usr_params

    def get_incl_list_params(self):
        input_lipid_path = self.input_lipid_le.text()
        input_epilipidome_path = self.input_epilipidome_le.text()
        output_incl_list_path = self.output_incl_list_le.text()
        usr_params = {
            "input_lipid_path": input_lipid_path,
            "input_epilipidome_path": input_epilipidome_path,
            "output_incl_list_path": output_incl_list_path,
        }
        return usr_params

    def get_identification_params(self):
        if self.dump_json_chkb.isChecked():
            save_json = True
        else:
            save_json = False

        dpi_txt = self.dpi_cmb.currentText()
        dpi_val = 300

        if re.search(r'300', dpi_txt):
            dpi_val = 300
        elif re.search(r'150', dpi_txt):
            dpi_val = 150
        elif re.search(r'96', dpi_txt):
            dpi_val = 96

        else:
            dpi_val = 300

        if self.Spectra_tab.currentIndex() == 0:
            spectra_type = "mzml"
            ident_params = {
                "input_epilipidome": self.input_ident_epilipidome_le.text(),
                "input_spectra": self.input_spectra_mzml_le.text(),
                "report_folder_path": self.report_folder_le.text(),
                "export_table_path": self.export_table_le.text(),
                "tiger_folder": self.tiger_folder_path,
                "spectra_type": spectra_type,
                "ms1_ppm": self.ms1_ppm_spb.value(),
                "ms2_ppm": self.ms2_ppm_spb.value(),
                "ms1_abs_threshold": self.ms1_threshold_spb.value(),
                "ms2_relative_threshold": self.ms2_threshold_spb.value() * 0.01,
                "rt_min": self.rt_min_spb.value(),
                "rt_max": self.rt_max_spb.value(),
                "mz_min": self.pr_mz_min_spb.value(),
                "mz_max": self.pr_mz_max_spb.value(),
                "selection_window": self.selection_window_spb.value(),
                "total_score_min": self.total_score_min_spb.value(),
                "isotope_score_min": self.isotope_score_min_spb.value(),
                "rank_score_min": self.rank_score_min_spb.value(),
                "save_json": save_json,
                "dpi": dpi_val,
                "worker_count": self.worker_count,
            }
        elif self.Spectra_tab.currentIndex() == 1:
            spectra_type = "star"
            ident_params = {
                "input_epilipidome": self.input_ident_epilipidome_le.text(),
                "input_spectra": self.input_spectra_mzml_le.text(),
                "report_folder_path": self.report_folder_le.text(),
                "export_table_path": self.export_table_le.text(),
                "tiger_folder": self.tiger_folder_path,
                "spectra_type": spectra_type,
                "ms1_ppm": self.ms1_ppm_spb.value(),
                "ms2_ppm": self.ms2_ppm_spb.value(),
                "ms1_abs_threshold": self.ms1_threshold_spb.value(),
                "ms2_relative_threshold": self.ms2_threshold_spb.value() * 0.01,
                "rt_min": self.rt_min_spb.value(),
                "rt_max": self.rt_max_spb.value(),
                "mz_min": self.pr_mz_min_spb.value(),
                "mz_max": self.pr_mz_max_spb.value(),
                "selection_window": self.selection_window_spb.value(),
                "total_score_min": self.total_score_min_spb.value(),
                "isotope_score_min": self.isotope_score_min_spb.value(),
                "rank_score_min": self.rank_score_min_spb.value(),
                "save_json": save_json,
                "dpi": dpi_val,
                "worker_count": self.worker_count,
            }
        else:
            ident_params = {}

        return ident_params

    def run_prediction(self):
        # cfg_params = self.load_cfg()
        run_params = self.get_prediction_params()
        self.run_prediction_status_te.clear()
        # print(run_params)

        is_ready_to_run = True
        for rk in run_params:
            if isinstance(run_params[rk], str):
                if len(run_params[rk]) > 0:
                    pass
                else:
                    if rk in ["refine_epilipidome_path"]:
                        pass
                    else:
                        is_ready_to_run = False
            elif isinstance(run_params[rk], int) or isinstance(run_params[rk], float):
                if run_params[rk] > 0:
                    pass
                else:
                    if rk in ["max_keto", "max_ooh", "max_ep"]:
                        pass
                    else:
                        is_ready_to_run = False
            elif isinstance(run_params[rk], list):
                if len(run_params[rk]) > 0:
                    pass
                else:
                    is_ready_to_run = False
            else:
                if rk in ["refine_epilipidome_path"]:
                    pass
                else:
                    is_ready_to_run = False
            # print(rk, is_ready_to_run)

        in_fa_path = run_params.get("input_fa_path", "")
        if os.path.isfile(in_fa_path):
            pass
        else:
            is_ready_to_run = False
            self.run_prediction_status_te.append("Can not load the input FA list.")

        print("run params:")
        print(run_params)
        print(f"is_ready_to_run: {is_ready_to_run}")
        if is_ready_to_run:
            print("start prediction...")
            self.run_prediction_status_te.append(f"start prediction for lipid class:")
            self.run_prediction_status_te.append(
                f"{', '.join(run_params.get('lipid_class_lst', [''])) }"
            )
            self.run_prediction_pb.setEnabled(False)
            run_info = get_predictions(run_params)
            self.run_prediction_status_te.append("FINISHED.")
            self.run_prediction_pb.setEnabled(True)
        else:
            self.run_prediction_status_te.append("Please check parameters...")

    def run_incl_list(self):
        self.run_incl_list_status_te.clear()
        run_params = self.get_incl_list_params()
        origins = run_params.get("input_lipid_path")
        space = run_params.get("input_epilipidome_path")
        export_file = run_params.get("output_incl_list_path")
        self.run_incl_list_status_te.append("Start to extract...")
        time.sleep(1)
        skipped_pr_lst = gen_incl(space, origins, export_file)
        if skipped_pr_lst:
            self.run_incl_list_status_te.append("Skipped:")
            self.run_incl_list_status_te.append("\n".join(skipped_pr_lst))
        else:
            pass

        if os.path.isfile(export_file):
            self.run_incl_list_status_te.append("File saved.")
        self.run_incl_list_status_te.append("\n>>> >>> >>> FINISHED <<< <<< <<<\n")

    def run_identification(self):
        ident_params = self.get_identification_params()

        print(ident_params)

        self.run_identification_status_te.clear()

        run_time = hunt(
            spectra_file_path=ident_params.get("input_spectra", ""),
            output_file_path=ident_params.get("export_table_path", ""),
            output_folder_path=ident_params.get("report_folder_path", ""),
            search_space=ident_params.get("input_epilipidome", ""),
            tiger_folder=ident_params.get("tiger_folder", ""),
            spectra_type=ident_params.get("spectra_type", "mzml"),
            mz_range=[
                ident_params.get("mz_min", 300),
                ident_params.get("mz_max", 1200),
            ],
            rt_range=[ident_params.get("rt_min", 3), ident_params.get("rt_max", 30)],
            selection_window=ident_params.get("selection_window", 0.5),
            ms1_ppm=ident_params.get("ms1_ppm", 10),
            ms2_ppm=ident_params.get("ms2_ppm", 20),
            ms1_abs_threshold=ident_params.get("ms1_abs_threshold", 1000),
            ms2_percent_threshold=ident_params.get("ms2_percent_threshold", 0.01),
            min_total_score=ident_params.get("total_score_min", 50),
            min_isotope_score=ident_params.get("isotope_score_min", 80),
            min_rank_score=ident_params.get("rank_score_min", 50),
            min_msp_score=ident_params.get("msp_score_min", 0),
            min_fp_score=ident_params.get("fp_score_min", 0),
            min_snr_score=ident_params.get("snr_score_min", 0),
            save_json=ident_params.get("save_json", False),
            dpi=ident_params.get("dpi", 300),
            worker=self.worker_count
        )
        msg = f"Finished hunter in {run_time}"
        print(msg)
        self.run_identification_status_te.setText(msg)

    def run_single_prediciton(self):
        run_params = self.get_prediction_params()
        self.run_prediction_status_te.clear()

        is_ready_to_run = True
        for rk in run_params:
            if isinstance(run_params[rk], str):
                if len(run_params[rk]) > 0:
                    pass
                else:
                    if rk in ["refine_epilipidome_path"]:
                        pass
                    else:
                        is_ready_to_run = False
            elif isinstance(run_params[rk], int) or isinstance(run_params[rk], float):
                if run_params[rk] > 0:
                    pass
                else:
                    if rk in ["max_keto", "max_ooh", "max_ep"]:
                        pass
                    else:
                        is_ready_to_run = False
            elif isinstance(run_params[rk], list):
                if len(run_params[rk]) > 0:
                    pass
                else:
                    is_ready_to_run = False
            else:
                if rk in ["refine_epilipidome_path"]:
                    pass
                else:
                    is_ready_to_run = False
            # print(rk, is_ready_to_run)

        in_fa_path = run_params.get("input_fa_path", "")
        if os.path.isfile(in_fa_path):
            pass
        else:
            is_ready_to_run = False
            self.run_prediction_status_te.append("Can not load the input FA list.")

        print("run params:")
        print(run_params)
        print(f"is_ready_to_run: {is_ready_to_run}")
        if is_ready_to_run:
            print("start prediction...")
            self.run_prediction_status_te.append(f"Start prediction for lipid class:")
            self.run_prediction_status_te.append(
                f"{', '.join(run_params.get('lipid_class_lst', ['']))}"
            )
            time.sleep(1)
            self.run_prediction_pb.setEnabled(False)
            # run_info = get_predictions(run_params)
            self.single_worker.request_work(run_params)
        else:
            self.run_prediction_status_te.append("Please check parameters...")

    def single_worker_info_update(self):
        back_info_str = self.single_worker.infoback()
        self.run_prediction_status_te.append(back_info_str)

    def single_worker_on_finish(self):
        self.single_thread.quit()
        print("!! single_worker stopped !!")
        self.run_prediction_pb.setEnabled(True)


class SingleWorker(QtCore.QObject):
    workRequested = QtCore.Signal()
    finished = QtCore.Signal()
    info_update = QtCore.Signal(str)

    def __init__(self, parent=None):
        super(SingleWorker, self).__init__(parent)
        self.params_dct = {}
        self.run_count = 0
        self.total_count = 0
        self.info_str = ""

    def request_work(self, params_dct):
        """
        Get parameters from Main window
        :param(dict) params_dct: a dict contains all parameters from UI
        """
        self.workRequested.emit()
        self.params_dct = params_dct

    def infoback(self):

        return self.info_str

    def run_prediction_task(self):

        print("Running Prediction Task...")
        time.sleep(1)
        self.info_str = f"Running Prediction Task..."
        self.infoback()
        self.info_update.emit(self.info_str)
        run_info = get_predictions(self.params_dct)

        if run_info is not False:

            if isinstance(run_info, float):
                time.sleep(1)
                self.info_str = (
                    f"\n>>> >>> >>> FINISHED in {run_info:.3f} Sec <<< <<< <<<\n"
                )
                self.infoback()
                self.info_update.emit(self.info_str)
            else:
                if isinstance(run_info, str):
                    time.sleep(1)
                    self.info_str = f"{run_info}\n>>> >>> >>> FINISHED <<< <<< <<<\n"
                    self.infoback()
                    self.info_update.emit(self.info_str)
                else:
                    time.sleep(1)
                    self.info_str = "\n!! Sorry, an error has occurred, please check your settings !!\n\n"
                    self.infoback()
                    self.info_update.emit(self.info_str)

        else:
            time.sleep(1)
            self.info_str = (
                "\n!! Sorry, an error has occurred, please check your settings !!\n\n"
            )
            self.infoback()
            self.info_update.emit(self.info_str)
            time.sleep(1)
            self.finished.emit()

        time.sleep(1)
        self.finished.emit()

    def run_identification_task(self):

        print("Running Identification Task...")
        time.sleep(1)
        self.info_str = f"Running Identification Task..."
        self.infoback()
        self.info_update.emit(self.info_str)
        run_info = ident(self.params_dct)

        if run_info is not False:

            if isinstance(run_info, float):
                time.sleep(1)
                self.info_str = (
                    f"\n>>> >>> >>> FINISHED in {run_info:.3f} Sec <<< <<< <<<\n"
                )
                self.infoback()
                self.info_update.emit(self.info_str)
            else:
                if isinstance(run_info, str):
                    time.sleep(1)
                    self.info_str = f"{run_info}\n>>> >>> >>> FINISHED <<< <<< <<<\n"
                    self.infoback()
                    self.info_update.emit(self.info_str)
                else:
                    time.sleep(1)
                    self.info_str = "\n!! Sorry, an error has occurred, please check your settings !!\n\n"
                    self.infoback()
                    self.info_update.emit(self.info_str)

        else:
            time.sleep(1)
            self.info_str = (
                "\n!! Sorry, an error has occurred, please check your settings !!\n\n"
            )
            self.infoback()
            self.info_update.emit(self.info_str)
            time.sleep(1)
            self.finished.emit()

        time.sleep(1)
        self.finished.emit()


if __name__ == "__main__":
    import sys

    # # Scale GUI on HiDPI monitors e.g. 2k, 4k resolution
    # QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    # QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps)
    # os.chdir('../..')
    # usr_tiger_folder_path = os.path.abspath(os.getcwd())

    multiprocessing.freeze_support()
    app = QtWidgets.QApplication(sys.argv)
    # window = Tiger(tiger_cwd=usr_tiger_folder_path)
    window = TigerMainWindow()
    window.show()
    app.exec()
