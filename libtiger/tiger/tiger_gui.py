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
################################################################################
## Form generated from reading UI file 'tiger_GUI.ui'
##
## Created by: Qt User Interface Compiler version 6.3.0
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (
    QCoreApplication,
    QDate,
    QDateTime,
    QLocale,
    QMetaObject,
    QObject,
    QPoint,
    QRect,
    QSize,
    QTime,
    QUrl,
    Qt,
)
from PySide6.QtGui import (
    QBrush,
    QColor,
    QConicalGradient,
    QCursor,
    QFont,
    QFontDatabase,
    QGradient,
    QIcon,
    QImage,
    QKeySequence,
    QLinearGradient,
    QPainter,
    QPalette,
    QPixmap,
    QRadialGradient,
    QTransform,
)
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QPushButton,
    QSizePolicy,
    QSpacerItem,
    QSpinBox,
    QStatusBar,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName("MainWindow")
        MainWindow.resize(900, 600)
        MainWindow.setMinimumSize(QSize(900, 600))
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.tabWidget = QTabWidget(self.centralwidget)
        self.tabWidget.setObjectName("tabWidget")
        self.tabWidget.setGeometry(QRect(10, 10, 881, 561))
        self.tab_a = QWidget()
        self.tab_a.setObjectName("tab_a")
        self.verticalLayoutWidget = QWidget(self.tab_a)
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayoutWidget.setGeometry(QRect(10, 10, 861, 401))
        self.verticalLayout = QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout = QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QLabel(self.verticalLayoutWidget)
        self.label.setObjectName("label")

        self.horizontalLayout.addWidget(self.label)

        self.input_fa_le = QLineEdit(self.verticalLayoutWidget)
        self.input_fa_le.setObjectName("input_fa_le")
        self.input_fa_le.setMinimumSize(QSize(0, 28))

        self.horizontalLayout.addWidget(self.input_fa_le)

        self.input_fa_pb = QPushButton(self.verticalLayoutWidget)
        self.input_fa_pb.setObjectName("input_fa_pb")

        self.horizontalLayout.addWidget(self.input_fa_pb)

        self.verticalLayout.addLayout(self.horizontalLayout)

        self.horizontalLayout_8 = QHBoxLayout()
        self.horizontalLayout_8.setObjectName("horizontalLayout_8")
        self.lpa_chkb = QCheckBox(self.verticalLayoutWidget)
        self.lpa_chkb.setObjectName("lpa_chkb")

        self.horizontalLayout_8.addWidget(self.lpa_chkb)

        self.lpc_chkb = QCheckBox(self.verticalLayoutWidget)
        self.lpc_chkb.setObjectName("lpc_chkb")
        self.lpc_chkb.setChecked(False)

        self.horizontalLayout_8.addWidget(self.lpc_chkb)

        self.lpe_chkb = QCheckBox(self.verticalLayoutWidget)
        self.lpe_chkb.setObjectName("lpe_chkb")
        self.lpe_chkb.setChecked(False)

        self.horizontalLayout_8.addWidget(self.lpe_chkb)

        self.lpg_chkb = QCheckBox(self.verticalLayoutWidget)
        self.lpg_chkb.setObjectName("lpg_chkb")

        self.horizontalLayout_8.addWidget(self.lpg_chkb)

        self.lpi_chkb = QCheckBox(self.verticalLayoutWidget)
        self.lpi_chkb.setObjectName("lpi_chkb")

        self.horizontalLayout_8.addWidget(self.lpi_chkb)

        self.lps_chkb = QCheckBox(self.verticalLayoutWidget)
        self.lps_chkb.setObjectName("lps_chkb")

        self.horizontalLayout_8.addWidget(self.lps_chkb)

        self.horizontalSpacer_4 = QSpacerItem(
            40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum
        )

        self.horizontalLayout_8.addItem(self.horizontalSpacer_4)

        self.verticalLayout.addLayout(self.horizontalLayout_8)

        self.horizontalLayout_9 = QHBoxLayout()
        self.horizontalLayout_9.setObjectName("horizontalLayout_9")
        self.pa_chkb = QCheckBox(self.verticalLayoutWidget)
        self.pa_chkb.setObjectName("pa_chkb")

        self.horizontalLayout_9.addWidget(self.pa_chkb)

        self.pc_chkb = QCheckBox(self.verticalLayoutWidget)
        self.pc_chkb.setObjectName("pc_chkb")
        self.pc_chkb.setChecked(True)

        self.horizontalLayout_9.addWidget(self.pc_chkb)

        self.pe_chkb = QCheckBox(self.verticalLayoutWidget)
        self.pe_chkb.setObjectName("pe_chkb")
        self.pe_chkb.setChecked(True)

        self.horizontalLayout_9.addWidget(self.pe_chkb)

        self.pg_chkb = QCheckBox(self.verticalLayoutWidget)
        self.pg_chkb.setObjectName("pg_chkb")

        self.horizontalLayout_9.addWidget(self.pg_chkb)

        self.pi_chkb = QCheckBox(self.verticalLayoutWidget)
        self.pi_chkb.setObjectName("pi_chkb")

        self.horizontalLayout_9.addWidget(self.pi_chkb)

        self.ps_chkb = QCheckBox(self.verticalLayoutWidget)
        self.ps_chkb.setObjectName("ps_chkb")

        self.horizontalLayout_9.addWidget(self.ps_chkb)

        self.horizontalSpacer_3 = QSpacerItem(
            40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum
        )

        self.horizontalLayout_9.addItem(self.horizontalSpacer_3)

        self.verticalLayout.addLayout(self.horizontalLayout_9)

        self.horizontalLayout_10 = QHBoxLayout()
        self.horizontalLayout_10.setObjectName("horizontalLayout_10")
        self.ce_chkb = QCheckBox(self.verticalLayoutWidget)
        self.ce_chkb.setObjectName("ce_chkb")
        self.ce_chkb.setEnabled(True)
        self.ce_chkb.setChecked(False)

        self.horizontalLayout_10.addWidget(self.ce_chkb)

        self.dg_chkb = QCheckBox(self.verticalLayoutWidget)
        self.dg_chkb.setObjectName("dg_chkb")

        self.horizontalLayout_10.addWidget(self.dg_chkb)

        self.tg_chkb = QCheckBox(self.verticalLayoutWidget)
        self.tg_chkb.setObjectName("tg_chkb")
        self.tg_chkb.setChecked(False)

        self.horizontalLayout_10.addWidget(self.tg_chkb)

        self.horizontalSpacer_2 = QSpacerItem(
            40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum
        )

        self.horizontalLayout_10.addItem(self.horizontalSpacer_2)

        self.verticalLayout.addLayout(self.horizontalLayout_10)

        self.line = QFrame(self.verticalLayoutWidget)
        self.line.setObjectName("line")
        self.line.setFrameShape(QFrame.HLine)
        self.line.setFrameShadow(QFrame.Sunken)

        self.verticalLayout.addWidget(self.line)

        self.horizontalLayout_5 = QHBoxLayout()
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.label_6 = QLabel(self.verticalLayoutWidget)
        self.label_6.setObjectName("label_6")

        self.horizontalLayout_5.addWidget(self.label_6)

        self.max_site_spb = QSpinBox(self.verticalLayoutWidget)
        self.max_site_spb.setObjectName("max_site_spb")
        self.max_site_spb.setMaximum(9)
        self.max_site_spb.setValue(2)

        self.horizontalLayout_5.addWidget(self.max_site_spb)

        self.label_7 = QLabel(self.verticalLayoutWidget)
        self.label_7.setObjectName("label_7")

        self.horizontalLayout_5.addWidget(self.label_7)

        self.max_o_spb = QSpinBox(self.verticalLayoutWidget)
        self.max_o_spb.setObjectName("max_o_spb")
        self.max_o_spb.setMaximum(9)
        self.max_o_spb.setValue(3)

        self.horizontalLayout_5.addWidget(self.max_o_spb)

        self.max_oh_lb = QLabel(self.verticalLayoutWidget)
        self.max_oh_lb.setObjectName("max_oh_lb")
        self.max_oh_lb.setAlignment(Qt.AlignRight | Qt.AlignTrailing | Qt.AlignVCenter)

        self.horizontalLayout_5.addWidget(self.max_oh_lb)

        self.max_oh_spb = QSpinBox(self.verticalLayoutWidget)
        self.max_oh_spb.setObjectName("max_oh_spb")
        self.max_oh_spb.setMaximum(9)
        self.max_oh_spb.setValue(2)

        self.horizontalLayout_5.addWidget(self.max_oh_spb)

        self.label_4 = QLabel(self.verticalLayoutWidget)
        self.label_4.setObjectName("label_4")
        self.label_4.setAlignment(Qt.AlignRight | Qt.AlignTrailing | Qt.AlignVCenter)

        self.horizontalLayout_5.addWidget(self.label_4)

        self.max_keto_spb = QSpinBox(self.verticalLayoutWidget)
        self.max_keto_spb.setObjectName("max_keto_spb")
        self.max_keto_spb.setMaximum(9)
        self.max_keto_spb.setValue(1)

        self.horizontalLayout_5.addWidget(self.max_keto_spb)

        self.label_5 = QLabel(self.verticalLayoutWidget)
        self.label_5.setObjectName("label_5")
        self.label_5.setAlignment(Qt.AlignRight | Qt.AlignTrailing | Qt.AlignVCenter)

        self.horizontalLayout_5.addWidget(self.label_5)

        self.max_ooh_spb = QSpinBox(self.verticalLayoutWidget)
        self.max_ooh_spb.setObjectName("max_ooh_spb")
        self.max_ooh_spb.setMaximum(9)
        self.max_ooh_spb.setValue(1)

        self.horizontalLayout_5.addWidget(self.max_ooh_spb)

        self.label_10 = QLabel(self.verticalLayoutWidget)
        self.label_10.setObjectName("label_10")
        self.label_10.setAlignment(Qt.AlignRight | Qt.AlignTrailing | Qt.AlignVCenter)

        self.horizontalLayout_5.addWidget(self.label_10)

        self.max_ep_spb = QSpinBox(self.verticalLayoutWidget)
        self.max_ep_spb.setObjectName("max_ep_spb")
        self.max_ep_spb.setMaximum(9)
        self.max_ep_spb.setSingleStep(1)
        self.max_ep_spb.setValue(0)

        self.horizontalLayout_5.addWidget(self.max_ep_spb)

        self.horizontalSpacer = QSpacerItem(
            40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum
        )

        self.horizontalLayout_5.addItem(self.horizontalSpacer)

        self.verticalLayout.addLayout(self.horizontalLayout_5)

        self.line_4 = QFrame(self.verticalLayoutWidget)
        self.line_4.setObjectName("line_4")
        self.line_4.setFrameShape(QFrame.HLine)
        self.line_4.setFrameShadow(QFrame.Sunken)

        self.verticalLayout.addWidget(self.line_4)

        self.horizontalLayout_11 = QHBoxLayout()
        self.horizontalLayout_11.setObjectName("horizontalLayout_11")
        self.line_9 = QFrame(self.verticalLayoutWidget)
        self.line_9.setObjectName("line_9")
        self.line_9.setFrameShape(QFrame.HLine)
        self.line_9.setFrameShadow(QFrame.Sunken)

        self.horizontalLayout_11.addWidget(self.line_9)

        self.refine_epilipidome_chkb = QCheckBox(self.verticalLayoutWidget)
        self.refine_epilipidome_chkb.setObjectName("refine_epilipidome_chkb")

        self.horizontalLayout_11.addWidget(self.refine_epilipidome_chkb)

        self.label_11 = QLabel(self.verticalLayoutWidget)
        self.label_11.setObjectName("label_11")

        self.horizontalLayout_11.addWidget(self.label_11)

        self.refine_epilipidome_le = QLineEdit(self.verticalLayoutWidget)
        self.refine_epilipidome_le.setObjectName("refine_epilipidome_le")
        self.refine_epilipidome_le.setMinimumSize(QSize(0, 28))

        self.horizontalLayout_11.addWidget(self.refine_epilipidome_le)

        self.refine_epilipidome_pb = QPushButton(self.verticalLayoutWidget)
        self.refine_epilipidome_pb.setObjectName("refine_epilipidome_pb")

        self.horizontalLayout_11.addWidget(self.refine_epilipidome_pb)

        self.verticalLayout.addLayout(self.horizontalLayout_11)

        self.line_10 = QFrame(self.verticalLayoutWidget)
        self.line_10.setObjectName("line_10")
        self.line_10.setFrameShape(QFrame.HLine)
        self.line_10.setFrameShadow(QFrame.Sunken)

        self.verticalLayout.addWidget(self.line_10)

        self.horizontalLayout_3 = QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.line_3 = QFrame(self.verticalLayoutWidget)
        self.line_3.setObjectName("line_3")
        self.line_3.setFrameShape(QFrame.HLine)
        self.line_3.setFrameShadow(QFrame.Sunken)

        self.horizontalLayout_3.addWidget(self.line_3)

        self.label_2 = QLabel(self.verticalLayoutWidget)
        self.label_2.setObjectName("label_2")

        self.horizontalLayout_3.addWidget(self.label_2)

        self.output_epilipidome_le = QLineEdit(self.verticalLayoutWidget)
        self.output_epilipidome_le.setObjectName("output_epilipidome_le")
        self.output_epilipidome_le.setMinimumSize(QSize(0, 28))

        self.horizontalLayout_3.addWidget(self.output_epilipidome_le)

        self.output_epilipidome_pb = QPushButton(self.verticalLayoutWidget)
        self.output_epilipidome_pb.setObjectName("output_epilipidome_pb")

        self.horizontalLayout_3.addWidget(self.output_epilipidome_pb)

        self.verticalLayout.addLayout(self.horizontalLayout_3)

        self.line_2 = QFrame(self.verticalLayoutWidget)
        self.line_2.setObjectName("line_2")
        self.line_2.setFrameShape(QFrame.HLine)
        self.line_2.setFrameShadow(QFrame.Sunken)

        self.verticalLayout.addWidget(self.line_2)

        self.run_prediction_pb = QPushButton(self.verticalLayoutWidget)
        self.run_prediction_pb.setObjectName("run_prediction_pb")

        self.verticalLayout.addWidget(self.run_prediction_pb)

        self.run_prediction_status_te = QTextEdit(self.tab_a)
        self.run_prediction_status_te.setObjectName("run_prediction_status_te")
        self.run_prediction_status_te.setGeometry(QRect(10, 410, 861, 111))
        self.tabWidget.addTab(self.tab_a, "")
        self.tab_b = QWidget()
        self.tab_b.setObjectName("tab_b")
        self.verticalLayoutWidget_3 = QWidget(self.tab_b)
        self.verticalLayoutWidget_3.setObjectName("verticalLayoutWidget_3")
        self.verticalLayoutWidget_3.setGeometry(QRect(10, 10, 851, 198))
        self.verticalLayout_3 = QVBoxLayout(self.verticalLayoutWidget_3)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.verticalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_2 = QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label_9 = QLabel(self.verticalLayoutWidget_3)
        self.label_9.setObjectName("label_9")

        self.horizontalLayout_2.addWidget(self.label_9)

        self.input_lipid_le = QLineEdit(self.verticalLayoutWidget_3)
        self.input_lipid_le.setObjectName("input_lipid_le")
        self.input_lipid_le.setMinimumSize(QSize(0, 28))

        self.horizontalLayout_2.addWidget(self.input_lipid_le)

        self.input_lipid_pb = QPushButton(self.verticalLayoutWidget_3)
        self.input_lipid_pb.setObjectName("input_lipid_pb")

        self.horizontalLayout_2.addWidget(self.input_lipid_pb)

        self.verticalLayout_3.addLayout(self.horizontalLayout_2)

        self.horizontalLayout_6 = QHBoxLayout()
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.line_7 = QFrame(self.verticalLayoutWidget_3)
        self.line_7.setObjectName("line_7")
        self.line_7.setFrameShape(QFrame.HLine)
        self.line_7.setFrameShadow(QFrame.Sunken)

        self.horizontalLayout_6.addWidget(self.line_7)

        self.label_3 = QLabel(self.verticalLayoutWidget_3)
        self.label_3.setObjectName("label_3")

        self.horizontalLayout_6.addWidget(self.label_3)

        self.input_epilipidome_le = QLineEdit(self.verticalLayoutWidget_3)
        self.input_epilipidome_le.setObjectName("input_epilipidome_le")
        self.input_epilipidome_le.setMinimumSize(QSize(0, 28))

        self.horizontalLayout_6.addWidget(self.input_epilipidome_le)

        self.input_epilipidome_pb = QPushButton(self.verticalLayoutWidget_3)
        self.input_epilipidome_pb.setObjectName("input_epilipidome_pb")

        self.horizontalLayout_6.addWidget(self.input_epilipidome_pb)

        self.verticalLayout_3.addLayout(self.horizontalLayout_6)

        self.line_8 = QFrame(self.verticalLayoutWidget_3)
        self.line_8.setObjectName("line_8")
        self.line_8.setFrameShape(QFrame.HLine)
        self.line_8.setFrameShadow(QFrame.Sunken)

        self.verticalLayout_3.addWidget(self.line_8)

        self.horizontalLayout_4 = QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.line_5 = QFrame(self.verticalLayoutWidget_3)
        self.line_5.setObjectName("line_5")
        self.line_5.setFrameShape(QFrame.HLine)
        self.line_5.setFrameShadow(QFrame.Sunken)

        self.horizontalLayout_4.addWidget(self.line_5)

        self.label_8 = QLabel(self.verticalLayoutWidget_3)
        self.label_8.setObjectName("label_8")

        self.horizontalLayout_4.addWidget(self.label_8)

        self.output_incl_list_le = QLineEdit(self.verticalLayoutWidget_3)
        self.output_incl_list_le.setObjectName("output_incl_list_le")
        self.output_incl_list_le.setMinimumSize(QSize(0, 28))

        self.horizontalLayout_4.addWidget(self.output_incl_list_le)

        self.output_incl_list_pb = QPushButton(self.verticalLayoutWidget_3)
        self.output_incl_list_pb.setObjectName("output_incl_list_pb")

        self.horizontalLayout_4.addWidget(self.output_incl_list_pb)

        self.verticalLayout_3.addLayout(self.horizontalLayout_4)

        self.run_incl_list_pb = QPushButton(self.verticalLayoutWidget_3)
        self.run_incl_list_pb.setObjectName("run_incl_list_pb")

        self.verticalLayout_3.addWidget(self.run_incl_list_pb)

        self.run_incl_list_status_te = QTextEdit(self.tab_b)
        self.run_incl_list_status_te.setObjectName("run_incl_list_status_te")
        self.run_incl_list_status_te.setGeometry(QRect(10, 230, 851, 291))
        self.tabWidget.addTab(self.tab_b, "")
        self.tab = QWidget()
        self.tab.setObjectName("tab")
        self.verticalLayoutWidget_4 = QWidget(self.tab)
        self.verticalLayoutWidget_4.setObjectName("verticalLayoutWidget_4")
        self.verticalLayoutWidget_4.setGeometry(QRect(10, 0, 851, 431))
        self.verticalLayout_4 = QVBoxLayout(self.verticalLayoutWidget_4)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.verticalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.Spectra_tab = QTabWidget(self.verticalLayoutWidget_4)
        self.Spectra_tab.setObjectName("Spectra_tab")
        self.Spectra_tab.setSizeIncrement(QSize(0, 0))
        self.tab_2 = QWidget()
        self.tab_2.setObjectName("tab_2")
        self.verticalLayoutWidget_6 = QWidget(self.tab_2)
        self.verticalLayoutWidget_6.setObjectName("verticalLayoutWidget_6")
        self.verticalLayoutWidget_6.setGeometry(QRect(10, 9, 831, 151))
        self.verticalLayout_6 = QVBoxLayout(self.verticalLayoutWidget_6)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.verticalLayout_6.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_13 = QHBoxLayout()
        self.horizontalLayout_13.setObjectName("horizontalLayout_13")
        self.label_17 = QLabel(self.verticalLayoutWidget_6)
        self.label_17.setObjectName("label_17")

        self.horizontalLayout_13.addWidget(self.label_17)

        self.input_spectra_mzml_le = QLineEdit(self.verticalLayoutWidget_6)
        self.input_spectra_mzml_le.setObjectName("input_spectra_mzml_le")
        self.input_spectra_mzml_le.setMinimumSize(QSize(0, 28))

        self.horizontalLayout_13.addWidget(self.input_spectra_mzml_le)

        self.input_spectra_mzml_pb = QPushButton(self.verticalLayoutWidget_6)
        self.input_spectra_mzml_pb.setObjectName("input_spectra_mzml_pb")

        self.horizontalLayout_13.addWidget(self.input_spectra_mzml_pb)

        self.verticalLayout_6.addLayout(self.horizontalLayout_13)

        self.horizontalLayout_16 = QHBoxLayout()
        self.horizontalLayout_16.setObjectName("horizontalLayout_16")
        self.label_14 = QLabel(self.verticalLayoutWidget_6)
        self.label_14.setObjectName("label_14")

        self.horizontalLayout_16.addWidget(self.label_14)

        self.ms1_ppm_spb = QSpinBox(self.verticalLayoutWidget_6)
        self.ms1_ppm_spb.setObjectName("ms1_ppm_spb")
        self.ms1_ppm_spb.setMinimumSize(QSize(0, 24))
        self.ms1_ppm_spb.setMaximum(500)
        self.ms1_ppm_spb.setValue(5)

        self.horizontalLayout_16.addWidget(self.ms1_ppm_spb)

        self.label_15 = QLabel(self.verticalLayoutWidget_6)
        self.label_15.setObjectName("label_15")

        self.horizontalLayout_16.addWidget(self.label_15)

        self.rt_min_spb = QDoubleSpinBox(self.verticalLayoutWidget_6)
        self.rt_min_spb.setObjectName("rt_min_spb")
        self.rt_min_spb.setMinimumSize(QSize(0, 24))

        self.horizontalLayout_16.addWidget(self.rt_min_spb)

        self.max_oh_lb_2 = QLabel(self.verticalLayoutWidget_6)
        self.max_oh_lb_2.setObjectName("max_oh_lb_2")
        self.max_oh_lb_2.setAlignment(
            Qt.AlignRight | Qt.AlignTrailing | Qt.AlignVCenter
        )

        self.horizontalLayout_16.addWidget(self.max_oh_lb_2)

        self.rt_max_spb = QDoubleSpinBox(self.verticalLayoutWidget_6)
        self.rt_max_spb.setObjectName("rt_max_spb")
        self.rt_max_spb.setMinimumSize(QSize(0, 24))
        self.rt_max_spb.setValue(20.000000000000000)

        self.horizontalLayout_16.addWidget(self.rt_max_spb)

        self.label_16 = QLabel(self.verticalLayoutWidget_6)
        self.label_16.setObjectName("label_16")

        self.horizontalLayout_16.addWidget(self.label_16)

        self.ms1_threshold_spb = QSpinBox(self.verticalLayoutWidget_6)
        self.ms1_threshold_spb.setObjectName("ms1_threshold_spb")
        self.ms1_threshold_spb.setMinimumSize(QSize(0, 24))
        self.ms1_threshold_spb.setMaximum(99999999)
        self.ms1_threshold_spb.setSingleStep(1000)
        self.ms1_threshold_spb.setValue(2000)

        self.horizontalLayout_16.addWidget(self.ms1_threshold_spb)

        self.horizontalSpacer_8 = QSpacerItem(
            40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum
        )

        self.horizontalLayout_16.addItem(self.horizontalSpacer_8)

        self.verticalLayout_6.addLayout(self.horizontalLayout_16)

        self.horizontalLayout_20 = QHBoxLayout()
        self.horizontalLayout_20.setObjectName("horizontalLayout_20")
        self.label_19 = QLabel(self.verticalLayoutWidget_6)
        self.label_19.setObjectName("label_19")

        self.horizontalLayout_20.addWidget(self.label_19)

        self.ms2_ppm_spb = QSpinBox(self.verticalLayoutWidget_6)
        self.ms2_ppm_spb.setObjectName("ms2_ppm_spb")
        self.ms2_ppm_spb.setMinimumSize(QSize(0, 24))
        self.ms2_ppm_spb.setMaximum(500)
        self.ms2_ppm_spb.setValue(10)

        self.horizontalLayout_20.addWidget(self.ms2_ppm_spb)

        self.label_18 = QLabel(self.verticalLayoutWidget_6)
        self.label_18.setObjectName("label_18")

        self.horizontalLayout_20.addWidget(self.label_18)

        self.pr_mz_min_spb = QDoubleSpinBox(self.verticalLayoutWidget_6)
        self.pr_mz_min_spb.setObjectName("pr_mz_min_spb")
        self.pr_mz_min_spb.setMinimumSize(QSize(0, 24))
        self.pr_mz_min_spb.setDecimals(2)
        self.pr_mz_min_spb.setMaximum(2000.000000000000000)
        self.pr_mz_min_spb.setSingleStep(50.000000000000000)
        self.pr_mz_min_spb.setValue(400.000000000000000)

        self.horizontalLayout_20.addWidget(self.pr_mz_min_spb)

        self.max_oh_lb_3 = QLabel(self.verticalLayoutWidget_6)
        self.max_oh_lb_3.setObjectName("max_oh_lb_3")
        self.max_oh_lb_3.setAlignment(
            Qt.AlignRight | Qt.AlignTrailing | Qt.AlignVCenter
        )

        self.horizontalLayout_20.addWidget(self.max_oh_lb_3)

        self.pr_mz_max_spb = QDoubleSpinBox(self.verticalLayoutWidget_6)
        self.pr_mz_max_spb.setObjectName("pr_mz_max_spb")
        self.pr_mz_max_spb.setMinimumSize(QSize(0, 24))
        self.pr_mz_max_spb.setDecimals(2)
        self.pr_mz_max_spb.setMaximum(2000.000000000000000)
        self.pr_mz_max_spb.setSingleStep(50.000000000000000)
        self.pr_mz_max_spb.setValue(1200.000000000000000)

        self.horizontalLayout_20.addWidget(self.pr_mz_max_spb)

        self.label_20 = QLabel(self.verticalLayoutWidget_6)
        self.label_20.setObjectName("label_20")

        self.horizontalLayout_20.addWidget(self.label_20)

        self.selection_window_spb = QDoubleSpinBox(self.verticalLayoutWidget_6)
        self.selection_window_spb.setObjectName("selection_window_spb")
        self.selection_window_spb.setMinimumSize(QSize(0, 24))
        self.selection_window_spb.setDecimals(2)
        self.selection_window_spb.setMaximum(5.000000000000000)
        self.selection_window_spb.setSingleStep(0.100000000000000)
        self.selection_window_spb.setValue(0.500000000000000)

        self.horizontalLayout_20.addWidget(self.selection_window_spb)

        self.horizontalSpacer_9 = QSpacerItem(
            40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum
        )

        self.horizontalLayout_20.addItem(self.horizontalSpacer_9)

        self.verticalLayout_6.addLayout(self.horizontalLayout_20)

        self.horizontalLayout_21 = QHBoxLayout()
        self.horizontalLayout_21.setObjectName("horizontalLayout_21")
        self.label_27 = QLabel(self.verticalLayoutWidget_6)
        self.label_27.setObjectName("label_27")

        self.horizontalLayout_21.addWidget(self.label_27)

        self.max_oh_lb_4 = QLabel(self.verticalLayoutWidget_6)
        self.max_oh_lb_4.setObjectName("max_oh_lb_4")
        self.max_oh_lb_4.setAlignment(
            Qt.AlignRight | Qt.AlignTrailing | Qt.AlignVCenter
        )

        self.horizontalLayout_21.addWidget(self.max_oh_lb_4)

        self.total_score_min_spb = QSpinBox(self.verticalLayoutWidget_6)
        self.total_score_min_spb.setObjectName("total_score_min_spb")
        self.total_score_min_spb.setMaximum(99)
        self.total_score_min_spb.setValue(50)

        self.horizontalLayout_21.addWidget(self.total_score_min_spb)

        self.label_32 = QLabel(self.verticalLayoutWidget_6)
        self.label_32.setObjectName("label_32")

        self.horizontalLayout_21.addWidget(self.label_32)

        self.isotope_score_min_spb = QSpinBox(self.verticalLayoutWidget_6)
        self.isotope_score_min_spb.setObjectName("isotope_score_min_spb")
        self.isotope_score_min_spb.setMaximum(99)
        self.isotope_score_min_spb.setValue(80)

        self.horizontalLayout_21.addWidget(self.isotope_score_min_spb)

        self.label_28 = QLabel(self.verticalLayoutWidget_6)
        self.label_28.setObjectName("label_28")

        self.horizontalLayout_21.addWidget(self.label_28)

        self.rank_score_min_spb = QSpinBox(self.verticalLayoutWidget_6)
        self.rank_score_min_spb.setObjectName("rank_score_min_spb")
        self.rank_score_min_spb.setMaximum(99)
        self.rank_score_min_spb.setValue(40)

        self.horizontalLayout_21.addWidget(self.rank_score_min_spb)

        self.label_23 = QLabel(self.verticalLayoutWidget_6)
        self.label_23.setObjectName("label_23")

        self.horizontalLayout_21.addWidget(self.label_23)

        self.ms2_threshold_spb = QDoubleSpinBox(self.verticalLayoutWidget_6)
        self.ms2_threshold_spb.setObjectName("ms2_threshold_spb")
        self.ms2_threshold_spb.setMinimumSize(QSize(0, 24))
        self.ms2_threshold_spb.setDecimals(1)
        self.ms2_threshold_spb.setMaximum(100.000000000000000)
        self.ms2_threshold_spb.setSingleStep(0.500000000000000)
        self.ms2_threshold_spb.setValue(1.000000000000000)

        self.horizontalLayout_21.addWidget(self.ms2_threshold_spb)

        self.horizontalSpacer_10 = QSpacerItem(
            40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum
        )

        self.horizontalLayout_21.addItem(self.horizontalSpacer_10)

        self.verticalLayout_6.addLayout(self.horizontalLayout_21)

        self.Spectra_tab.addTab(self.tab_2, "")
        self.tab_3 = QWidget()
        self.tab_3.setObjectName("tab_3")
        self.verticalLayoutWidget_7 = QWidget(self.tab_3)
        self.verticalLayoutWidget_7.setObjectName("verticalLayoutWidget_7")
        self.verticalLayoutWidget_7.setGeometry(QRect(10, 10, 831, 151))
        self.verticalLayout_7 = QVBoxLayout(self.verticalLayoutWidget_7)
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.verticalLayout_7.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_14 = QHBoxLayout()
        self.horizontalLayout_14.setObjectName("horizontalLayout_14")
        self.label_25 = QLabel(self.verticalLayoutWidget_7)
        self.label_25.setObjectName("label_25")

        self.horizontalLayout_14.addWidget(self.label_25)

        self.input_spectra_star_le = QLineEdit(self.verticalLayoutWidget_7)
        self.input_spectra_star_le.setObjectName("input_spectra_star_le")
        self.input_spectra_star_le.setMinimumSize(QSize(0, 28))

        self.horizontalLayout_14.addWidget(self.input_spectra_star_le)

        self.input_spectra_star_pb = QPushButton(self.verticalLayoutWidget_7)
        self.input_spectra_star_pb.setObjectName("input_spectra_star_pb")

        self.horizontalLayout_14.addWidget(self.input_spectra_star_pb)

        self.verticalLayout_7.addLayout(self.horizontalLayout_14)

        self.horizontalLayout_22 = QHBoxLayout()
        self.horizontalLayout_22.setObjectName("horizontalLayout_22")
        self.label_36 = QLabel(self.verticalLayoutWidget_7)
        self.label_36.setObjectName("label_36")

        self.horizontalLayout_22.addWidget(self.label_36)

        self.star_ms2_ppm_spb = QSpinBox(self.verticalLayoutWidget_7)
        self.star_ms2_ppm_spb.setObjectName("star_ms2_ppm_spb")
        self.star_ms2_ppm_spb.setMinimumSize(QSize(0, 24))
        self.star_ms2_ppm_spb.setMaximum(500)
        self.star_ms2_ppm_spb.setValue(10)

        self.horizontalLayout_22.addWidget(self.star_ms2_ppm_spb)

        self.label_37 = QLabel(self.verticalLayoutWidget_7)
        self.label_37.setObjectName("label_37")

        self.horizontalLayout_22.addWidget(self.label_37)

        self.star_pr_mz_min_spb = QDoubleSpinBox(self.verticalLayoutWidget_7)
        self.star_pr_mz_min_spb.setObjectName("star_pr_mz_min_spb")
        self.star_pr_mz_min_spb.setMinimumSize(QSize(0, 24))
        self.star_pr_mz_min_spb.setDecimals(2)
        self.star_pr_mz_min_spb.setMaximum(2000.000000000000000)
        self.star_pr_mz_min_spb.setSingleStep(50.000000000000000)
        self.star_pr_mz_min_spb.setValue(400.000000000000000)

        self.horizontalLayout_22.addWidget(self.star_pr_mz_min_spb)

        self.max_oh_lb_5 = QLabel(self.verticalLayoutWidget_7)
        self.max_oh_lb_5.setObjectName("max_oh_lb_5")
        self.max_oh_lb_5.setAlignment(
            Qt.AlignRight | Qt.AlignTrailing | Qt.AlignVCenter
        )

        self.horizontalLayout_22.addWidget(self.max_oh_lb_5)

        self.star_pr_mz_max_spb = QDoubleSpinBox(self.verticalLayoutWidget_7)
        self.star_pr_mz_max_spb.setObjectName("star_pr_mz_max_spb")
        self.star_pr_mz_max_spb.setMinimumSize(QSize(0, 24))
        self.star_pr_mz_max_spb.setDecimals(2)
        self.star_pr_mz_max_spb.setMaximum(2000.000000000000000)
        self.star_pr_mz_max_spb.setSingleStep(50.000000000000000)
        self.star_pr_mz_max_spb.setValue(1200.000000000000000)

        self.horizontalLayout_22.addWidget(self.star_pr_mz_max_spb)

        self.horizontalSpacer_11 = QSpacerItem(
            40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum
        )

        self.horizontalLayout_22.addItem(self.horizontalSpacer_11)

        self.verticalLayout_7.addLayout(self.horizontalLayout_22)

        self.horizontalLayout_23 = QHBoxLayout()
        self.horizontalLayout_23.setObjectName("horizontalLayout_23")
        self.label_39 = QLabel(self.verticalLayoutWidget_7)
        self.label_39.setObjectName("label_39")

        self.horizontalLayout_23.addWidget(self.label_39)

        self.max_oh_lb_6 = QLabel(self.verticalLayoutWidget_7)
        self.max_oh_lb_6.setObjectName("max_oh_lb_6")
        self.max_oh_lb_6.setAlignment(
            Qt.AlignRight | Qt.AlignTrailing | Qt.AlignVCenter
        )

        self.horizontalLayout_23.addWidget(self.max_oh_lb_6)

        self.star_total_score_min_spb = QSpinBox(self.verticalLayoutWidget_7)
        self.star_total_score_min_spb.setObjectName("star_total_score_min_spb")
        self.star_total_score_min_spb.setMaximum(99)
        self.star_total_score_min_spb.setValue(50)

        self.horizontalLayout_23.addWidget(self.star_total_score_min_spb)

        self.label_41 = QLabel(self.verticalLayoutWidget_7)
        self.label_41.setObjectName("label_41")

        self.horizontalLayout_23.addWidget(self.label_41)

        self.star_rank_score_min_spb = QSpinBox(self.verticalLayoutWidget_7)
        self.star_rank_score_min_spb.setObjectName("star_rank_score_min_spb")
        self.star_rank_score_min_spb.setMaximum(99)
        self.star_rank_score_min_spb.setValue(40)

        self.horizontalLayout_23.addWidget(self.star_rank_score_min_spb)

        self.horizontalSpacer_12 = QSpacerItem(
            40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum
        )

        self.horizontalLayout_23.addItem(self.horizontalSpacer_12)

        self.verticalLayout_7.addLayout(self.horizontalLayout_23)

        self.verticalSpacer = QSpacerItem(
            20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding
        )

        self.verticalLayout_7.addItem(self.verticalSpacer)

        self.Spectra_tab.addTab(self.tab_3, "")

        self.verticalLayout_4.addWidget(self.Spectra_tab)

        self.horizontalLayout_12 = QHBoxLayout()
        self.horizontalLayout_12.setObjectName("horizontalLayout_12")
        self.label_13 = QLabel(self.verticalLayoutWidget_4)
        self.label_13.setObjectName("label_13")

        self.horizontalLayout_12.addWidget(self.label_13)

        self.input_ident_epilipidome_le = QLineEdit(self.verticalLayoutWidget_4)
        self.input_ident_epilipidome_le.setObjectName("input_ident_epilipidome_le")
        self.input_ident_epilipidome_le.setMinimumSize(QSize(0, 28))

        self.horizontalLayout_12.addWidget(self.input_ident_epilipidome_le)

        self.input_ident_epilipidome_pb = QPushButton(self.verticalLayoutWidget_4)
        self.input_ident_epilipidome_pb.setObjectName("input_ident_epilipidome_pb")

        self.horizontalLayout_12.addWidget(self.input_ident_epilipidome_pb)

        self.verticalLayout_4.addLayout(self.horizontalLayout_12)

        self.line_12 = QFrame(self.verticalLayoutWidget_4)
        self.line_12.setObjectName("line_12")
        self.line_12.setFrameShape(QFrame.HLine)
        self.line_12.setFrameShadow(QFrame.Sunken)

        self.verticalLayout_4.addWidget(self.line_12)

        self.horizontalLayout_17 = QHBoxLayout()
        self.horizontalLayout_17.setObjectName("horizontalLayout_17")
        self.line_13 = QFrame(self.verticalLayoutWidget_4)
        self.line_13.setObjectName("line_13")
        self.line_13.setFrameShape(QFrame.HLine)
        self.line_13.setFrameShadow(QFrame.Sunken)

        self.horizontalLayout_17.addWidget(self.line_13)

        self.label_21 = QLabel(self.verticalLayoutWidget_4)
        self.label_21.setObjectName("label_21")

        self.horizontalLayout_17.addWidget(self.label_21)

        self.report_folder_le = QLineEdit(self.verticalLayoutWidget_4)
        self.report_folder_le.setObjectName("report_folder_le")
        self.report_folder_le.setMinimumSize(QSize(0, 28))

        self.horizontalLayout_17.addWidget(self.report_folder_le)

        self.report_folder_pb = QPushButton(self.verticalLayoutWidget_4)
        self.report_folder_pb.setObjectName("report_folder_pb")

        self.horizontalLayout_17.addWidget(self.report_folder_pb)

        self.label_24 = QLabel(self.verticalLayoutWidget_4)
        self.label_24.setObjectName("label_24")

        self.horizontalLayout_17.addWidget(self.label_24)

        self.dpi_cmb = QComboBox(self.verticalLayoutWidget_4)
        self.dpi_cmb.addItem("")
        self.dpi_cmb.addItem("")
        self.dpi_cmb.addItem("")
        self.dpi_cmb.setObjectName("dpi_cmb")
        self.dpi_cmb.setMinimumSize(QSize(150, 24))

        self.horizontalLayout_17.addWidget(self.dpi_cmb)

        self.verticalLayout_4.addLayout(self.horizontalLayout_17)

        self.horizontalLayout_19 = QHBoxLayout()
        self.horizontalLayout_19.setObjectName("horizontalLayout_19")
        self.line_17 = QFrame(self.verticalLayoutWidget_4)
        self.line_17.setObjectName("line_17")
        self.line_17.setFrameShape(QFrame.HLine)
        self.line_17.setFrameShadow(QFrame.Sunken)

        self.horizontalLayout_19.addWidget(self.line_17)

        self.label_22 = QLabel(self.verticalLayoutWidget_4)
        self.label_22.setObjectName("label_22")

        self.horizontalLayout_19.addWidget(self.label_22)

        self.export_table_le = QLineEdit(self.verticalLayoutWidget_4)
        self.export_table_le.setObjectName("export_table_le")
        self.export_table_le.setMinimumSize(QSize(0, 28))

        self.horizontalLayout_19.addWidget(self.export_table_le)

        self.export_table_pb = QPushButton(self.verticalLayoutWidget_4)
        self.export_table_pb.setObjectName("export_table_pb")

        self.horizontalLayout_19.addWidget(self.export_table_pb)

        self.dump_json_chkb = QCheckBox(self.verticalLayoutWidget_4)
        self.dump_json_chkb.setObjectName("dump_json_chkb")
        self.dump_json_chkb.setChecked(True)

        self.horizontalLayout_19.addWidget(self.dump_json_chkb)

        self.verticalLayout_4.addLayout(self.horizontalLayout_19)

        self.line_14 = QFrame(self.verticalLayoutWidget_4)
        self.line_14.setObjectName("line_14")
        self.line_14.setFrameShape(QFrame.HLine)
        self.line_14.setFrameShadow(QFrame.Sunken)

        self.verticalLayout_4.addWidget(self.line_14)

        self.run_identification_pb = QPushButton(self.verticalLayoutWidget_4)
        self.run_identification_pb.setObjectName("run_identification_pb")

        self.verticalLayout_4.addWidget(self.run_identification_pb)

        self.run_identification_status_te = QTextEdit(self.tab)
        self.run_identification_status_te.setObjectName("run_identification_status_te")
        self.run_identification_status_te.setGeometry(QRect(10, 430, 851, 91))
        self.tabWidget.addTab(self.tab, "")
        self.tab_c = QWidget()
        self.tab_c.setObjectName("tab_c")
        self.verticalLayoutWidget_2 = QWidget(self.tab_c)
        self.verticalLayoutWidget_2.setObjectName("verticalLayoutWidget_2")
        self.verticalLayoutWidget_2.setGeometry(QRect(10, 10, 851, 91))
        self.verticalLayout_2 = QVBoxLayout(self.verticalLayoutWidget_2)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_7 = QHBoxLayout()
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.label_12 = QLabel(self.verticalLayoutWidget_2)
        self.label_12.setObjectName("label_12")

        self.horizontalLayout_7.addWidget(self.label_12)

        self.lineEdit_4 = QLineEdit(self.verticalLayoutWidget_2)
        self.lineEdit_4.setObjectName("lineEdit_4")

        self.horizontalLayout_7.addWidget(self.lineEdit_4)

        self.pushButton_5 = QPushButton(self.verticalLayoutWidget_2)
        self.pushButton_5.setObjectName("pushButton_5")

        self.horizontalLayout_7.addWidget(self.pushButton_5)

        self.verticalLayout_2.addLayout(self.horizontalLayout_7)

        self.line_6 = QFrame(self.verticalLayoutWidget_2)
        self.line_6.setObjectName("line_6")
        self.line_6.setFrameShape(QFrame.HLine)
        self.line_6.setFrameShadow(QFrame.Sunken)

        self.verticalLayout_2.addWidget(self.line_6)

        self.pushButton_8 = QPushButton(self.verticalLayoutWidget_2)
        self.pushButton_8.setObjectName("pushButton_8")

        self.verticalLayout_2.addWidget(self.pushButton_8)

        self.tabWidget.addTab(self.tab_c, "")
        self.tab_z = QWidget()
        self.tab_z.setObjectName("tab_z")
        self.verticalLayoutWidget_5 = QWidget(self.tab_z)
        self.verticalLayoutWidget_5.setObjectName("verticalLayoutWidget_5")
        self.verticalLayoutWidget_5.setGeometry(QRect(20, 10, 852, 501))
        self.verticalLayout_5 = QVBoxLayout(self.verticalLayoutWidget_5)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.verticalLayout_5.setContentsMargins(0, 0, 0, 0)
        self.label_33 = QLabel(self.verticalLayoutWidget_5)
        self.label_33.setObjectName("label_33")

        self.verticalLayout_5.addWidget(self.label_33)

        self.horizontalLayout_15 = QHBoxLayout()
        self.horizontalLayout_15.setObjectName("horizontalLayout_15")
        self.label_29 = QLabel(self.verticalLayoutWidget_5)
        self.label_29.setObjectName("label_29")

        self.horizontalLayout_15.addWidget(self.label_29)

        self.horizontalSpacer_6 = QSpacerItem(
            40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum
        )

        self.horizontalLayout_15.addItem(self.horizontalSpacer_6)

        self.label_31 = QLabel(self.verticalLayoutWidget_5)
        self.label_31.setObjectName("label_31")

        self.horizontalLayout_15.addWidget(self.label_31)

        self.horizontalSpacer_5 = QSpacerItem(
            40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum
        )

        self.horizontalLayout_15.addItem(self.horizontalSpacer_5)

        self.verticalLayout_5.addLayout(self.horizontalLayout_15)

        self.label_26 = QLabel(self.verticalLayoutWidget_5)
        self.label_26.setObjectName("label_26")

        self.verticalLayout_5.addWidget(self.label_26)

        self.label_30 = QLabel(self.verticalLayoutWidget_5)
        self.label_30.setObjectName("label_30")

        self.verticalLayout_5.addWidget(self.label_30)

        self.label_35 = QLabel(self.verticalLayoutWidget_5)
        self.label_35.setObjectName("label_35")

        self.verticalLayout_5.addWidget(self.label_35)

        self.tabWidget.addTab(self.tab_z, "")
        self.label_34 = QLabel(self.centralwidget)
        self.label_34.setObjectName("label_34")
        self.label_34.setGeometry(QRect(20, 575, 861, 21))
        MainWindow.setCentralWidget(self.centralwidget)

        self.retranslateUi(MainWindow)

        self.tabWidget.setCurrentIndex(0)
        self.Spectra_tab.setCurrentIndex(0)
        self.dpi_cmb.setCurrentIndex(1)

        QMetaObject.connectSlotsByName(MainWindow)

    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(
            QCoreApplication.translate(
                "MainWindow", "LPPtiger 2 (Beta) | LMAI@TU Dresden", None
            )
        )
        self.label.setText(
            QCoreApplication.translate("MainWindow", "List of FA:", None)
        )
        self.input_fa_pb.setText(QCoreApplication.translate("MainWindow", "Open", None))
        # if QT_CONFIG(tooltip)
        self.lpa_chkb.setToolTip(
            QCoreApplication.translate(
                "MainWindow", "<html><head/><body><p>LPA [M-H]-</p></body></html>", None
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.lpa_chkb.setText(QCoreApplication.translate("MainWindow", "LPA", None))
        # if QT_CONFIG(tooltip)
        self.lpc_chkb.setToolTip(
            QCoreApplication.translate(
                "MainWindow",
                "<html><head/><body><p>LPC [M+HCOO]-</p></body></html>",
                None,
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.lpc_chkb.setText(QCoreApplication.translate("MainWindow", "LPC", None))
        # if QT_CONFIG(tooltip)
        self.lpe_chkb.setToolTip(
            QCoreApplication.translate(
                "MainWindow", "<html><head/><body><p>LPE [M-H]-</p></body></html>", None
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.lpe_chkb.setText(QCoreApplication.translate("MainWindow", "LPE", None))
        # if QT_CONFIG(tooltip)
        self.lpg_chkb.setToolTip(
            QCoreApplication.translate(
                "MainWindow", "<html><head/><body><p>LPG [M-H]-</p></body></html>", None
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.lpg_chkb.setText(QCoreApplication.translate("MainWindow", "LPG", None))
        # if QT_CONFIG(tooltip)
        self.lpi_chkb.setToolTip(
            QCoreApplication.translate(
                "MainWindow", "<html><head/><body><p>LPI [M-H]-</p></body></html>", None
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.lpi_chkb.setText(QCoreApplication.translate("MainWindow", "LPI", None))
        # if QT_CONFIG(tooltip)
        self.lps_chkb.setToolTip(
            QCoreApplication.translate(
                "MainWindow", "<html><head/><body><p>LPS [M-H]-</p></body></html>", None
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.lps_chkb.setText(QCoreApplication.translate("MainWindow", "LPS", None))
        # if QT_CONFIG(tooltip)
        self.pa_chkb.setToolTip(
            QCoreApplication.translate(
                "MainWindow", "<html><head/><body><p>PA [M-H]-</p></body></html>", None
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.pa_chkb.setText(QCoreApplication.translate("MainWindow", "PA", None))
        # if QT_CONFIG(tooltip)
        self.pc_chkb.setToolTip(
            QCoreApplication.translate(
                "MainWindow",
                "<html><head/><body><p>PC [M+HCOO]-</p></body></html>",
                None,
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.pc_chkb.setText(QCoreApplication.translate("MainWindow", "PC", None))
        # if QT_CONFIG(tooltip)
        self.pe_chkb.setToolTip(
            QCoreApplication.translate(
                "MainWindow", "<html><head/><body><p>PE [M-H]-</p></body></html>", None
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.pe_chkb.setText(QCoreApplication.translate("MainWindow", "PE", None))
        # if QT_CONFIG(tooltip)
        self.pg_chkb.setToolTip(
            QCoreApplication.translate(
                "MainWindow", "<html><head/><body><p>PG [M-H]-</p></body></html>", None
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.pg_chkb.setText(QCoreApplication.translate("MainWindow", "PG", None))
        # if QT_CONFIG(tooltip)
        self.pi_chkb.setToolTip(
            QCoreApplication.translate(
                "MainWindow", "<html><head/><body><p>PI [M-H]-</p></body></html>", None
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.pi_chkb.setText(QCoreApplication.translate("MainWindow", "PI", None))
        # if QT_CONFIG(tooltip)
        self.ps_chkb.setToolTip(
            QCoreApplication.translate(
                "MainWindow", "<html><head/><body><p>PS [M-H]-</p></body></html>", None
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.ps_chkb.setText(QCoreApplication.translate("MainWindow", "PS", None))
        # if QT_CONFIG(tooltip)
        self.ce_chkb.setToolTip(
            QCoreApplication.translate(
                "MainWindow", "<html><head/><body><p>CE [M+Na]+</p></body></html>", None
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.ce_chkb.setText(QCoreApplication.translate("MainWindow", "CE", None))
        # if QT_CONFIG(tooltip)
        self.dg_chkb.setToolTip(
            QCoreApplication.translate(
                "MainWindow", "<html><head/><body><p>DG [M+Na]+</p></body></html>", None
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.dg_chkb.setText(QCoreApplication.translate("MainWindow", "DG", None))
        # if QT_CONFIG(tooltip)
        self.tg_chkb.setToolTip(
            QCoreApplication.translate(
                "MainWindow", "<html><head/><body><p>TG [M+Na]+</p></body></html>", None
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.tg_chkb.setText(QCoreApplication.translate("MainWindow", "TG", None))
        self.label_6.setText(
            QCoreApplication.translate("MainWindow", "Max modification site:", None)
        )
        self.label_7.setText(
            QCoreApplication.translate("MainWindow", "Max total O:", None)
        )
        self.max_oh_lb.setText(QCoreApplication.translate("MainWindow", "max OH", None))
        self.label_4.setText(QCoreApplication.translate("MainWindow", "max keto", None))
        self.label_5.setText(QCoreApplication.translate("MainWindow", "max OOH", None))
        self.label_10.setText(
            QCoreApplication.translate("MainWindow", "max epoxy", None)
        )
        # if QT_CONFIG(tooltip)
        self.refine_epilipidome_chkb.setToolTip(
            QCoreApplication.translate(
                "MainWindow",
                "<html><head/><body><p>Enable this option to submit a list of unmodified lipids to extract predicted epilipidome from given precursors. </p><p>Skip this option to obtain the whole epilipidome predicted using all FA and lipid classes from settings above.</p></body></html>",
                None,
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.refine_epilipidome_chkb.setText(
            QCoreApplication.translate(
                "MainWindow", "Refine epilipidome by unmodified lipids:", None
            )
        )
        self.label_11.setText("")
        # if QT_CONFIG(tooltip)
        self.refine_epilipidome_le.setToolTip(
            QCoreApplication.translate(
                "MainWindow",
                "<html><head/><body><p>Enable this option to submit a list of unmodified lipids to extract predicted epilipidome from given precursors. </p><p>Skip this option to obtain the whole epilipidome predicted using all FA and lipid classes from settings above.</p></body></html>",
                None,
            )
        )
        # endif // QT_CONFIG(tooltip)
        # if QT_CONFIG(tooltip)
        self.refine_epilipidome_pb.setToolTip(
            QCoreApplication.translate(
                "MainWindow",
                "<html><head/><body><p>Enable this option to submit a list of unmodified lipids to extract predicted epilipidome from given precursors. </p><p>Skip this option to obtain the whole epilipidome predicted using all FA and lipid classes from settings above.</p></body></html>",
                None,
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.refine_epilipidome_pb.setText(
            QCoreApplication.translate("MainWindow", "Choose path", None)
        )
        self.label_2.setText(
            QCoreApplication.translate("MainWindow", "Export epilipidome ", None)
        )
        self.output_epilipidome_pb.setText(
            QCoreApplication.translate("MainWindow", "Choose path", None)
        )
        self.run_prediction_pb.setText(
            QCoreApplication.translate("MainWindow", "Run prediction", None)
        )
        self.tabWidget.setTabText(
            self.tabWidget.indexOf(self.tab_a),
            QCoreApplication.translate("MainWindow", "in silico prediction", None),
        )
        self.label_9.setText(
            QCoreApplication.translate("MainWindow", "List of unmodified lipids:", None)
        )
        self.input_lipid_pb.setText(
            QCoreApplication.translate("MainWindow", "Open", None)
        )
        self.label_3.setText(
            QCoreApplication.translate("MainWindow", "Load epilipidome ", None)
        )
        self.input_epilipidome_pb.setText(
            QCoreApplication.translate("MainWindow", "Open", None)
        )
        self.label_8.setText(
            QCoreApplication.translate("MainWindow", "Export inclusion list", None)
        )
        self.output_incl_list_pb.setText(
            QCoreApplication.translate("MainWindow", "Choose path", None)
        )
        self.run_incl_list_pb.setText(
            QCoreApplication.translate("MainWindow", "Generate inclusion list", None)
        )
        self.tabWidget.setTabText(
            self.tabWidget.indexOf(self.tab_b),
            QCoreApplication.translate("MainWindow", "Inclusion list generation", None),
        )
        # if QT_CONFIG(tooltip)
        self.Spectra_tab.setToolTip(
            QCoreApplication.translate(
                "MainWindow",
                "<html><head/><body><p>LPPtiger can read spectra file from .mzML format and can load the exported .csv table from Lipostar 2.</p></body></html>",
                None,
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.label_17.setText(
            QCoreApplication.translate("MainWindow", "Spectra (.mzML)", None)
        )
        self.input_spectra_mzml_pb.setText(
            QCoreApplication.translate("MainWindow", "Open", None)
        )
        self.label_14.setText(
            QCoreApplication.translate("MainWindow", "MS1 tolerance", None)
        )
        self.ms1_ppm_spb.setSuffix(
            QCoreApplication.translate("MainWindow", " ppm", None)
        )
        self.ms1_ppm_spb.setPrefix(
            QCoreApplication.translate("MainWindow", "+/- ", None)
        )
        self.label_15.setText(QCoreApplication.translate("MainWindow", "RT from", None))
        self.rt_min_spb.setSuffix(
            QCoreApplication.translate("MainWindow", " min", None)
        )
        self.max_oh_lb_2.setText(QCoreApplication.translate("MainWindow", "to", None))
        self.rt_max_spb.setSuffix(
            QCoreApplication.translate("MainWindow", " min", None)
        )
        self.label_16.setText(
            QCoreApplication.translate("MainWindow", "Precursor intensity \u2265", None)
        )
        self.ms1_threshold_spb.setSuffix("")
        self.ms1_threshold_spb.setPrefix("")
        self.label_19.setText(
            QCoreApplication.translate("MainWindow", "MS2 tolerance", None)
        )
        self.ms2_ppm_spb.setSuffix(
            QCoreApplication.translate("MainWindow", " ppm", None)
        )
        self.ms2_ppm_spb.setPrefix(
            QCoreApplication.translate("MainWindow", "+/- ", None)
        )
        self.label_18.setText(
            QCoreApplication.translate("MainWindow", "Precursor m/z from", None)
        )
        self.pr_mz_min_spb.setPrefix(
            QCoreApplication.translate("MainWindow", "m/z ", None)
        )
        self.pr_mz_min_spb.setSuffix("")
        self.max_oh_lb_3.setText(QCoreApplication.translate("MainWindow", "to", None))
        self.pr_mz_max_spb.setPrefix(
            QCoreApplication.translate("MainWindow", "m/z ", None)
        )
        self.pr_mz_max_spb.setSuffix("")
        self.label_20.setText(
            QCoreApplication.translate("MainWindow", "Selection window", None)
        )
        self.selection_window_spb.setPrefix(
            QCoreApplication.translate("MainWindow", "m/z +/- ", None)
        )
        self.selection_window_spb.setSuffix("")
        self.label_27.setText(
            QCoreApplication.translate("MainWindow", "Score filters:", None)
        )
        self.max_oh_lb_4.setText(
            QCoreApplication.translate("MainWindow", "Total score \u2265", None)
        )
        self.label_32.setText(
            QCoreApplication.translate("MainWindow", "Isotope score \u2265", None)
        )
        self.label_28.setText(
            QCoreApplication.translate("MainWindow", "Rank score \u2265", None)
        )
        # if QT_CONFIG(tooltip)
        self.label_23.setToolTip(
            QCoreApplication.translate(
                "MainWindow",
                "<html><head/><body><p>Only consider peaks with relative intensity above threshold for score generation. Signals lower than threshold will still be plotted but will not be used to generate scores.</p></body></html>",
                None,
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.label_23.setText(
            QCoreApplication.translate(
                "MainWindow", "Consider peaks with intensity \u2265", None
            )
        )
        # if QT_CONFIG(tooltip)
        self.ms2_threshold_spb.setToolTip(
            QCoreApplication.translate(
                "MainWindow",
                "<html><head/><body><p>Only consider peaks with relative intensity above threshold for score generation. Signals lower than threshold will still be plotted but will not be used to generate scores.</p></body></html>",
                None,
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.ms2_threshold_spb.setPrefix("")
        self.ms2_threshold_spb.setSuffix(
            QCoreApplication.translate("MainWindow", " %", None)
        )
        self.Spectra_tab.setTabText(
            self.Spectra_tab.indexOf(self.tab_2),
            QCoreApplication.translate("MainWindow", "mzML", None),
        )
        self.label_25.setText(
            QCoreApplication.translate("MainWindow", "Lipostar 2 (.csv)", None)
        )
        self.input_spectra_star_pb.setText(
            QCoreApplication.translate("MainWindow", "Open", None)
        )
        self.label_36.setText(
            QCoreApplication.translate("MainWindow", "MS2 tolerance", None)
        )
        self.star_ms2_ppm_spb.setSuffix(
            QCoreApplication.translate("MainWindow", " ppm", None)
        )
        self.star_ms2_ppm_spb.setPrefix(
            QCoreApplication.translate("MainWindow", "+/- ", None)
        )
        self.label_37.setText(
            QCoreApplication.translate("MainWindow", "Precursor m/z from", None)
        )
        self.star_pr_mz_min_spb.setPrefix(
            QCoreApplication.translate("MainWindow", "m/z ", None)
        )
        self.star_pr_mz_min_spb.setSuffix("")
        self.max_oh_lb_5.setText(QCoreApplication.translate("MainWindow", "to", None))
        self.star_pr_mz_max_spb.setPrefix(
            QCoreApplication.translate("MainWindow", "m/z ", None)
        )
        self.star_pr_mz_max_spb.setSuffix("")
        self.label_39.setText(
            QCoreApplication.translate("MainWindow", "Score filters:", None)
        )
        self.max_oh_lb_6.setText(
            QCoreApplication.translate("MainWindow", "Total score \u2265", None)
        )
        self.label_41.setText(
            QCoreApplication.translate("MainWindow", "Rank score \u2265", None)
        )
        self.Spectra_tab.setTabText(
            self.Spectra_tab.indexOf(self.tab_3),
            QCoreApplication.translate("MainWindow", "Lipostar2 (.csv)", None),
        )
        self.label_13.setText(
            QCoreApplication.translate(
                "MainWindow", "Predicted epilipidome (.json)", None
            )
        )
        self.input_ident_epilipidome_pb.setText(
            QCoreApplication.translate("MainWindow", "Open", None)
        )
        self.label_21.setText(
            QCoreApplication.translate("MainWindow", "Report folder", None)
        )
        self.report_folder_pb.setText(
            QCoreApplication.translate("MainWindow", "Choose path", None)
        )
        self.label_24.setText(
            QCoreApplication.translate("MainWindow", "Figure Quality", None)
        )
        self.dpi_cmb.setItemText(
            0, QCoreApplication.translate("MainWindow", "High (300dpi)", None)
        )
        self.dpi_cmb.setItemText(
            1, QCoreApplication.translate("MainWindow", "Medium (150 dpi)", None)
        )
        self.dpi_cmb.setItemText(
            2, QCoreApplication.translate("MainWindow", "Low (96dpi)", None)
        )

        self.label_22.setText(
            QCoreApplication.translate("MainWindow", "Export Table", None)
        )
        self.export_table_pb.setText(
            QCoreApplication.translate("MainWindow", "Choose path", None)
        )
        # if QT_CONFIG(tooltip)
        self.dump_json_chkb.setToolTip(
            QCoreApplication.translate(
                "MainWindow",
                "<html><head/><body><p>Save all data of identified lipids including observed peak lists and scoring details in a .json file.</p></body></html>",
                None,
            )
        )
        # endif // QT_CONFIG(tooltip)
        self.dump_json_chkb.setText(
            QCoreApplication.translate("MainWindow", "Dump all data (.json)", None)
        )
        self.run_identification_pb.setText(
            QCoreApplication.translate("MainWindow", "Hunt lipids!", None)
        )
        self.tabWidget.setTabText(
            self.tabWidget.indexOf(self.tab),
            QCoreApplication.translate("MainWindow", "Lipid identification", None),
        )
        self.label_12.setText(
            QCoreApplication.translate("MainWindow", "LipidLynxX cli_lynx path:", None)
        )
        self.pushButton_5.setText(
            QCoreApplication.translate("MainWindow", "Open", None)
        )
        self.pushButton_8.setText(
            QCoreApplication.translate("MainWindow", "Save", None)
        )
        self.tabWidget.setTabText(
            self.tabWidget.indexOf(self.tab_c),
            QCoreApplication.translate("MainWindow", "Settings", None),
        )
        self.label_33.setText(
            QCoreApplication.translate(
                "MainWindow",
                '<html><head/><body><p><span style=" font-size:18pt; font-weight:700;">LPPtiger 2 (Beta)</span></p></body></html>',
                None,
            )
        )
        self.label_29.setText(
            QCoreApplication.translate(
                "MainWindow",
                '<html><head/><body><p><span style=" font-size:10pt;">Copyright (C) 2021-2022 LMAI_team @ TU Dresden</span></p><p><span style=" font-size:10pt;">+ LMAI_team: Zhixu Ni, Maria Fedorova</span></p><p><span style=" font-size:10pt;">Copyright (C) 2016-2021 SysMedOs_team @ University of Leipzig</span></p><p><span style=" font-size:10pt;">+ SysMedOs_team: Zhixu Ni, Georgia Angelidou, Maria Fedorova</span></p></body></html>',
                None,
            )
        )
        self.label_31.setText(
            QCoreApplication.translate(
                "MainWindow",
                '<html><head/><body><p><a href="https://tu-dresden.de/med/mf/zml/forschungsgruppen/fedorova"><span style=" font-size:10pt; text-decoration: underline; color:#094fd1;">https://tu-dresden.de/med/mf/zml/forschungsgruppen/fedorova</span></a></p><p><span style=" font-size:10pt;">For more info please contact:</span></p><p><span style=" font-size:10pt;">+ Developer Zhixu Ni: zhixu.ni@tu-dresden.de</span></p><p><span style=" font-size:10pt;">+ LPPtiger2 repository: </span><a href="https://github.com/LMAI-TUD/lpptiger2"><span style=" font-size:10pt; text-decoration: underline; color:#094fd1;">https://github.com/LMAI-TUD/lpptiger2</span></a></p></body></html>',
                None,
            )
        )
        self.label_26.setText(
            QCoreApplication.translate(
                "MainWindow",
                '<html><head/><body><p><span style=" font-size:10pt;">LPPtiger2 is Dual-licensed:</span></p><p><span style=" font-size:10pt;">+ For academic and non-commercial use: GNU Affero General Public License (AGPL v3 License) </span></p><p><span style=" font-size:10pt;">       + https://www.gnu.org/licenses/agpl-3.0.en.html </span></p><p><span style=" font-size:10pt;">+ For commercial use: please contact the LMAI_team by email.</span></p></body></html>',
                None,
            )
        )
        self.label_30.setText(
            QCoreApplication.translate(
                "MainWindow",
                '<html><head/><body><p><span style=" font-size:10pt;">Please cite our publication in an appropriate form:</span></p><p><span style=" font-size:10pt;">Ni, Zhixu, Georgia Angelidou, Ralf Hoffmann, and Maria Fedorova.</span></p><p><span style=" font-size:10pt;">LPPtiger software for lipidome-specific prediction and identification of oxidized phospholipids from LC-MS datasets</span></p><p><span style=" font-size:10pt;">Scientific Reports 7, Article number: 15138 (2017). DOI: 10.1038/s41598-017-15363-z</span></p></body></html>',
                None,
            )
        )
        self.label_35.setText(
            QCoreApplication.translate(
                "MainWindow",
                '<html><head/><body><p><span style=" font-size:10pt;">*. LPPtiger 2 is powered by opensource projects: Matplotlib, Natsort, Numpy, Pandas, pymzml, PySide6, rdkit</span></p></body></html>',
                None,
            )
        )
        self.tabWidget.setTabText(
            self.tabWidget.indexOf(self.tab_z),
            QCoreApplication.translate("MainWindow", "About", None),
        )
        self.label_34.setText(
            QCoreApplication.translate(
                "MainWindow",
                '<html><head/><body><p><span style=" font-size:10pt;">LPPtiger 2 | Copyright (C) 2021-2022 LMAI_team @ TU Dresden  # https://github.com/LMAI-TUD/lpptiger2</span></p></body></html>',
                None,
            )
        )

    # retranslateUi
