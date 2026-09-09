# uncompyle6 version 3.9.3
# Python bytecode version base 3.7.0 (3394)
# Decompiled from: Python 3.12.10 (tags/v3.12.10:0cc8128, Apr  8 2025, 12:21:36) [MSC v.1943 64 bit (AMD64)]
# Embedded file name: phylogeny_tools/partition_make.py
import traceback
from PySide2.QtWidgets import QApplication, QWidget, QPushButton, QHBoxLayout, QVBoxLayout, QLabel, QComboBox, QFileDialog, QHeaderView, QRadioButton
from PySide2.QtGui import QIcon, QPalette, QFont
from PySide2.QtCore import Qt, QThread, Signal
import config
from custom_libraries.gui_class import QLineEdit_getDir, PopTipWindow, TableWithCopy, ErrorLogGUI, QProgressBar_beautify, Start_Button, InputFile_Button, MyGroupBox

class PartitionMake(QWidget):

    def __init__(self, parent=None):
        super(PartitionMake, self).__init__(parent)
        self.setAcceptDrops(True)
        desktop = QApplication.desktop()
        des_w = desktop.width()
        des_h = desktop.height()
        palette = QPalette()
        palette.setColor(QPalette.Background, Qt.white)
        self.setPalette(palette)
        self.setWindowIcon(QIcon(":/partition.png"))
        self.setWindowTitle("Partition-Model Make")
        win_w = des_w * 0.35
        win_h = des_h * 0.5
        self.resize(int(win_w), int(win_h))
        self.initUI()

    def initUI(self):
        self.setToolTip("提示：<b>此功能用于制作序列的分区信息文件</b>，可用于Mrbayes和IQ-tree的多基因联合建树。")
        self.table = TableWithCopy()
        self.table.setRowCount(4)
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(["Partition name", "Start", "End", "Model"])
        self.table.setFont(QFont(config.app_font, config.app_font_size))
        self.table.setSortingEnabled(True)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.verticalHeader().setVisible(True)
        self.add_row_button = QPushButton("Add row")
        self.add_row_button.clicked.connect(self.add_new_row)
        self.delte_row_button = QPushButton("Delete row")
        self.delte_row_button.clicked.connect(self.delte_row)
        tablebox = MyGroupBox("Input information of partition and model", self)
        tablebox.setToolTip("The example format of <b>model column</b>: GTR+G or HKY+G+F or GTR+G+I+F")
        file_hbox = QHBoxLayout()
        file_hbox.addWidget(self.table)
        row_change = QHBoxLayout()
        row_change.addWidget(self.add_row_button)
        row_change.addWidget(self.delte_row_button)
        table_vbox = QVBoxLayout()
        table_vbox.addLayout(file_hbox)
        table_vbox.addLayout(row_change)
        tablebox.setLayout(table_vbox)
        software_label = QLabel("Partition file for:")
        self.model_for = QComboBox(self)
        self.model_for.addItem("IQ-Tree")
        self.model_for.addItem("MrBayes")
        data_type = QLabel("Data type:")
        self.nt_select = QRadioButton("Nucleotide")
        self.aa_select = QRadioButton("Amino acid")
        para_box = MyGroupBox("Parameter", self)
        para_h = QHBoxLayout()
        para_h.addWidget(software_label)
        para_h.addWidget(self.model_for)
        data_h = QHBoxLayout()
        data_h.addWidget(data_type)
        data_h.addWidget(self.nt_select)
        data_h.addWidget(self.aa_select)
        outdir_label = QLabel("Output directory:")
        self.outdir_edit = QLineEdit_getDir()
        self.outdir_edit.setToolTip("Drag a <b>folder</b> (for export result) and drop here")
        self.opendir_button = InputFile_Button()
        self.opendir_button.setToolTip("Or click to select a <b>folder</b>")
        self.opendir_button.clicked.connect(self.opendir)
        outfile_hbox = QHBoxLayout()
        outfile_hbox.addWidget(outdir_label)
        outfile_hbox.addWidget(self.outdir_edit)
        outfile_hbox.addWidget(self.opendir_button)
        para_V = QVBoxLayout()
        para_V.addLayout(para_h)
        para_V.addLayout(data_h)
        para_V.addLayout(outfile_hbox)
        para_box.setLayout(para_V)
        self.progressBar = QProgressBar_beautify()
        self.progressBar.setValue(0)
        self.start_button = Start_Button()
        self.start_button.clicked.connect(self.Work)
        progressBar_hbox = QHBoxLayout()
        progressBar_hbox.addWidget(self.progressBar)
        progressBar_hbox.addWidget(self.start_button)
        progressBar_box = MyGroupBox("&Run and Progress", self)
        progressBar_box.setLayout(progressBar_hbox)
        win_layout_V = QVBoxLayout()
        win_layout_V.setSpacing(20)
        win_layout_V.addWidget(tablebox)
        win_layout_V.addWidget(para_box)
        win_layout_V.addWidget(progressBar_box)
        self.setLayout(win_layout_V)

    def add_new_row(self):
        cur_row_count = self.table.rowCount()
        cur_row_count += 1
        self.table.setRowCount(cur_row_count)

    def delte_row(self):
        cur_row_count = self.table.rowCount()
        cur_row_count -= 1
        self.table.setRowCount(cur_row_count)

    def opendir(self):
        """
        按钮加载文件夹路径
        :return:
        """
        self.outdir_edit.setText("")
        loaddir = QFileDialog.getExistingDirectory(self, "Select a directory", "/home")
        self.outdir_edit.setText(loaddir.rstrip("/"))

    def run_finshed(self, fassave_dir):
        """
        运行结束时弹窗
        :param fassave_filepath: 储存文件路径
        :return:
        """
        PopTipWindow.success("Successfully produced!", "<b>Saved in: </b>" + fassave_dir)

    def set_progressbar_value(self, value):
        """
        进度条,如果收到“1314”信号，则弹出error log窗口
        :param value: 进度条值
        :return:
        """
        self.progressBar.setValue(value)
        if value == 100:
            self.start_button.setEnabled(True)
            self.start_button.setText("Start")
            return
        if value == 1314:
            self.error = ErrorLogGUI()
            self.error.show()
            self.start_button.setEnabled(True)
            self.start_button.setText("Start")
            return

    def closeEvent(self, event):
        self.progressBar.setValue(0)
        self.start_button.setEnabled(True)
        self.start_button.setText("Start")
        self.close()

    def Work(self):
        out_dir_path = self.outdir_edit.text().replace("\\", "/")
        if out_dir_path == "":
            return
        row_count = self.table.rowCount()
        colum_count = self.table.columnCount()
        part_list = []
        for i in range(row_count):
            if self.table.item(i, 0) != None:
                new_list = []
                for j in range(0, colum_count):
                    if self.table.item(i, j) == None:
                        PopTipWindow.error("Error!", "Missing data in row " + str(i + 1) + " and column " + str(j + 1))
                        return
                    new_list.append(self.table.item(i, j).text())

                part_list.append(new_list)

        if part_list == []:
            return
        used_sorftware = self.model_for.currentText()
        if used_sorftware == "MrBayes":
            support_model = [
             'JC', 'F81', 'K80', 'K2P', 'HKY', 'SYM', 
             'GTR', 
             'JTT', 'WAG', 'LG', 
             'VT', 'BLOSUM62', 'BLOSUM', 'EQUALIN', 
             'MIXED', 
             'POISSON', 'CPREV', 'RTREV', 'DAYHOFF', 'MTMAM', 
             'MTREV']
            for part in part_list:
                key = part[0]
                model_str = part[-1].upper()
                if model_str.count("+") == 0:
                    if model_str not in support_model:
                        PopTipWindow.error("Error!", model_str + " model of partition " + key + " is not supported in mrbayes!")
                        return
                else:
                    model_str_list = model_str.split("+")
                    for each_model in model_str_list:
                        if len(each_model) >= 2 and each_model not in support_model:
                            PopTipWindow.error("Error!", each_model + " model of partition " + key + " is not supported in mrbayes!")
                            return

        data_type = "nt"
        if self.aa_select.isChecked():
            data_type = "aa"
        self.thread_partmake = RunThread_partmake(part_list, used_sorftware, data_type, out_dir_path)
        self.thread_partmake.single_progress.connect(self.set_progressbar_value)
        self.thread_partmake.single_save.connect(self.run_finshed)
        self.start_button.setEnabled(False)
        self.start_button.setText("Run...")
        self.thread_partmake.start()


class RunThread_partmake(QThread):
    single_progress = Signal(int)
    single_save = Signal(str)

    def __init__(self, part_list, used_sorftware, data_type, out_dir_path):
        super(RunThread_partmake, self).__init__()
        self.part_list = part_list
        self.used_sorftware = used_sorftware
        self.data_type = data_type
        self.out_dir_path = out_dir_path

    def run(self):
        try:
            try:
                contain_list = []
                if self.used_sorftware == "IQ-Tree":
                    contain_list.append("#nexus")
                    contain_list.append("begin sets;")
                    model_infor = []
                    for part in self.part_list:
                        key = part[0]
                        contain_list.append("charset " + key + " = " + str(part[1]) + "-" + str(part[2]) + ";")
                        model_infor.append(part[-1].upper() + ":" + key)

                    contain_list.append("charpartition partitions = " + ", ".join(model_infor) + ";")
                    contain_list.append("end;")
                else:
                    if self.used_sorftware == "MrBayes":
                        if self.data_type == "nt":
                            lset_list = []
                            num_id = 1
                            part_name_list = []
                            for part in self.part_list:
                                key = part[0]
                                part_name_list.append(key)
                                contain_list.append("charset " + key + " = " + str(part[1]) + "-" + str(part[2]) + ";")
                                model_str = part[-1].upper()
                                model_set, rates_set = self.model_resolution(model_str, "nt")
                                lset_list.append("lset applyto=(" + str(num_id) + ") " + model_set + rates_set + ";")
                                num_id += 1

                            contain_list.append("partition Names = " + str(len(self.part_list)) + ":" + ", ".join(part_name_list) + ";")
                            contain_list.append("set partition=Names;")
                            for line in lset_list:
                                contain_list.append(line)

                            contain_list.append("prset applyto=(all) ratepr=variable;")
                            contain_list.append("unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all) tratio=(all);")
                        else:
                            if self.data_type == "aa":
                                num_ids = 1
                                model_para_list = []
                                part_name_list = []
                                for part in self.part_list:
                                    key = part[0]
                                    part_name_list.append(key)
                                    contain_list.append("charset " + key + " = " + str(part[1]) + "-" + str(part[2]) + ";")
                                    model_str = part[-1].upper()
                                    model_set, rates_set = self.model_resolution(model_str, "aa")
                                    model_para_list.append("lset applyto=(" + str(num_ids) + ")" + rates_set + ";")
                                    model_para_list.append("prset applyto=(" + str(num_ids) + ")" + " " + model_set + ";")
                                    if model_str.find("+F") != -1:
                                        model_para_list.append("prset applyto=(" + str(num_ids) + ")" + " statefreqpr=fixed(empirical);")
                                    num_ids += 1

                                contain_list.append("partition Names = " + str(len(self.part_list)) + ":" + ", ".join(part_name_list) + ";")
                                contain_list.append("set partition=Names;")
                                for line in model_para_list:
                                    contain_list.append(line)

                                contain_list.append("prset applyto=(all) ratepr=variable;")
                                contain_list.append("unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all) tratio=(all);")
                    with open((self.out_dir_path + "/partition_file.txt"), "w", encoding="utf-8") as outfile:
                        outfile.write("\n".join(contain_list))
                    self.single_progress.emit(100)
                    self.single_save.emit(self.out_dir_path + "/partition_file.txt")
            except Exception as e:
                try:
                    config.run_error_log = "traceback.format_exc():\n%s" % traceback.format_exc()
                    self.single_progress.emit(1314)
                finally:
                    e = None
                    del e

        finally:
            self.exit(0)

    def model_resolution(self, model_str, data_type):
        nt_model_dic = {
         'JC': '"nst=1"', 
         'F81': '"nst=1"', 
         'K80': '"nst=2"', 
         'K2P': '"nst=2"', 
         'HKY': '"nst=2"', 
         'SYM': '"nst=6"', 
         'GTR': '"nst=6"'}
        aa_model_dic = {
         'JTT': '"Aamodelpr=fixed(jones)"', 
         'WAG': '"Aamodelpr=fixed(wag)"', 
         'LG': '"Aamodelpr=fixed(lg)"', 
         'VT': '"Aamodelpr=fixed(vt)"', 
         'BLOSUM62': '"Aamodelpr=fixed(blosum62)"', 
         'BLOSUM': '"Aamodelpr=fixed(blosum)"', 
         'EQUALIN': '"Aamodelpr=fixed(equalin)"', 
         'GTR': '"Aamodelpr=fixed(gtr)"', 
         'MIXED': '"Aamodelpr=fixed(mixed)"', 
         'POISSON': '"Aamodelpr=fixed(poisson)"', 
         'CPREV': '"Aamodelpr=fixed(cprev)"', 
         'RTREV': '"Aamodelpr=fixed(rtrev)"', 
         'DAYHOFF': '"Aamodelpr=fixed(dayhoff)"', 
         'MTMAM': '"Aamodelpr=fixed(mtmam)"', 
         'MTREV': '"Aamodelpr=fixed(mtrev)"'}
        model_set = ""
        rates_set = ""
        if data_type.upper() == "NT":
            if model_str.count("+") == 0:
                model_set = nt_model_dic[model_str]
            else:
                model_str_list = model_str.split("+")
                for each_model in model_str_list:
                    if each_model in ('JC', 'F81', 'K80', 'K2P', 'HKY', 'SYM', 'GTR'):
                        model_set = nt_model_dic[each_model]

            if model_str.find("+I") != -1 and model_str.find("+G") != -1:
                rates_set = " rates=invgamma"
            else:
                if model_str.find("+G") != -1:
                    rates_set = " rates=gamma"
                else:
                    if model_str.find("+I") != -1:
                        rates_set = " rates=propinv"
            return (
             model_set, rates_set)
        if data_type.upper() == "AA":
            if model_str.find("+I") != -1 and model_str.find("+G") != -1:
                rates_set = " rates=invgamma"
            else:
                if model_str.find("+G") != -1:
                    rates_set = " rates=gamma"
                else:
                    if model_str.find("+I") != -1:
                        rates_set = " rates=propinv"
                    elif model_str.count("+") == 0:
                        model_set = aa_model_dic[model_str]
                    else:
                        model_str_list = model_str.split("+")
                        for each_model in model_str_list:
                            if each_model in ('JTT', 'WAG', 'LG', 'VT', 'BLOSUM62',
                                              'BLOSUM', 'EQUALIN', 'MIXED', 'POISSON',
                                              'CPREV', 'RTREV', 'DAYHOFF', 'MTMAM',
                                              'MTREV', 'GTR'):
                                model_set = aa_model_dic[each_model]

                    return (
                     model_set, rates_set)
