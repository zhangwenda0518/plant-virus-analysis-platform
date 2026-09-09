# uncompyle6 version 3.9.3
# Python bytecode version base 3.7.0 (3394)
# Decompiled from: Python 3.12.10 (tags/v3.12.10:0cc8128, Apr  8 2025, 12:21:36) [MSC v.1943 64 bit (AMD64)]
# Embedded file name: mutation_tools/site_counter.py
"""
Author: Zhou Zhi-Jian

"""
import os, traceback, config, datetime, pandas as pd, numpy as np, platform
from PySide2.QtWidgets import QApplication, QWidget, QRadioButton, QVBoxLayout, QHBoxLayout, QFileDialog, QMessageBox
from PySide2.QtGui import QIcon, QPalette
from PySide2.QtCore import Qt, QThread, Signal
from custom_libraries.gui_class import QLineEdit_get, PopTipWindow, ErrorLogGUI, QProgressBar_beautify, Start_Button, InputFile_Button, MyGroupBox
from custom_libraries.public_functions import resolve_file_path, make_dir

class SiteCounts(QWidget):
    __doc__ = "\n    位点计数器GUI\n    "

    def __init__(self, parent=None):
        super(SiteCounts, self).__init__(parent)
        self.setAcceptDrops(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, Qt.white)
        self.setPalette(palette)
        self.setWindowIcon(QIcon(":/site_count.png"))
        self.setWindowTitle("Site Counter")
        desktop = QApplication.desktop()
        des_w = desktop.width()
        des_h = desktop.height()
        win_w = des_w * 0.3
        win_h = des_h * 0.35
        self.resize(int(win_w), int(win_h))
        self.initUI()

    def initUI(self):
        self.setToolTip("提示：<b>此功能用于计算对齐的序列集中每个位点的核苷酸/氨基酸种类、数量和最丰度核苷酸/氨基酸占比，并生成同一性序列。</b>")
        self.file_edit = QLineEdit_get()
        self.file_button = InputFile_Button()
        self.file_button.setToolTip("Or click to load a <b> alignment sequence-file with fasta format </b>")
        self.file_button.clicked.connect(self.openfile)
        box1 = MyGroupBox("&Input file", self)
        box1.setMaximumHeight(100)
        box1.setToolTip("Drag a <b> alignment sequence-file with fasta format </b> and drop here")
        hbox1 = QHBoxLayout()
        hbox1.addWidget(self.file_edit)
        hbox1.addWidget(self.file_button)
        box1.setLayout(hbox1)
        self.data_type_box = MyGroupBox("&Data type", self)
        self.data_type_box.setMaximumHeight(100)
        self.data_hbox = QHBoxLayout()
        self.Radio1 = QRadioButton("Nucleotide")
        self.Radio1.setChecked(True)
        self.Radio2 = QRadioButton("Protein")
        self.data_hbox.addWidget(self.Radio1)
        self.data_hbox.addWidget(self.Radio2)
        self.data_type_box.setLayout(self.data_hbox)
        self.progressBar = QProgressBar_beautify()
        self.progressBar.setValue(0)
        self.start_button = Start_Button()
        self.start_button.clicked.connect(self.works)
        progressBar_hbox = QHBoxLayout()
        progressBar_hbox.addWidget(self.progressBar)
        progressBar_hbox.addWidget(self.start_button)
        progressBar_box = MyGroupBox("&Run and Progress", self)
        progressBar_box.setMaximumHeight(100)
        progressBar_box.setLayout(progressBar_hbox)
        wins_V_layout = QVBoxLayout()
        wins_V_layout.setSpacing(20)
        wins_V_layout.addWidget(box1)
        wins_V_layout.addWidget(self.data_type_box)
        wins_V_layout.addWidget(progressBar_box)
        self.setLayout(wins_V_layout)

    def openfile(self):
        """
        按钮文件加载文件路径
        :return: 
        """
        self.fname = QFileDialog.getOpenFileName(self, "Select an aligned sequence file", "/home")
        self.file_path = self.fname[0]
        if self.file_path != "":
            self.file_edit.setText(self.file_path)

    def works(self):
        """
        运行前检查参数和配置
        :return: 
        """
        now = datetime.datetime.now()
        time_now = now.strftime("%Y-%m-%d %H:%M:%S")
        run_id = time_now.replace("-", "").replace(" ", "_").replace(":", "")
        fasta_file_path = self.file_edit.text().replace("\\", "/")
        if fasta_file_path == "" or os.path.isdir(fasta_file_path):
            return
        else:
            input_data_dir, out_prefix = resolve_file_path(fasta_file_path)
            save_file_dir = input_data_dir + "/site_count_" + run_id + "/" + out_prefix
            save_dir = os.path.dirname(fasta_file_path) + "/site_count_" + run_id
            make_dir(save_dir)
            select = ""
            if self.Radio1.isChecked():
                select = "Radio1"
            else:
                if self.Radio2.isChecked():
                    select = "Radio2"
        self.start_button.setEnabled(False)
        self.start_button.setText("Run...")
        self.thread_e = Thread_site(fasta_file_path, save_file_dir, save_dir, select)
        self.thread_e.site_trigger.connect(self.set_progressbar_value)
        self.thread_e.trigger_save.connect(self.run_finshed)
        self.thread_e.start()

    def set_progressbar_value(self, value):
        """
        进度条属性
        :param value:
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

    def run_finshed(self, savedir):
        PopTipWindow.success("Successfully produced!", "<b>Saved in: </b>" + savedir)

    def closeEvent(self, QCloseEvent):
        """
        重写窗口关闭事件
        :param QCloseEvent:
        :return:
        """
        if self.start_button.text() == "Run...":
            a = QMessageBox.question(self, "Exit?", "The task is running. Are you sure to exit?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if a == QMessageBox.Yes:
                pass
            else:
                QCloseEvent.ignore()
                return
        self.progressBar.setValue(0)
        self.start_button.setEnabled(True)
        self.start_button.setText("Start")
        if platform.system().lower() == "windows":
            try:
                self.thread_e.terminate()
            except:
                pass

            try:
                self.error.close()
            except:
                pass

        self.close()


class Thread_site(QThread):
    __doc__ = "\n    位点计数运行线程\n    "
    site_trigger = Signal(int)
    trigger_save = Signal(str)

    def __init__(self, fasta_file_path, save_file_dir, save_dir, select):
        super(Thread_site, self).__init__()
        self.fasta_file_path = fasta_file_path
        self.save_file_dir = save_file_dir
        self.save_dir = save_dir
        self.select = select

    def run(self):
        try:
            try:
                self.site_trigger.emit(10)
                sequence = []
                seq_dic = {}
                seq_name_list = []
                seq_name = ""
                with open((self.fasta_file_path), "r", encoding="utf-8") as gen_file:
                    for line in gen_file:
                        line = line.strip()
                        if line.startswith(">"):
                            seq_name = line.strip(">")
                            seq_name_list.append(seq_name)
                            seq_dic[seq_name] = []
                        else:
                            seq_dic[seq_name].append(line)

                    for each_seqname in seq_name_list:
                        seqs = "".join(seq_dic[each_seqname])
                        sequence.append(seqs.upper())

                self.site_trigger.emit(50)
                type_kinds = ['A', 'G', 'C', 'T', '-']
                if self.select == "Radio1" and sequence[0].find("U") != -1:
                    type_kinds = [
                     'A', 'G', 'C', 'U', 
                     '-']
                else:
                    if self.select == "Radio2":
                        type_kinds = [
                         'A', 'C', 
                         'D', 'E', 'F', 'G', 'H', 
                         'I', 'K', 'L', 'M', 'N', 
                         'P', 'Q', 'R', 'S', 'T', 
                         'V', 'W', 'Y', '-']
                    site_count = len(sequence[0].strip())
                    d_site = []
                    for site in range(1, site_count + 1):
                        d_site.append("site" + str(site))

                    arr1 = np.zeros((len(type_kinds), site_count))
                    for site in range(1, site_count + 1):
                        for line in sequence:
                            s = line[(site - 1)[:site]]
                            for i in range(len(type_kinds)):
                                if s == type_kinds[i]:
                                    arr1[i][site - 1] += 1

                    data1 = pd.DataFrame(arr1, columns=d_site)
                    data1.index = type_kinds
                    data1.T.to_csv(self.save_file_dir + "_site_count.csv")
                    arr2 = np.zeros((len(type_kinds), site_count))
                    for row in range(len(type_kinds)):
                        for column in range(site_count):
                            arr2[row][column] = arr1[row][column] / len(sequence)

                    data2 = pd.DataFrame(arr2, columns=d_site)
                    data2.index = type_kinds
                    data2.T.to_csv(self.save_file_dir + "_site_proportion.csv")
                    self.site_trigger.emit(70)
                    con = []
                    max_probability = np.zeros((site_count, 1))
                    for column in range(site_count):
                        c = []
                        for row in range(len(type_kinds)):
                            c.append(arr2[row][column])

                        max_number = int(c.index(max(c)))
                        con.append(type_kinds[max_number])
                        max_probability[column][0] = max(c)

                    self.site_trigger.emit(80)
                    out_consensus = open((self.save_file_dir + "_consensus.fasta"), "w",
                      encoding="utf-8")
                    out_consensus.write(">consensus seq\n" + "".join(con))
                    data3 = pd.DataFrame(max_probability, columns=["proportion"])
                    col_name = data3.columns.tolist()
                    col_name.insert(0, "type_kinds")
                    data3 = data3.reindex(columns=col_name)
                    data3["type_kinds"] = con
                    e_site = []
                    for si in d_site:
                        si = si.replace("site", "")
                        e_site.append(si)

                    data3.index = e_site
                    data3.to_csv(self.save_file_dir + "_site_max_proportion.csv")
                    out_consensus.close()
                    self.site_trigger.emit(100)
                    self.trigger_save.emit(self.save_dir)
            except Exception as e:
                try:
                    config.run_error_log = "traceback.format_exc():\n%s" % traceback.format_exc()
                    self.site_trigger.emit(1314)
                finally:
                    e = None
                    del e

        finally:
            self.exit(0)
