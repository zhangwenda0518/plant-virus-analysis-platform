# uncompyle6 version 3.9.3
# Python bytecode version base 3.7.0 (3394)
# Decompiled from: Python 3.12.10 (tags/v3.12.10:0cc8128, Apr  8 2025, 12:21:36) [MSC v.1943 64 bit (AMD64)]
# Embedded file name: mutation_tools/site_screen.py
"""
Author: Zhou Zhi-Jian

"""
import os, traceback, config, platform
from PySide2.QtWidgets import QApplication, QWidget, QPushButton, QVBoxLayout, QHBoxLayout, QFileDialog, QHeaderView, QMessageBox
from PySide2.QtGui import QIcon, QPalette, QFont
from PySide2.QtCore import Qt, QThread, Signal
from custom_libraries.gui_class import QLineEdit_get, PopTipWindow, ErrorLogGUI, QProgressBar_beautify, TableWithCopy, Start_Button, InputFile_Button, MyGroupBox
from custom_libraries.public_functions import fast_read_fasta, resolve_file_path

class SiteScreening(QWidget):
    __doc__ = "\n    位点扫描运行线程\n    "

    def __init__(self):
        super(SiteScreening, self).__init__()
        self.setAcceptDrops(True)
        desktop = QApplication.desktop()
        des_w = desktop.width()
        des_h = desktop.height()
        palette = QPalette()
        palette.setColor(QPalette.Background, Qt.white)
        self.setPalette(palette)
        self.setWindowIcon(QIcon(":/linked.png"))
        self.setWindowTitle("Site Scree")
        win_w = des_w * 0.35
        win_h = des_h * 0.5
        self.resize(int(win_w), int(win_h))
        self.initUI()

    def initUI(self):
        self.setToolTip("提示：<b>此功能用于提取在多个位点同时包含指定的碱基或氨基酸的序列</b>，可用于探索不同位点之间的关联突变。")
        self.file_edit = QLineEdit_get()
        self.file_button = InputFile_Button()
        self.file_button.setToolTip("Or click to load a <b>sequence file with fasta format</b> ")
        self.file_button.clicked.connect(self.openfile)
        self.file_button.setMinimumSize(30, 25)
        self.file_button.setMaximumSize(35, 30)
        self.fasfile_hbox1 = QHBoxLayout()
        self.fasfile_hbox1.addWidget(self.file_edit)
        self.fasfile_hbox1.addWidget(self.file_button)
        self.fasfile_box1 = MyGroupBox("&Input file", self)
        self.fasfile_box1.setToolTip("Drag a <b>sequence file with fasta format</b>  and drop here")
        self.fasfile_box1.setLayout(self.fasfile_hbox1)
        self.table = TableWithCopy()
        self.table.setRowCount(4)
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(["Position", "Label (nt or aa)"])
        self.table.setFont(QFont(config.app_font, 12))
        self.table.setSortingEnabled(True)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.verticalHeader().setVisible(True)
        self.add_row_button = QPushButton("Add row")
        self.add_row_button.clicked.connect(self.add_new_row)
        self.delte_row_button = QPushButton("Delete row")
        self.delte_row_button.clicked.connect(self.delte_row)
        table_hbox = QHBoxLayout()
        table_hbox.addWidget(self.table)
        row_change = QHBoxLayout()
        row_change.addWidget(self.add_row_button)
        row_change.addWidget(self.delte_row_button)
        table_vbox = QVBoxLayout()
        table_vbox.addLayout(table_hbox)
        table_vbox.addLayout(row_change)
        self.table_box = MyGroupBox("&Information of site(s)", self)
        self.table_box.setLayout(table_vbox)
        self.table_box.setToolTip("Please enter or pasta <b>position and nt(or aa) label</b>")
        self.progressBar = QProgressBar_beautify()
        self.progressBar.setValue(0)
        self.button = Start_Button()
        self.button.clicked.connect(self.Work)
        self.progressBar_hbox = QHBoxLayout()
        self.progressBar_hbox.addWidget(self.progressBar)
        self.progressBar_hbox.addWidget(self.button)
        self.progressBar_box = MyGroupBox("&Run and Progress", self)
        self.progressBar_box.setLayout(self.progressBar_hbox)
        self.win_layout = QVBoxLayout()
        self.win_layout.setSpacing(20)
        self.win_layout.addWidget(self.fasfile_box1)
        self.win_layout.addWidget(self.table_box)
        self.win_layout.addWidget(self.progressBar_box)
        self.setLayout(self.win_layout)

    def add_new_row(self):
        cur_row_count = self.table.rowCount()
        cur_row_count += 1
        self.table.setRowCount(cur_row_count)

    def delte_row(self):
        cur_row_count = self.table.rowCount()
        cur_row_count -= 1
        self.table.setRowCount(cur_row_count)

    def openfile(self):
        """
        按钮加载文件路径
        :return:
        """
        self.fasname = QFileDialog.getOpenFileName(self, "Select an aligned sequence file", "/home")
        self.fasfile_path = self.fasname[0]
        if self.fasfile_path != "":
            self.file_edit.setText(self.fasfile_path)

    def closeEvent(self, QCloseEvent):
        """
        重写窗口关闭事件
        :param QCloseEvent:
        :return:
        """
        if self.button.text() == "Run...":
            a = QMessageBox.question(self, "Exit?", "The task is running. Are you sure to exit?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if a == QMessageBox.Yes:
                pass
            else:
                QCloseEvent.ignore()
                return
        self.progressBar.setValue(0)
        self.button.setEnabled(True)
        self.button.setText("Start")
        if platform.system().lower() == "windows":
            try:
                self.scree_gene_thread.terminate()
            except:
                pass

            try:
                self.error.close()
            except:
                pass

        self.close()

    def Work(self):
        """
        运行前检查参数和配置
        :return:
        """
        self.progressBar.setValue(0)
        fasfile_path = self.file_edit.text()
        if fasfile_path == "" or os.path.isdir(fasfile_path):
            return
        else:
            row_count = self.table.rowCount()
            colum_count = self.table.columnCount()
            qury_str_table = []
            for i in range(row_count):
                if self.table.item(i, 0) and self.table.item(i, 1):
                    new_list = []
                    for j in range(colum_count):
                        new_list.append(self.table.item(i, j).text())

                    qury_str_table.append("\t".join(new_list))

            return qury_str_table or None
        self.button.setEnabled(False)
        self.button.setText("Run...")
        self.scree_gene_thread = RunThread_i(fasfile_path, qury_str_table)
        self.scree_gene_thread.single_progress.connect(self.set_progressbar_value)
        self.scree_gene_thread.single_save.connect(self.run_finshed)
        self.scree_gene_thread.start()

    def set_progressbar_value(self, value):
        """
        设置进度条
        :param value:
        :return:
        """
        self.progressBar.setValue(value)
        if value == 100:
            self.button.setEnabled(True)
            self.button.setText("Start")
            return
        if value == 1314:
            self.error = ErrorLogGUI()
            self.error.show()
            self.button.setEnabled(True)
            self.button.setText("Start")
            return

    def run_finshed(self, save_filepath):
        PopTipWindow.success("Successfully produced!", "<b>Saved in: </b>" + save_filepath)


class RunThread_i(QThread):
    __doc__ = "\n    位点扫描运行线程\n    "
    single_progress = Signal(int)
    single_save = Signal(str)

    def __init__(self, fasfile_path, qury_str_table):
        super(RunThread_i, self).__init__()
        self.fasfile_path = fasfile_path.replace("\\", "/")
        self.qury_str_table = qury_str_table

    def run(self):
        try:
            try:
                self.single_progress.emit(10)
                input_data_dir, out_prefix = resolve_file_path(self.fasfile_path)
                fassave_dir_path = input_data_dir + "/" + out_prefix
                seq_lis = fast_read_fasta(self.fasfile_path)
                seq_count = len(seq_lis)
                offset = int(seq_count / 2)
                self.single_progress.emit(20)
                len_tab = len(self.qury_str_table)
                out_name = []
                for infor in self.qury_str_table:
                    posi, nt_aa = infor.split("\t")
                    out_name.append(posi + "-" + nt_aa)

                save_file_path = fassave_dir_path + "_".join(out_name) + ".fasta"
                out = open(save_file_path, "w", encoding="utf-8")
                for n in range(len(seq_lis)):
                    seq_name = seq_lis[n][0]
                    seq = seq_lis[n][1]
                    k = 0
                    for line in self.qury_str_table:
                        site, zifu = line.split("\t")
                        site = int(site)
                        if seq[(site - 1)[:site]].upper() == zifu.upper():
                            k += 1

                    if k == len_tab:
                        out.write(">" + seq_name + "\n" + seq + "\n")
                    offset = offset + 1
                    proess = offset / seq_count * 50 - 1
                    self.single_progress.emit(int(proess))

                self.single_progress.emit(100)
                self.single_save.emit(save_file_path)
            except Exception as e:
                try:
                    config.run_error_log = "traceback.format_exc():\n%s" % traceback.format_exc()
                    self.single_progress.emit(1314)
                finally:
                    e = None
                    del e

        finally:
            self.exit(0)
