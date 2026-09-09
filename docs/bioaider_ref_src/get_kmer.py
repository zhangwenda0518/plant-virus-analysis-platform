# uncompyle6 version 3.9.3
# Python bytecode version base 3.7.0 (3394)
# Decompiled from: Python 3.12.10 (tags/v3.12.10:0cc8128, Apr  8 2025, 12:21:36) [MSC v.1943 64 bit (AMD64)]
# Embedded file name: clustering_tools/get_kmer.py
import os, traceback, platform, multiprocessing
from multiprocessing import Pool
from functools import reduce
from collections import Counter
from PySide2.QtWidgets import QApplication, QWidget, QLineEdit, QHBoxLayout, QFileDialog, QComboBox, QLabel, QVBoxLayout, QMessageBox
from PySide2.QtGui import QIcon, QPalette
from PySide2.QtCore import Qt, QThread, Signal
import config
from custom_libraries.gui_class import QLineEdit_get, PopTipWindow, ErrorLogGUI, QProgressBar_beautify, Start_Button, InputFile_Button, MyGroupBox
from custom_libraries.public_functions import resolve_file_path, fast_read_fasta

def piecewise_read_seq(input_file_path, max_num_seq, run_nums):
    """
    分段读取数据到内存进行计算
    :param input_file_path: 输入序列名称
    :param max_num_seq: 每次最大读取的序列条数
    :param run_nums: 第几次读取
    :return:
    """
    seq_tab = []
    seq_dic = {}
    seq_name_list = []
    seq_name = ""
    with open(input_file_path, "r", encoding="utf-8") as gen_file:
        for line in gen_file:
            line = line.strip()
            if line.startswith(">"):
                seq_name = line.strip(">")
                seq_name_list.append(seq_name)
                seq_dic[seq_name] = []
            elif line != "":
                seq_dic[seq_name].append(line)

        n_start = max_num_seq * (run_nums - 1)
        n_end = max_num_seq * run_nums
        m = len(seq_name_list) - 1
        for n in range(len(seq_name_list)):
            if n_start <= n < n_end:
                each_seqname = seq_name_list[n]
                seq = "".join(seq_dic[each_seqname]).upper()
                seq_tab.append([each_seqname, seq])

        if m < n_end:
            run_flage = False
        else:
            run_flage = True
    return (
     m + 1, seq_tab, run_flage)


def calculate_kmers(seq_name, seq, k_size, all_kmer_type):
    """

    :param seq: 传入一条序列
    :param k_size: 输入k-mer对应的K值
    :return: 返回一个计数的字典
    """
    kmers = []
    n_kmers = len(seq) - k_size + 1
    for i in range(n_kmers):
        kmer = seq[i[:i + k_size]]
        kmers.append(kmer)

    kmers_stat = dict(Counter(kmers))
    out_line = [
     seq_name]
    for each_kind in all_kmer_type:
        if each_kind in kmers_stat:
            out_line.append(str(kmers_stat[each_kind]))
        else:
            out_line.append("0")

    return "\t".join(out_line)


class GetKmer(QWidget):
    __doc__ = "\n    获取序列的k-mer矩阵\n    "

    def __init__(self, parent=None):
        super(GetKmer, self).__init__(parent)
        self.setAcceptDrops(True)
        palette = QPalette()
        palette.setColor(QPalette.Background, Qt.white)
        self.setPalette(palette)
        self.setWindowIcon(QIcon(":/kmer.png"))
        self.setWindowTitle("K-mer Matrix of Sequences")
        desktop = QApplication.desktop()
        des_w = desktop.width()
        des_h = desktop.height()
        win_w = des_w * 0.35
        win_h = des_h * 0.3
        self.resize(int(win_w), int(win_h))
        self.initUI()

    def initUI(self):
        self.setToolTip("提示：<b>此功能用于获取输入序列集的k-mer矩阵</b>。注：输入序列无需对齐，且软件会自动剔除序列集中的gap(-)。")
        self.file_edit = QLineEdit_get()
        self.file_edit.setToolTip("Drag a <b>sequence-file with fasta format</b> and drop here")
        self.file_button = InputFile_Button()
        self.file_button.setToolTip("Or click to load a <b>sequence-file with fasta format</b>")
        self.file_button.clicked.connect(self.openfile)
        fasfile_box = MyGroupBox("&Input file", self)
        fasfile_box.setMaximumHeight(100)
        hfasfile_box = QHBoxLayout()
        hfasfile_box.addWidget(self.file_edit)
        hfasfile_box.addWidget(self.file_button)
        fasfile_box.setLayout(hfasfile_box)
        kmer_label = QLabel("Value of K:")
        self.k_value_edit = QLineEdit("4")
        data_type_label = QLabel("Data type:")
        self.data_type_select = QComboBox()
        self.data_type_select.addItem("DNA")
        self.data_type_select.addItem("Amino acid")
        self.data_type_select.setToolTip("<b>DNA</b>: A, C, G, T; <b>Amino acid</b>: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y")
        part_load_label = QLabel("Max.num.seqs in iterative: ", self)
        self.max_seq_nums = QLineEdit("20000", self)
        self.max_seq_nums.setToolTip("The strategy of loading sequences in segments. <b>Smaller value consumes less memory, but increases some time in loading sequences.</b>")
        thread_label = QLabel("Threads used: ", self)
        self.align_thread = QComboBox(self)
        for x in range(0, int(multiprocessing.cpu_count())):
            self.align_thread.addItem(str(x + 1))

        self.align_thread.setToolTip("Larger value of thread consume more CPU, but shortens the computation time.")
        para_box = MyGroupBox("&Parameter", self)
        para_h = QHBoxLayout()
        para_h.addWidget(kmer_label)
        para_h.addWidget(self.k_value_edit)
        para_h.addWidget(data_type_label)
        para_h.addWidget(self.data_type_select)
        parax_Hlayout = QHBoxLayout()
        parax_Hlayout.addWidget(part_load_label)
        parax_Hlayout.addWidget(self.max_seq_nums)
        parax_Hlayout.addWidget(thread_label)
        parax_Hlayout.addWidget(self.align_thread)
        para_layout = QVBoxLayout()
        para_layout.addLayout(para_h)
        para_layout.addLayout(parax_Hlayout)
        para_box.setLayout(para_layout)
        self.progressBar = QProgressBar_beautify()
        self.progressBar.setValue(0)
        self.start_button = Start_Button()
        self.start_button.clicked.connect(self.Work)
        progressBar_hbox = QHBoxLayout()
        progressBar_hbox.addWidget(self.progressBar)
        progressBar_hbox.addWidget(self.start_button)
        progressBar_box = MyGroupBox("&Run and Progress", self)
        progressBar_box.setLayout(progressBar_hbox)
        win_V_layout = QVBoxLayout()
        win_V_layout.setSpacing(20)
        win_V_layout.addWidget(fasfile_box)
        win_V_layout.addWidget(para_box)
        win_V_layout.addWidget(progressBar_box)
        self.setLayout(win_V_layout)

    def openfile(self):
        """
        按钮加载文件路径
        :return:
        """
        fasname = QFileDialog.getOpenFileName(self, "Select a sequence file", "/home")
        fasfile_path = fasname[0]
        if fasfile_path != "":
            self.file_edit.setText(fasfile_path)

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

    def run_finshed(self, fassave_dir):
        """
        运行结束时弹窗
        :param fassave_filepath: 储存文件路径
        :return:
        """
        PopTipWindow.success("Successfully produced!", "<b>Saved in: </b>" + fassave_dir)

    def closeEvent(self, QCloseEvent):
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
                self.thread_kmer.terminate()
            except:
                pass

            try:
                self.error.close()
            except:
                pass

        self.close()

    def Work(self):
        self.progressBar.setValue(0)
        seqfile_path = self.file_edit.text().replace("\\", "/")
        if seqfile_path == "" or os.path.isdir(seqfile_path):
            return
        k_value = self.k_value_edit.text()
        data_types = self.data_type_select.currentText()
        thread_number = self.align_thread.currentText().strip()
        max_num_seq = self.max_seq_nums.text().strip()
        self.start_button.setEnabled(False)
        self.start_button.setText("Run...")
        self.thread_kmer = RunThread_kmer(seqfile_path, k_value, data_types, thread_number, max_num_seq)
        self.thread_kmer.single_progress.connect(self.set_progressbar_value)
        self.thread_kmer.single_save.connect(self.run_finshed)
        self.thread_kmer.start()


class RunThread_kmer(QThread):
    single_progress = Signal(int)
    single_save = Signal(str)

    def __init__(self, seqfile_path, k_value, data_types, thread_number, max_num_seq):
        """

         :param seqfile_path: 序列文件路径
         :param k_value: k-mer大小
         :param data_types: 数据类型
         """
        super(RunThread_kmer, self).__init__()
        self.seqfile_path = seqfile_path
        self.kmer_value = int(k_value)
        self.data_types = data_types
        self.thread_number = int(thread_number)
        self.max_num_seq = int(max_num_seq)

    def run(self):
        try:
            try:
                self.single_progress.emit(20)
                input_data_dir, out_prefix = resolve_file_path(self.seqfile_path)
                out_matrix_path = input_data_dir + "/" + out_prefix + "kmer_matrix.txt"
                str_kind = []
                if self.data_types == "Amino acid":
                    str_kind = [
                     'A', 'C', 'D', 'E', 
                     'F', 'G', 'H', 'I', 'K', 'L', 
                     'M', 
                     'N', 'P', 'Q', 'R', 'S', 'T', 'V', 
                     'W', 'Y']
                else:
                    if self.data_types == "DNA":
                        str_kind = [
                         "A", "C", "G", "T"]
                all_kmer_type = reduce((lambda x, y: [i + j for i in x for j in iter(y)]), [str_kind] * self.kmer_value)
                colum_name_list = [
                 "Seqname"]
                for each_kind in all_kmer_type:
                    colum_name_list.append(each_kind)

                first_line = "\t".join(colum_name_list) + "\n"
                out_matrix = open(out_matrix_path, "a", encoding="utf-8")
                out_matrix.write(first_line)
                run_flage = True
                run_nums = 0
                while run_flage:
                    run_nums += 1
                    all_seq_num, seq_list, run_flage = piecewise_read_seq(self.seqfile_path, self.max_num_seq, run_nums)
                    kemer_pool = Pool(self.thread_number)
                    seq_kmer_result = []
                    for n in range(len(seq_list)):
                        seq_name = seq_list[n][0]
                        seq_contain = seq_list[n][1].replace("-", "")
                        seq_kmer = kemer_pool.apply_async(calculate_kmers, args=(
                         seq_name,
                         seq_contain,
                         self.kmer_value,
                         all_kmer_type))
                        seq_kmer_result.append(seq_kmer)

                    kemer_pool.close()
                    kemer_pool.join()
                    for each_seq in seq_kmer_result:
                        out_matrix.write(each_seq.get() + "\n")

                    jindu = run_nums / (int(all_seq_num / self.max_num_seq) + 1)
                    self.single_progress.emit(45 + int(40 * jindu))

                self.single_progress.emit(100)
                self.single_save.emit(out_matrix_path)
            except Exception as e:
                try:
                    config.run_error_log = "traceback.format_exc():\n%s" % traceback.format_exc()
                    self.single_progress.emit(1314)
                finally:
                    e = None
                    del e

        finally:
            self.exit(0)
