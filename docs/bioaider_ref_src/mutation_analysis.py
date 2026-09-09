# uncompyle6 version 3.9.3
# Python bytecode version base 3.7.0 (3394)
# Decompiled from: Python 3.12.10 (tags/v3.12.10:0cc8128, Apr  8 2025, 12:21:36) [MSC v.1943 64 bit (AMD64)]
# Embedded file name: mutation_tools/mutation_analysis.py
"""
Author: Zhou Zhi-Jian
update Time: 2023/2/14 16:58

"""
import os, datetime, traceback, config
from math import ceil
import multiprocessing
from matplotlib import rcParams
from multiprocessing import Pool
import webbrowser as web, platform
from PySide2.QtWidgets import QApplication, QWidget, QPushButton, QPlainTextEdit, QVBoxLayout, QRadioButton, QComboBox, QHBoxLayout, QFileDialog, QCheckBox, QLabel, QMessageBox, QLineEdit
from PySide2.QtGui import QIcon, QPalette, QCursor, QFont
from PySide2.QtCore import Qt, QThread, Signal
from custom_libraries.gui_class import QLineEdit_get, PopTipWindow, ErrorLogGUI, QProgressBar_beautify, Start_Button, InputFile_Button, MyGroupBox
from custom_libraries.public_functions import aa_kind_dis, mutation_index, make_dir, resolve_file_path
from custom_libraries.codon_trans import DNA_genetic_godes_transle
from mutation_tools.mua_fuc import link_mutation_run1, link_mutation_run2, single_mutation_nt, single_mutation_codon, single_mutation_aa, nuc_aa_plot, codon_plot, piecewise_read, lollipop_data, save_codon_log
import psutil

def calculate_memory():
    pid = os.getpid()
    p = psutil.Process(pid)
    info = p.memory_full_info()
    memory = info.uss / 1024 / 1024
    return memory


class MutationAnalysis(QWidget):
    __doc__ = "\n    突变分析GUI\n    "

    def __init__(self, parent=None):
        super(MutationAnalysis, self).__init__(parent)
        self.setAcceptDrops(True)
        desktop = QApplication.desktop()
        des_w = desktop.width()
        des_h = desktop.height()
        palette = QPalette()
        palette.setColor(QPalette.Background, Qt.white)
        self.setPalette(palette)
        self.setWindowIcon(QIcon(":/M.png"))
        self.setWindowTitle("Mutation Analysis")
        win_w = des_w * 0.3
        win_h = des_h * 0.6
        self.resize(int(win_w), int(win_h))
        self.initUI()

    def initUI(self):
        self.setToolTip("提示：<b>此功能用于序列对齐后的突变分析</b>。将序列集中第一条序列视为参考序列，计算突变位点、突变类型和突变频率等。<b>支持核苷酸、编码基因和氨基酸序列数据，支持多核并行计算。</b>")
        self.file_edit1 = QLineEdit_get()
        self.file_button = InputFile_Button()
        self.file_button.setToolTip("Or click to load a <b>alignment sequence-file with fasta format</b>")
        self.file_button.clicked.connect(self.openfile)
        las_tip = QLabel("<b>Tip: </b>The <b>first sequence in the data set</b> will be used as the <b>reference sequence</b> for mutation analysis", self)
        las_tip.setTextInteractionFlags(Qt.TextSelectableByMouse)
        las_tip.setWordWrap(True)
        file_box1 = MyGroupBox("&Input file", self)
        file_box1.setToolTip("Drag a <b> alignment sequence-file with fasta format </b> and drop here")
        file_hbox1 = QHBoxLayout()
        file_hbox1.addWidget(self.file_edit1)
        file_hbox1.addWidget(self.file_button)
        las_H = QHBoxLayout()
        las_H.addWidget(las_tip)
        file_box_labs = QVBoxLayout()
        file_box_labs.addLayout(file_hbox1)
        file_box_labs.addLayout(las_H)
        file_box1.setLayout(file_box_labs)
        mode_box = MyGroupBox("&Choose an analysis method", self)
        self.singele_mode = QRadioButton("Single Mutation")
        self.singele_mode.clicked.connect(self.selcet_singmode)
        self.singele_mode.setChecked(True)
        self.linked_mode = QRadioButton("Linked (Multiply) Mutation")
        self.linked_mode.setToolTip("<b>Tip:</b> List all mutations (including insertion and deletion) of each sequence compared with the reference sequence.")
        self.linked_mode.clicked.connect(self.selcet_mulmode)
        modelH = QHBoxLayout()
        modelH.addWidget(self.singele_mode)
        modelH.addWidget(self.linked_mode)
        mode_box.setLayout(modelH)
        self.parameter_box = MyGroupBox("&Parameter", self)
        label_data = QLabel("Datas type: ", self)
        self.Radio_nucleotide = QRadioButton("Nucleotide")
        self.Radio_nucleotide.setChecked(True)
        self.Radio_nucleotide.clicked.connect(self.selcet_nucleotide)
        self.Radio_codon = QRadioButton("Codon")
        self.Radio_codon.clicked.connect(self.selcet_codon)
        self.Radio_codon.setToolTip("The sequence sets should <b>be aligned based on codon method beforehand</b> if you select the <b>'Codon'</b>.")
        self.Radio_protein = QRadioButton("Protein")
        self.Radio_protein.clicked.connect(self.selcet_protein)
        label_codon = QLabel("Codon Table: ", self)
        self.selcet = QComboBox(self)
        self.selcet.addItem("1.The Standard Code")
        self.selcet.addItem("2.The Vertebrate Mitochondrial Code")
        self.selcet.addItem("3.The Yeast Mitochondrial Code")
        self.selcet.addItem("4.The Mold, Protozoan,... and the Mycoplasma/Spiroplasma Code")
        self.selcet.addItem("5.The Invertebrate Mitochondrial Code")
        self.selcet.addItem("6.The Ciliate, Dasycladacean and Hexamita Nuclear Code")
        self.selcet.addItem("9.The Echinoderm and Flatworm Mitochondrial Code")
        self.selcet.addItem("10.The Euplotid Nuclear Code")
        self.selcet.addItem("11.The Bacterial, Archaeal and Plant Plastid Code")
        self.selcet.addItem("12.The Alternative Yeast Nuclear Code")
        self.selcet.addItem("13.The Ascidian Mitochondrial Code")
        self.selcet.addItem("14.The Alternative Flatworm Mitochondrial Code")
        self.selcet.addItem("16.Chlorophycean Mitochondrial Code")
        self.selcet.addItem("21.Trematode Mitochondrial Code")
        self.selcet.addItem("22.Scenedesmus obliquus Mitochondrial Code")
        self.selcet.addItem("23.Thraustochytrium Mitochondrial Code")
        self.selcet.addItem("24.Rhabdopleuridae Mitochondrial Code")
        self.selcet.addItem("25.Candidate Division SR1 and Gracilibacteria Code")
        self.selcet.addItem("26.Pachysolen tannophilus Nuclear Code")
        self.selcet.addItem("27.Karyorelict Nuclear Code")
        self.selcet.addItem("28.Condylostoma Nuclear Code")
        self.selcet.addItem("29.Mesodinium Nuclear Code")
        self.selcet.addItem("30.Peritrich Nuclear Code")
        self.selcet.addItem("31.Blastocrithidia Nuclear Code")
        self.selcet.addItem("33.Cephalodiscidae Mitochondrial UAA-Tyr Code")
        self.selcet.setEnabled(False)
        doubt_button = QPushButton("", self)
        doubt_button.setToolTip("Click to learn more about the genetic codes")
        doubt_button.setStyleSheet("QPushButton{border-image: url(:/doubt1.png)}QPushButton:hover{border-image: url(:/doubt2.png)}QPushButton:pressed{border-image: url(:/ico/doubt2.png)}")
        doubt_button.setMinimumSize(25, 24)
        doubt_button.setMaximumSize(25, 24)
        doubt_button.setCursor(QCursor(Qt.PointingHandCursor))
        doubt_button.clicked.connect(self.openUrl)
        self.com_selc = QCheckBox("Output substitution frequency distribution", self)
        self.com_selc.stateChanged.connect(self.selc_choose)
        self.sy_nonsy_delte = QCheckBox("Delete both Sys-Nonsys nt sites", self)
        self.sy_nonsy_delte.setChecked(False)
        self.sy_nonsy_delte.setEnabled(False)
        self.sy_nonsy_delte.setToolTip("Whether to delete the nt sites with both synonymous and non-synonymous substitutions?")
        part_load_label = QLabel("Max.num.seqs in iterative: ", self)
        self.max_seq_nums = QLineEdit("20000", self)
        self.max_seq_nums.setToolTip("The strategy of loading sequences in segments. <b>Smaller value consumes less memory, but increases some time in loading sequences.</b>")
        thread_label = QLabel("Threads used: ", self)
        self.align_thread = QComboBox(self)
        for x in range(0, int(multiprocessing.cpu_count())):
            self.align_thread.addItem(str(x + 1))

        self.align_thread.setToolTip("Larger value of thread consume more CPU, but shortens the computation time.")
        self.save_log = QCheckBox("Save the run logs?", self)
        self.save_log.setChecked(True)
        self.save_log.setToolTip("It will take extra time to generate run logs.")
        data_type_H = QHBoxLayout()
        data_type_H.addWidget(label_data)
        data_type_H.addWidget(self.Radio_nucleotide)
        data_type_H.addWidget(self.Radio_codon)
        data_type_H.addWidget(self.Radio_protein)
        H_condon = QHBoxLayout()
        H_condon.addWidget(label_codon)
        H_condon.addWidget(self.selcet)
        H_condon.addWidget(doubt_button)
        parax_Hlayout = QHBoxLayout()
        parax_Hlayout.addWidget(thread_label)
        parax_Hlayout.addWidget(self.align_thread)
        parax_Hlayout.addWidget(part_load_label)
        parax_Hlayout.addWidget(self.max_seq_nums)
        parax_Hlayout.addWidget(self.save_log)
        parameter_V_box = QVBoxLayout()
        parameter_V_box.addLayout(data_type_H)
        parameter_V_box.addLayout(H_condon)
        parameter_V_box.addLayout(parax_Hlayout)
        self.parameter_box.setLayout(parameter_V_box)
        com_selc_hbox = QHBoxLayout()
        com_selc_hbox.addWidget(self.com_selc)
        com_selc_hbox.addWidget(self.sy_nonsy_delte)
        self.taxon = QPlainTextEdit(self)
        self.taxon.setFont(QFont(config.app_font, config.app_font_size))
        self.taxon.setEnabled(False)
        self.taxon.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.taxon.setStyleSheet("border: 1px solid #8B8B7A;border-radius:2px")
        self.taxon.setToolTip("Please specify the groups of substitution frequency. <b>Tip: </b>Insertions, deletions and degenerate bases are not used to plot.")
        self.taxon_hbox = QHBoxLayout()
        self.taxon_hbox.addWidget(self.taxon)
        self.taxon_box = MyGroupBox("&Plot substitution frequency", self)
        self.taxon_vbox = QVBoxLayout()
        self.taxon_vbox.addLayout(com_selc_hbox)
        self.taxon_vbox.addLayout(self.taxon_hbox)
        self.taxon_box.setLayout(self.taxon_vbox)
        self.progressBar = QProgressBar_beautify()
        self.progressBar.setValue(0)
        self.make_button = Start_Button()
        self.make_button.clicked.connect(self.Work)
        progressBar_hbox = QHBoxLayout()
        progressBar_hbox.addWidget(self.progressBar)
        progressBar_hbox.addWidget(self.make_button)
        progressBar_box = MyGroupBox("&Run and Progress", self)
        progressBar_box.setLayout(progressBar_hbox)
        win_H_layout = QVBoxLayout()
        win_H_layout.setSpacing(20)
        win_H_layout.addWidget(file_box1)
        win_H_layout.addWidget(mode_box)
        win_H_layout.addWidget(self.parameter_box)
        win_H_layout.addWidget(self.taxon_box)
        win_H_layout.addWidget(progressBar_box)
        self.setLayout(win_H_layout)

    def selcet_singmode(self):
        self.com_selc.setEnabled(True)

    def selcet_mulmode(self):
        self.selcet.setEnabled(False)
        self.sy_nonsy_delte.setChecked(False)
        self.sy_nonsy_delte.setEnabled(False)
        self.com_selc.setChecked(False)
        self.com_selc.setEnabled(False)

    def selc_choose(self, state):
        """
        是否生成突变频率分布图 复选框状态
        :param state: 复选框状态
        :return:
        """
        if state == Qt.Checked:
            self.taxon.setEnabled(True)
            if self.Radio_codon.isChecked():
                self.sy_nonsy_delte.setEnabled(True)
        elif state == Qt.Unchecked:
            self.taxon.setEnabled(False)
            self.sy_nonsy_delte.setEnabled(False)
            self.sy_nonsy_delte.setChecked(False)

    def selcet_codon(self):
        """
        选择编码基因
        :return:
        """
        if self.linked_mode.isChecked():
            return
        self.selcet.setEnabled(True)
        if self.com_selc.isChecked():
            self.sy_nonsy_delte.setEnabled(True)

    def selcet_nucleotide(self):
        """
        选择核苷酸
        :return:

        """
        self.selcet.setEnabled(False)
        self.sy_nonsy_delte.setEnabled(False)
        self.sy_nonsy_delte.setChecked(False)
        self.taxon_box.setToolTip("Please pasta the frequency groups of substitution sites (<b>Tip: </b>Insertions, deletions and degenerate bases are not considered for calculation!)")

    def selcet_protein(self):
        """
        选择蛋白质序列
        :return:
        """
        self.selcet.setEnabled(False)
        self.sy_nonsy_delte.setEnabled(False)
        self.sy_nonsy_delte.setChecked(False)
        self.taxon_box.setToolTip("Please pasta the frequency groups of substitution sites (<b>Tip: </b>Insertions, deletions, terminator and unknown amino acid are not considered for calculation!)")

    def openUrl(self):
        """
        打开网页链接
        :return:
        """
        web.open_new_tab("https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG11")

    def Work(self):
        """
        运行前检查和参数设置
        :return:
        """
        file_path = self.file_edit1.text()
        file_path = file_path.replace("\\", "/")
        if file_path == "" or os.path.isdir(file_path):
            return
        para_dic = {}
        para_dic["input_file_path"] = file_path
        if self.singele_mode.isChecked():
            para_dic["work_modes"] = "single_modes"
        else:
            if self.linked_mode.isChecked():
                para_dic["work_modes"] = "linked_modes"
            else:
                para_dic["thread_num"] = self.align_thread.currentText()
                para_dic["max_num_seq"] = int(self.max_seq_nums.text().strip())
                if self.save_log.isChecked():
                    para_dic["save_log"] = True
                else:
                    para_dic["save_log"] = False
            para_dic["codon_table"] = ""
        if self.Radio_nucleotide.isChecked():
            para_dic["muta_data_type"] = "nucleotide"
        else:
            if self.Radio_codon.isChecked():
                para_dic["muta_data_type"] = "codon"
                para_dic["codon_table"] = self.selcet.currentText()
            else:
                if self.Radio_protein.isChecked():
                    para_dic["muta_data_type"] = "protein"
                else:
                    para_dic["group_tab"] = []
                    if self.com_selc.isChecked():
                        if self.taxon.toPlainText() == "":
                            return
                            para_dic["draw_frequency_distribution"] = "Yes"
                            para_dic["group_tab"] = self.taxon.toPlainText().strip().split("\n")
                        else:
                            para_dic["draw_frequency_distribution"] = "No"
                        if self.sy_nonsy_delte.isChecked():
                            para_dic["sy_nonsy_select_value"] = "True"
                    else:
                        para_dic["sy_nonsy_select_value"] = "Flase"
                self.make_button.setEnabled(False)
                self.make_button.setText("Run...")
                self.thread_e = RunThread_e(para_dic)
                self.thread_e.trigger.connect(self.set_progressbar_value)
                self.thread_e.trigger_save.connect(self.run_finshed)
                self.thread_e.start()

    def openfile(self):
        """
        按钮加载文件路径
        :return:
        """
        self.fname = QFileDialog.getOpenFileName(self, "Select an aligned sequence file", "/home")
        input_file_path = self.fname[0]
        if input_file_path != "":
            self.file_edit1.setText(input_file_path)

    def set_progressbar_value(self, value):
        """
        进度条,如果收到“1314”信号，则弹出error log窗口
        :param value: 进度条值
        :return:
        """
        self.progressBar.setValue(value)
        if value == 100:
            self.make_button.setEnabled(True)
            self.make_button.setText("Start")
            return
        if value == 1314:
            self.error = ErrorLogGUI()
            self.error.show()
            self.make_button.setEnabled(True)
            self.make_button.setText("Start")
            return

    def run_finshed(self, savesdir):
        """
        运行结束时弹窗
        :param savesdir: 储存文件路径
        :return:
        """
        PopTipWindow.success("Successfully produced!", "<b>Saved in: </b>" + savesdir)

    def closeEvent(self, QCloseEvent):
        """
        重写窗口关闭事件
        :param QCloseEvent:
        :return:
        """
        if self.make_button.text() == "Run...":
            a = QMessageBox.question(self, "Exit?", "The task is running. Are you sure to exit?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if a == QMessageBox.Yes:
                pass
            else:
                QCloseEvent.ignore()
                return
        self.progressBar.setValue(0)
        self.make_button.setEnabled(True)
        self.make_button.setText("Start")
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


class RunThread_e(QThread):
    __doc__ = "\n    突变分析运行线程\n    "
    trigger = Signal(int)
    trigger_save = Signal(str)

    def __init__(self, para_dic):
        """
        :param file_path: 输入文件路径
        :param  para_dic: 参数字典
        """
        super(RunThread_e, self).__init__()
        self.para_dic = para_dic
        rcParams["font.sans-serif"] = [
         config.app_font]

    def run(self):
        start_memory = calculate_memory()
        try:
            try:
                input_file_path = self.para_dic["input_file_path"]
                work_modes = self.para_dic["work_modes"]
                muta_data_type = self.para_dic["muta_data_type"]
                codon_table = self.para_dic["codon_table"]
                draw_frequency_dis = self.para_dic["draw_frequency_distribution"]
                frequency_group_tab = self.para_dic["group_tab"]
                sy_nonsy_delte_value = self.para_dic["sy_nonsy_select_value"]
                input_data_dir, out_prefix = resolve_file_path(input_file_path)
                now = datetime.datetime.now()
                time_now = now.strftime("%Y-%m-%d %H:%M:%S")
                run_id = time_now.replace("-", "").replace(" ", "_").replace(":", "")
                saves_dir = input_data_dir + "/" + "Mutation_analysis_" + run_id
                make_dir(saves_dir)
                out_prefix = out_prefix
                savefile_path = saves_dir + "/" + out_prefix
                self.trigger.emit(40)
                if work_modes == "single_modes":
                    if muta_data_type == "nucleotide":
                        run_flage = True
                        run_nums = 0
                        max_num_seq = self.para_dic["max_num_seq"]
                        seq_lenth = 0
                        qury_seq_name = ""
                        qury_seq = ""
                        mutation_count_dic = {}
                        mutation_jx_count = {}
                        run_log_list = []
                        while run_flage:
                            run_nums += 1
                            all_seq_num, seq_tab, run_flage = piecewise_read(input_file_path, max_num_seq, run_nums)
                            seq_count = len(seq_tab)
                            run_pool = Pool(int(self.para_dic["thread_num"]))
                            mutation_results_list = []
                            for n in range(seq_count):
                                seq_index = seq_tab[n][0]
                                seq_name = seq_tab[n][1]
                                seq_contain = seq_tab[n][2].replace("U", "T").strip()
                                if seq_index == 1:
                                    qury_seq_name = seq_name
                                    qury_seq = seq_contain
                                    seq_lenth = len(qury_seq)
                                else:
                                    mutation_results = run_pool.apply_async(single_mutation_nt, args=(
                                     seq_name, qury_seq, seq_contain, seq_lenth))
                                    mutation_results_list.append(mutation_results)

                            run_pool.close()
                            run_pool.join()
                            for each_seq in mutation_results_list:
                                seq_result = each_seq.get()
                                site_mution_list = seq_result[1]
                                site_mution_JX_list = seq_result[2]
                                if site_mution_list != []:
                                    for site in site_mution_list:
                                        new_site = site[1] + "\t" + "".join(site)
                                        if new_site not in mutation_count_dic:
                                            mutation_count_dic[new_site] = 1
                                        else:
                                            mutation_count_dic[new_site] += 1

                                if site_mution_JX_list != []:
                                    for site in site_mution_JX_list:
                                        new_site = site[1] + "\t" + "".join(site)
                                        if new_site not in mutation_jx_count:
                                            mutation_jx_count[new_site] = 1
                                        else:
                                            mutation_jx_count[new_site] += 1

                            if self.para_dic["save_log"]:
                                for each_seq in mutation_results_list:
                                    seq_result = each_seq.get()
                                    seq_name = seq_result[0]
                                    site_mution_list = seq_result[1]
                                    if site_mution_list != []:
                                        seq_mution = []
                                        for each_site in site_mution_list:
                                            seq_mution.append("".join(each_site))

                                        run_log_list.append(seq_name + "\t" + ", ".join(seq_mution))

                            jindu = run_nums / ceil(all_seq_num / max_num_seq)
                            self.trigger.emit(45 + int(40 * jindu))

                        self.trigger.emit(90)
                        if run_log_list:
                            with open((savefile_path + "_run_logs.txt"), "w", encoding="utf-8") as log_txt:
                                log_txt.write("Reference sequence: " + qury_seq_name + "\n\n" + "sequence name" + "\t" + "mutation nt" + "\n")
                                for each_record in run_log_list:
                                    log_txt.write(each_record + "\n")

                        self.trigger.emit(95)
                        summary_txt = open((savefile_path + "_mutation_site_summary.txt"),
                          "w", encoding="utf-8")
                        summary_txt.write("Reference sequence: " + qury_seq_name + "\n\n" + "Site (nt)" + "\t" + "mutation" + "\t" + "sequence count" + "\n")
                        for key in mutation_count_dic:
                            summary_txt.write(key + "\t" + str(mutation_count_dic[key]) + "\n")

                        summary_txt.close()
                        sus_site_no_others_count = open((savefile_path + "_mutation_site_(not_gap_degenerate_bases).txt"),
                          "w",
                          encoding="utf-8")
                        sus_site_no_others_count.write("Reference sequence: " + qury_seq_name + "\n\n" + "Site (nt)" + "\t" + "mutation" + "\t" + "sequence count" + "\n")
                        for key in mutation_jx_count:
                            sus_site_no_others_count.write(key + "\t" + str(mutation_jx_count[key]) + "\n")

                        sus_site_no_others_count.close()
                        self.trigger.emit(95)
                        if draw_frequency_dis == "Yes":
                            site_fru = {}
                            for key in mutation_jx_count:
                                site = key.split("\t")[0]
                                if site not in site_fru:
                                    site_fru[site] = mutation_jx_count[key]
                                else:
                                    site_fru[site] += mutation_jx_count[key]

                            nuc_aa_plot(site_fru, frequency_group_tab, savefile_path)
                elif muta_data_type == "codon":
                    codons_statistics_dir = saves_dir + "/" + "Statistics_in_codon"
                    base_statistics_dir = saves_dir + "/" + "Statistics_in_base(nt)"
                    make_dir(codons_statistics_dir)
                    make_dir(base_statistics_dir)
                    translate = DNA_genetic_godes_transle
                    d_site = []
                    self.trigger.emit(50)
                    mutation_summary_dic = {}
                    mutation_type_dic = {}
                    run_log_list = []
                    run_flage = True
                    run_nums = 0
                    max_num_seq = self.para_dic["max_num_seq"]
                    all_seq_count = 0
                    codon_count = 0
                    qury_seq_name = ""
                    qury_seq = ""
                    while run_flage:
                        run_nums += 1
                        all_seq_num, seq_tab, run_flage = piecewise_read(input_file_path, max_num_seq, run_nums)
                        seq_count = len(seq_tab)
                        all_seq_count += seq_count
                        run_pool = Pool(int(self.para_dic["thread_num"]))
                        mutation_results_list = []
                        for n in range(len(seq_tab)):
                            seq_index = seq_tab[n][0]
                            seq_name = seq_tab[n][1]
                            seq = seq_tab[n][2].replace("U", "T").strip()
                            if seq_index == 1:
                                qury_seq_name = seq_name
                                qury_seq = seq
                                codon_count = int(len(qury_seq) / 3)
                            else:
                                mutation_results = run_pool.apply_async(single_mutation_codon, args=(
                                 qury_seq, seq_name, seq, codon_count,
                                 codon_table, translate))
                                mutation_results_list.append(mutation_results)

                        run_pool.close()
                        run_pool.join()
                        jindu = run_nums / ceil(all_seq_num / max_num_seq)
                        self.trigger.emit(45 + int(40 * jindu))
                        for each_seq in mutation_results_list:
                            seq_result = each_seq.get()
                            mutation_list = seq_result[1]
                            if mutation_list != []:
                                for each_mutation in mutation_list:
                                    if each_mutation not in mutation_summary_dic:
                                        mutation_summary_dic[each_mutation] = 1
                                    else:
                                        mutation_summary_dic[each_mutation] += 1
                                    mutation_codon_list = each_mutation.split("\t")
                                    site_mutation_type = "\t".join([mutation_codon_list[1],
                                     mutation_codon_list[3]])
                                    if site_mutation_type not in mutation_type_dic:
                                        mutation_type_dic[site_mutation_type] = 1
                                    else:
                                        mutation_type_dic[site_mutation_type] += 1

                        if self.para_dic["save_log"]:
                            for each_seq in mutation_results_list:
                                seq_result = each_seq.get()
                                seq_name = seq_result[0]
                                mutation_list = seq_result[1]
                                if mutation_list != []:
                                    for each_mutation in mutation_list:
                                        mutation_codon_list = each_mutation.split("\t")
                                        mutation_nt_location = mutation_index(mutation_codon_list[0], mutation_codon_list[2])
                                        mutation_str = "".join([
                                         mutation_codon_list[0], " to ", mutation_codon_list[2],
                                         " (", mutation_codon_list[-2], mutation_codon_list[1],
                                         mutation_codon_list[-1], ")"])
                                        logs = "\t".join([seq_name, mutation_codon_list[1],
                                         mutation_nt_location,
                                         mutation_codon_list[3],
                                         mutation_str])
                                        run_log_list.append(logs)

                    if run_log_list:
                        with open((savefile_path + "_run_logs.txt"), "w", encoding="utf-8") as log_txt:
                            log_txt.write("Reference sequence: " + qury_seq_name + "\n\n" + "Sequence name" + "\t" + " Site (codon)" + "\t" + "Index of nt in codon" + "\t" + "mutation type" + "\t" + "mutation" + "\n")
                            for each_record in run_log_list:
                                log_txt.write(each_record + "\n")

                    mutation_type_list = [
                     'Synonymous', 'Nonsynonymous', 
                     'Termination', 
                     'Insertion', 'Deletion', 'Unknown']
                    for site in range(codon_count):
                        d_site.append("codon_site" + str(site + 1))

                    mutation_site_list = []
                    for key in mutation_type_dic.keys():
                        codon_site, mu_type = key.split("\t")
                        if int(codon_site) not in mutation_site_list:
                            mutation_site_list.append(int(codon_site))

                    mutation_site_list.sort()
                    with open(codons_statistics_dir + "/" + out_prefix + "_mutation_codon_site.csv", "w") as mutation_count_file:
                        mutation_count_file.write("Site (codon)," + ",".join(mutation_type_list) + "\n")
                        for site in mutation_site_list:
                            site = str(site)
                            site_summary_list = []
                            for mutation_type in mutation_type_list:
                                query_key = site + "\t" + mutation_type
                                if query_key in mutation_type_dic:
                                    site_summary_list.append(str(mutation_type_dic[query_key]))
                                else:
                                    site_summary_list.append("0")

                            mutation_count_file.write(site + "," + ",".join(site_summary_list) + "\n")

                    concise_infor_list = []
                    lollipop_infro = []
                    summary_txt = open((savefile_path + "_mutation_site_summary.txt"), "w",
                      encoding="utf-8")
                    summary_txt.write("Reference sequence: " + qury_seq_name + "\n\n" + "Site (codon)" + "\t" + "Index of mutation nt" + "\t" + "mutation type" + "\t" + "mutation" + "\t" + "changed in properties of aa?" + "\t" + "changed type of properties" + "\t" + "sequence count" + "\n")
                    for key in mutation_summary_dic:
                        muta_list = key.split("\t")
                        mutation_nt_location = mutation_index(muta_list[0], muta_list[2])
                        tik = "\t"
                        if muta_list[3] == "Nonsynonymous":
                            tik = aa_kind_dis(muta_list[-2], muta_list[-1])
                        mutation_str = "".join([muta_list[0], " to ", muta_list[2],
                         " (", muta_list[-2], muta_list[1], muta_list[-1], ")"])
                        numbs = str(mutation_summary_dic[key])
                        concise_infor_list.append([muta_list[1], mutation_nt_location, muta_list[3], numbs])
                        lollipop_infro.append([muta_list[1], mutation_nt_location,
                         muta_list[3], mutation_str, numbs])
                        summary_txt.write("\t".join([muta_list[1], mutation_nt_location,
                         muta_list[3], mutation_str, tik, numbs]) + "\n")

                    self.trigger.emit(90)
                    singe_nt_site = {}
                    two_three_musite = []
                    stop_nt_site = {}
                    Synonymous_nt_site = {}
                    Nonsynonymous_nt_site = {}
                    self.trigger.emit(95)
                    for line in concise_infor_list:
                        condons, base_inx, kinds, seqcou = line
                        if kinds == "Synonymous" or kinds == "Nonsynonymous":
                            if base_inx.find(",") != -1:
                                two_three_musite.append("\t".join(line))
                            else:
                                key = condons + "-" + base_inx + "\t" + kinds
                                if key not in singe_nt_site:
                                    singe_nt_site[key] = int(seqcou)
                                else:
                                    singe_nt_site[key] += int(seqcou)
                        elif kinds == "Termination":
                            if base_inx.find(",") != -1:
                                d_cout = base_inx.count(",")
                                bb = base_inx.split(",")
                                for xx in range(d_cout + 1):
                                    key = condons + "-" + bb[xx] + "\t" + kinds
                                    if key not in stop_nt_site:
                                        stop_nt_site[key] = int(seqcou)
                                    else:
                                        stop_nt_site[key] += int(seqcou)

                            else:
                                key = condons + "-" + base_inx + "\t" + kinds
                                if key not in stop_nt_site:
                                    stop_nt_site[key] = int(seqcou)
                                else:
                                    stop_nt_site[key] += int(seqcou)

                    updata_singe_nt_site = singe_nt_site
                    for mul_subs in two_three_musite:
                        condons, base_inx, kinds, seqcou = mul_subs.split("\t")
                        base_cout = base_inx.count(",") + 1
                        base_list = base_inx.split(",")
                        for m in range(base_cout):
                            condon_nt_kinds = condons + "-" + base_list[m] + "\t" + kinds
                            if condon_nt_kinds in updata_singe_nt_site:
                                updata_singe_nt_site[condon_nt_kinds] += int(seqcou)
                            else:
                                updata_singe_nt_site[condon_nt_kinds] = int(seqcou)

                    for keys in updata_singe_nt_site:
                        if keys.split("\t")[1] == "Synonymous":
                            Synonymous_nt_site[keys] = updata_singe_nt_site[keys]

                    sy_nonsy = 0
                    syn_nonsy_site = {**Synonymous_nt_site, **Nonsynonymous_nt_site}
                    with open((base_statistics_dir + "/" + out_prefix + "_syn_or_non-syn_substitution_nt_sites.csv"),
                      "w", encoding="utf-8") as syn_nonsy_sites:
                        syn_nonsy_sites.write("Site (codon),Index of base in codon,Nucleotide site,Type,Substitution frequency\n")
                        for key in syn_nonsy_site:
                            uf1, uf2 = key.split("\t")
                            codon_site, nt_index = uf1.split("-")
                            nt_num_index = 0
                            if nt_index == "1":
                                nt_num_index = 3 * (int(codon_site) - 1) + 1
                            else:
                                if nt_index == "2":
                                    nt_num_index = 3 * (int(codon_site) - 1) + 2
                                else:
                                    if nt_index == "3":
                                        nt_num_index = 3 * int(codon_site)
                                    syn_nonsy_sites.write(codon_site + "," + nt_index + "," + str(nt_num_index) + "," + uf2 + "," + str(syn_nonsy_site[key]) + "\n")

                    sy_nonsy_biaoji_list = []
                    for sy_key in Synonymous_nt_site:
                        uu1 = sy_key.split("\t")[0]
                        for non_sy_key in Nonsynonymous_nt_site:
                            uu2 = non_sy_key.split("\t")[0]
                            if uu1 == uu2:
                                sy_nonsy += 1
                                sy_nonsy_biaoji_list.append(uu1)
                                if sy_nonsy_delte_value == "True":
                                    del syn_nonsy_site[sy_key]
                                    del syn_nonsy_site[non_sy_key]
                                break

                    sy_nonsy_biaoji = ", ".join(sy_nonsy_biaoji_list)
                    count_sy = len(Synonymous_nt_site) - sy_nonsy
                    count_nonsy = len(Nonsynonymous_nt_site) - sy_nonsy
                    count_stop_mu = len(stop_nt_site)
                    stop_mu_list = []
                    for kes in stop_nt_site:
                        stop_mu_list.append(kes.split("\t")[0])

                    stop_mu_all = ", ".join(stop_mu_list)
                    with open(base_statistics_dir + "/" + out_prefix + "_summary_of_mutant_bases.txt", "w") as nt_summary_txt:
                        only_sy = "Only synonymous mutation\t" + str(count_sy)
                        only_nonsy = "Only non-synonymous mutation\t" + str(count_nonsy)
                        both_sy_nonsy = "Both synonymous and non-synonymous mutation\t" + str(sy_nonsy) + "\t" + "codon-base: " + sy_nonsy_biaoji
                        termi_site = "Termination mutation\t" + str(count_stop_mu) + "\t" + "codon-base: " + stop_mu_all
                        nt_summary_txt.write("sites(nucleotide) count of different types of mutation\n\n")
                        nt_summary_txt.write("Mutation types\tSites count\tDescription\n")
                        nt_summary_txt.write("\n".join([only_sy, only_nonsy, both_sy_nonsy, termi_site]))
                    lollipop_single_nt = lollipop_data(lollipop_infro)
                    with open((savefile_path + "_mutation_site_summary_used_for_lollipop.txt"), "w",
                      encoding="utf-8") as lollipop_summary_txt:
                        lollipop_summary_txt.write("Site(nt)\tType\tMutation\tCount\n")
                        for key in lollipop_single_nt:
                            lollipop_summary_txt.write(key + "\t" + str(lollipop_single_nt[key]) + "\n")

                    if draw_frequency_dis == "Yes":
                        codon_plot(saves_dir, syn_nonsy_site, frequency_group_tab, out_prefix)
                elif muta_data_type == "protein":
                    max_num_seq = self.para_dic["max_num_seq"]
                    run_flage = True
                    run_nums = 0
                    seq_lenth = 0
                    qury_seq_name = ""
                    qury_seq = ""
                    mutation_count_dic = {}
                    mutation_jx_count = {}
                    run_log_list = []
                    while run_flage:
                        run_nums += 1
                        all_seq_num, seq_tab, run_flage = piecewise_read(input_file_path, max_num_seq, run_nums)
                        seq_count = len(seq_tab)
                        run_pool = Pool(int(self.para_dic["thread_num"]))
                        mutation_results_list = []
                        for n in range(seq_count):
                            seq_index = seq_tab[n][0]
                            seq_name = seq_tab[n][1]
                            seq_contain = seq_tab[n][2]
                            if seq_index == 1:
                                qury_seq_name = seq_name
                                qury_seq = seq_contain
                                seq_lenth = len(qury_seq)
                            else:
                                mutation_results = run_pool.apply_async(single_mutation_aa, args=(
                                 qury_seq, seq_lenth, seq_name, seq_contain))
                                mutation_results_list.append(mutation_results)

                        run_pool.close()
                        run_pool.join()
                        jindu = run_nums / ceil(all_seq_num / max_num_seq)
                        self.trigger.emit(45 + int(40 * jindu))
                        for each_seq in mutation_results_list:
                            seq_result = each_seq.get()
                            aa_mutation_list = seq_result[1]
                            aa_mutation_JX_list = seq_result[2]
                            if aa_mutation_list != []:
                                for site in aa_mutation_list:
                                    new_site = site[1] + "\t" + "".join([site[0], site[1], site[2]]) + "\t" + site[3]
                                    if new_site not in mutation_count_dic:
                                        mutation_count_dic[new_site] = 1
                                    else:
                                        mutation_count_dic[new_site] += 1

                            if aa_mutation_JX_list != []:
                                for site in aa_mutation_JX_list:
                                    new_site = site[1] + "\t" + "".join([site[0], site[1], site[2]]) + "\t" + site[3]
                                    if new_site not in mutation_jx_count:
                                        mutation_jx_count[new_site] = 1
                                    else:
                                        mutation_jx_count[new_site] += 1

                        if self.para_dic["save_log"]:
                            for each_seq in mutation_results_list:
                                seq_result = each_seq.get()
                                seq_name = seq_result[0]
                                aa_mutation_list = seq_result[1]
                                if aa_mutation_list != []:
                                    muatation_summary = ["".join([site[0], site[1], site[2]]) for site in aa_mutation_list]
                                    run_log_list.append(seq_name + "\t" + ", ".join(muatation_summary))

                    if run_log_list:
                        with open((savefile_path + "_run_logs.txt"), "w", encoding="utf-8") as log_txt:
                            log_txt.write("Reference sequence: " + qury_seq_name + "\n\n" + "sequence name" + "\t" + "mutation aa" + "\n")
                            for each_record in run_log_list:
                                log_txt.write(each_record + "\n")

                            log_txt.write("\nTip: 'Unknown or ?' indicates that the site maybe degenerate aa.")
                    with open((savefile_path + "_mutation_site_summary.txt"),
                      "w",
                      encoding="utf-8") as summary_txt:
                        summary_txt.write("Reference sequence: " + qury_seq_name + "\n\n" + "Site (aa)" + "\t" + "mutation" + "\t" + "changed in properties of aa?" + "\t" + "changed type of properties" + "\t" + "sequence count" + "\n")
                        for key in mutation_count_dic:
                            summary_txt.write(key + "\t" + str(mutation_count_dic[key]) + "\n")

                        summary_txt.write("\nTip: 'Unknown or ?' indicates that the site maybe degenerate bases.")
                    summary_2_file = savefile_path + "_mutation_site_(not_gap_termination_and_unknown_aa).txt"
                    with open(summary_2_file, "w", encoding="utf-8") as sus_site_no_others_count:
                        sus_site_no_others_count.write("Reference sequence: " + qury_seq_name + "\n\n" + "Site" + "\t" + "mutation" + "\t" + "changed in properties of aa?" + "\t" + "changed type of properties" + "\t" + "sequence count" + "\n")
                        for key in mutation_jx_count:
                            sus_site_no_others_count.write(key + "\t" + str(mutation_jx_count[key]) + "\n")

                    self.trigger.emit(95)
                    if draw_frequency_dis == "Yes":
                        site_fru = {}
                        for key in mutation_jx_count:
                            site = key.split("\t")[0]
                            if site not in site_fru:
                                site_fru[site] = mutation_jx_count[key]
                            else:
                                site_fru[site] += mutation_jx_count[key]

                        nuc_aa_plot(site_fru, frequency_group_tab, savefile_path)
                else:
                    if work_modes == "linked_modes":
                        self.trigger.emit(40)
                        max_num_seq = self.para_dic["max_num_seq"]
                        qury_seq = ""
                        qury_seq_name = ""
                        count_dic = {}
                        if muta_data_type == "nucleotide":
                            seq_lenth = 0
                            run_flage = True
                            run_nums = 0
                            run_log_list = []
                            while run_flage:
                                run_nums += 1
                                all_seq_num, seq_tab, run_flage = piecewise_read(input_file_path, max_num_seq, run_nums)
                                seq_count = len(seq_tab)
                                p2 = Pool(int(self.para_dic["thread_num"]))
                                seq_mu_list = []
                                for n in range(seq_count):
                                    seq_index = seq_tab[n][0]
                                    seq_name = seq_tab[n][1]
                                    seq_contain = seq_tab[n][2].replace("U", "T")
                                    if seq_index == 1:
                                        qury_seq_name = seq_name
                                        qury_seq = seq_contain
                                        seq_lenth = len(qury_seq)
                                    else:
                                        seq_mu = p2.apply_async(link_mutation_run1, args=(
                                         seq_name, seq_lenth, seq_contain, qury_seq))
                                        seq_mu_list.append(seq_mu)

                                p2.close()
                                p2.join()
                                jindu = run_nums / ceil(all_seq_num / max_num_seq)
                                self.trigger.emit(45 + int(40 * jindu))
                                for each_seq in seq_mu_list:
                                    each_seq_mutation = each_seq.get()
                                    if each_seq_mutation:
                                        mutation_label = each_seq_mutation[-1]
                                        if mutation_label not in count_dic:
                                            count_dic[mutation_label] = 1
                                        else:
                                            count_dic[mutation_label] += 1

                                if self.para_dic["save_log"]:
                                    for each_seq in seq_mu_list:
                                        each_seq_mutation = each_seq.get()
                                        if each_seq_mutation:
                                            run_log_list.append(each_seq_mutation[0] + "\t" + each_seq_mutation[1])

                            if run_log_list:
                                with open((savefile_path + "_run_logs.txt"), "w", encoding="utf-8") as run_logs:
                                    run_logs.write("Reference sequence: " + qury_seq_name + "\n\n" + "Seq_names" + "\t" + "Mutation nt" + "\n")
                                    for each_record in run_log_list:
                                        run_logs.write(each_record + "\n")

                        elif muta_data_type == "protein":
                            seq_lenth = 0
                            run_flage = True
                            run_nums = 0
                            run_log_list = []
                            while run_flage:
                                run_nums += 1
                                all_seq_num, seq_tab, run_flage = piecewise_read(input_file_path, max_num_seq, run_nums)
                                seq_count = len(seq_tab)
                                p2 = Pool(int(self.para_dic["thread_num"]))
                                seq_mu_list = []
                                for n in range(seq_count):
                                    seq_index = seq_tab[n][0]
                                    seq_name = seq_tab[n][1]
                                    seq_contain = seq_tab[n][2]
                                    if seq_index == 1:
                                        qury_seq_name = seq_name
                                        qury_seq = seq_contain
                                        seq_lenth = len(qury_seq)
                                    else:
                                        seq_mu = p2.apply_async(link_mutation_run1, args=(
                                         seq_name, seq_lenth, seq_contain, qury_seq))
                                        seq_mu_list.append(seq_mu)

                                p2.close()
                                p2.join()
                                jindu = run_nums / ceil(all_seq_num / max_num_seq)
                                self.trigger.emit(45 + int(40 * jindu))
                                for each_seq in seq_mu_list:
                                    each_seq_mutation = each_seq.get()
                                    if each_seq_mutation:
                                        mutation_label = each_seq_mutation[-1]
                                        if mutation_label not in count_dic:
                                            count_dic[mutation_label] = 1
                                        else:
                                            count_dic[mutation_label] += 1

                                if self.para_dic["save_log"]:
                                    for each_seq in seq_mu_list:
                                        each_seq_mutation = each_seq.get()
                                        if each_seq_mutation:
                                            run_log_list.append(each_seq_mutation[0] + "\t" + each_seq_mutation[1])

                            if run_log_list:
                                with open((savefile_path + "_run_logs.txt"), "w", encoding="utf-8") as run_logs:
                                    run_logs.write("Reference sequence: " + qury_seq_name + "\n\n" + "Seq_names" + "\t" + "Mutation aa" + "\n")
                                    for each_record in run_log_list:
                                        run_logs.write(each_record + "\n")

                        elif muta_data_type == "codon":
                            run_flage = True
                            run_nums = 0
                            codon_count = 0
                            run_log_list = []
                            while run_flage:
                                run_nums += 1
                                all_seq_num, seq_tab, run_flage = piecewise_read(input_file_path, max_num_seq, run_nums)
                                seq_count = len(seq_tab)
                                p2 = Pool(int(self.para_dic["thread_num"]))
                                seq_mu_list = []
                                for n in range(seq_count):
                                    seq_index = seq_tab[n][0]
                                    seq_name = seq_tab[n][1]
                                    seq_contain = seq_tab[n][2].replace("U", "T").strip()
                                    if seq_index == 1:
                                        qury_seq_name = seq_name
                                        qury_seq = seq_contain
                                        codon_count = int(len(qury_seq) / 3)
                                    else:
                                        seq_mu = p2.apply_async(link_mutation_run2, args=(
                                         seq_name, seq_contain, qury_seq, codon_count))
                                        seq_mu_list.append(seq_mu)

                                p2.close()
                                p2.join()
                                jindu = run_nums / ceil(all_seq_num / max_num_seq)
                                self.trigger.emit(45 + int(40 * jindu))
                                for each_seq in seq_mu_list:
                                    each_seq_mutation = each_seq.get()
                                    if each_seq_mutation:
                                        mutation_label = each_seq_mutation[-1]
                                        if mutation_label not in count_dic:
                                            count_dic[mutation_label] = 1
                                        else:
                                            count_dic[mutation_label] += 1

                                if self.para_dic["save_log"]:
                                    for each_seq in seq_mu_list:
                                        each_seq_mutation = each_seq.get()
                                        if each_seq_mutation:
                                            run_log_list.append(each_seq_mutation[0] + "\t" + each_seq_mutation[1])

                            if run_log_list:
                                with open((savefile_path + "_run_logs.txt"), "w", encoding="utf-8") as run_logs:
                                    run_logs.write("Reference sequence: " + qury_seq_name + "\n\n" + "Seq_names" + "\t" + "Mutation codon" + "\n")
                                    for each_record in run_log_list:
                                        run_logs.write(each_record + "\n")

                        self.trigger.emit(90)
                        with open((savefile_path + "_summary.txt"), "w",
                          encoding="utf-8") as summay_file:
                            summay_file.write("Reference sequence: " + qury_seq_name + "\n\n" + "Mutation" + "\t" + "Sequence count" + "\n")
                            for key in count_dic:
                                summay_file.write(key + "\t" + str(count_dic[key]) + "\n")

                            summay_file.write("\nTip: the string '-' indicates insertion or deletion.")
                self.trigger.emit(100)
                self.trigger_save.emit(saves_dir)
            except Exception as e:
                try:
                    config.run_error_log = "traceback.format_exc():\n%s" % traceback.format_exc()
                    self.trigger.emit(1314)
                finally:
                    e = None
                    del e

        finally:
            self.exit(0)
