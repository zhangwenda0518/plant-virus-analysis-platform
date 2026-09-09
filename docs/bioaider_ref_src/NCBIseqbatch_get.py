# uncompyle6 version 3.9.3
# Python bytecode version base 3.7.0 (3394)
# Decompiled from: Python 3.12.10 (tags/v3.12.10:0cc8128, Apr  8 2025, 12:21:36) [MSC v.1943 64 bit (AMD64)]
# Embedded file name: download_tools/NCBIseqbatch_get.py
"""
Author: Zhou Zhi-Jian
Time: 2020/9/16 15:26

"""
import webbrowser as web, traceback, config, platform
from PySide2.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QComboBox, QPlainTextEdit
from PySide2.QtGui import QIcon, QPalette
from PySide2.QtCore import Qt, QThread, Signal
from custom_libraries.gui_class import ErrorLogGUI, QProgressBar_beautify, TextEdit_contain, Start_Button, MyGroupBox

class SeqBatchDownload(QWidget):
    __doc__ = "\n    NCBI序列批量下载器\n    "

    def __init__(self, parent=None):
        super(SeqBatchDownload, self).__init__(parent)
        self.setAcceptDrops(True)
        desktop = QApplication.desktop()
        des_w = desktop.width()
        des_h = desktop.height()
        palette = QPalette()
        palette.setColor(QPalette.Background, Qt.white)
        self.setPalette(palette)
        self.setWindowIcon(QIcon(":/sailing-ship.png"))
        self.setWindowTitle("Sequence Batch Download(NCBI)")
        win_w = des_w * 0.3
        win_h = des_h * 0.6
        self.resize(int(win_w), int(win_h))
        self.initUI()

    def initUI(self):
        self.taxon = TextEdit_contain()
        self.taxon.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.taxon.setStyleSheet("border: 1px solid #8B8B7A;border-radius:2px")
        self.taxon.setPlainText("AY567487\nAY994055\nCQ870486\nCS124012")
        taxon_hbox = QHBoxLayout()
        taxon_hbox.addWidget(self.taxon)
        taxon_box = MyGroupBox("&STEP 1 - Enter a list of sequence accession", self)
        taxon_box.setLayout(taxon_hbox)
        taxon_box.setToolTip("Please pasta one or multiple accession of sequence(s) and separated by newlines. <b>Note: No more than 180 accessions!</b>")
        self.selcet = QComboBox(self)
        self.selcet.addItem("Nucleotide")
        self.selcet.addItem("Protein")
        self.selcet.addItem("Gene")
        self.selcet.addItem("Genome")
        self.selcet.addItem("All databases")
        data_bases = MyGroupBox("&STEP 2 - Select a database", self)
        data_bases_hbox = QHBoxLayout()
        data_bases_hbox.addWidget(self.selcet)
        data_bases.setLayout(data_bases_hbox)
        self.progressBar = QProgressBar_beautify()
        self.strat_button = Start_Button()
        self.strat_button.clicked.connect(self.Work)
        progressBar_hbox = QHBoxLayout()
        progressBar_hbox.addWidget(self.progressBar)
        progressBar_hbox.addWidget(self.strat_button)
        progressBar_box = MyGroupBox("&STEP 3 - Run and Progress", self)
        progressBar_box.setLayout(progressBar_hbox)
        wins_V_layout = QVBoxLayout()
        wins_V_layout.setSpacing(30)
        wins_V_layout.addWidget(taxon_box)
        wins_V_layout.addWidget(data_bases)
        wins_V_layout.addWidget(progressBar_box)
        self.setLayout(wins_V_layout)

    def closeEvent(self, QCloseEvent):
        """
        重写窗口关闭事件
        :param QCloseEvent:
        :return:
        """
        self.progressBar.setValue(0)
        self.strat_button.setEnabled(True)
        self.strat_button.setText("Start")
        if platform.system().lower() == "windows":
            try:
                self.run_thread.terminate()
            except:
                pass

            try:
                self.error.close()
            except:
                pass

        self.close()

    def Work(self):
        """
        运行前检查和参数配置
        :return:
        """
        self.progressBar.setValue(0)
        seq_ID = self.taxon.toPlainText().strip()
        if seq_ID == "":
            return
        database_select = self.selcet.currentText()
        self.strat_button.setEnabled(False)
        self.strat_button.setText("Run...")
        self.run_thread = RunThread_s2(seq_ID, database_select)
        self.run_thread.single_progress.connect(self.set_progressbar_value)
        self.run_thread.start()

    def set_progressbar_value(self, value):
        """
        进度条,如果收到“1314”信号，则弹出error log窗口
        :param value: 进度条值
        :return:
        """
        self.progressBar.setValue(value)
        if value == 100:
            self.strat_button.setEnabled(True)
            self.strat_button.setText("Start")
            return
        if value == 1314:
            self.error = ErrorLogGUI()
            self.error.show()
            self.strat_button.setEnabled(True)
            self.strat_button.setText("Start")
            return


class RunThread_s2(QThread):
    __doc__ = "\n    寻找重复序列运行线程\n    "
    single_progress = Signal(int)

    def __init__(self, seq_taxon, database_select):
        """
        :param seq_taxon: 序列ID
        :param database_select: 选择的数据库

        """
        super(RunThread_s2, self).__init__()
        self.seq_taxon = seq_taxon
        self.database_select = database_select

    def run(self):
        try:
            try:
                self.single_progress.emit(20)
                seqID_list = self.seq_taxon.strip().split("\n")
                search_head = ""
                if self.database_select == "Nucleotide":
                    search_head = "https://www.ncbi.nlm.nih.gov/nuccore/?term="
                else:
                    if self.database_select == "Protein":
                        search_head = "https://www.ncbi.nlm.nih.gov/protein/?term="
                    else:
                        if self.database_select == "Gene":
                            search_head = "https://www.ncbi.nlm.nih.gov/gene/?term="
                        else:
                            if self.database_select == "Genome":
                                search_head = "https://www.ncbi.nlm.nih.gov/genome/?term="
                            else:
                                if self.database_select == "All databases":
                                    search_head = "https://www.ncbi.nlm.nih.gov/search/all/?term="
                search_link = search_head + "|".join(seqID_list)
                web.open_new_tab(search_link)
                self.single_progress.emit(100)
            except Exception as e:
                try:
                    config.run_error_log = "traceback.format_exc():\n%s" % traceback.format_exc()
                    self.single_progress.emit(1314)
                finally:
                    e = None
                    del e

        finally:
            self.exit(0)
