# uncompyle6 version 3.9.3
# Python bytecode version base 3.7.0 (3394)
# Decompiled from: Python 3.12.10 (tags/v3.12.10:0cc8128, Apr  8 2025, 12:21:36) [MSC v.1943 64 bit (AMD64)]
# Embedded file name: convert_tools/file_merge.py
"""
Author: Zhou Zhi-Jian
Time: 2020/9/15 20:10

"""
import os, traceback, config, platform
from PySide2.QtWidgets import QApplication, QWidget, QPushButton, QPlainTextEdit, QHBoxLayout, QVBoxLayout, QMessageBox
from PySide2.QtGui import QIcon, QPalette, QCursor, QFont
from PySide2.QtCore import Qt, QThread, Signal
from custom_libraries.gui_class import QLineEdit_getDir, PopTipWindow, ErrorLogGUI, QProgressBar_beautify, Start_Button, MyGroupBox
from custom_libraries.public_functions import get_all_path

class FileMerge(QWidget):
    __doc__ = "\n    多个文件合并\n    "

    def __init__(self, parent=None):
        super(FileMerge, self).__init__(parent)
        self.setAcceptDrops(True)
        desktop = QApplication.desktop()
        des_w = desktop.width()
        des_h = desktop.height()
        palette = QPalette()
        palette.setColor(QPalette.Background, Qt.white)
        self.setPalette(palette)
        self.setWindowIcon(QIcon(":/merge.png"))
        self.setWindowTitle("File Merge")
        win_w = des_w * 0.35
        win_h = des_h * 0.5
        self.resize(int(win_w), int(win_h))
        self.initUI()

    def initUI(self):
        self.setToolTip("提示：<b>此功能将指定文件夹中的多个文本文件合并为一个</b>。注：输入文件夹中不要放与合并无关的其他任何文件。")
        self.file_edit = QLineEdit_getDir()
        self.file_edit.setFont(QFont(config.app_font))
        self.file_button = QPushButton("Preview file", self)
        self.file_button.setIcon(QIcon(":/eye.png"))
        self.file_button.setFont(QFont(config.app_font, config.app_font_size, QFont.Bold))
        self.file_button.setCursor(QCursor(Qt.PointingHandCursor))
        self.file_button.setToolTip("You can click button (or ignore it) to <b>preview</b> file from the folder")
        self.file_button.setMinimumSize(130, 40)
        self.file_button.setMaximumSize(160, 45)
        self.file_button.clicked.connect(self.load_file)
        file_box = MyGroupBox("&Input a folder", self)
        file_box.setToolTip("Drag a <b>folder which contains multiple files only</b> to be merged and drop here")
        file_hbox = QHBoxLayout()
        file_hbox.addWidget(self.file_edit)
        file_hbox.addWidget(self.file_button)
        file_box.setLayout(file_hbox)
        self.filepath = QPlainTextEdit(self)
        self.filepath.setFont(QFont(config.app_font, config.app_font_size))
        self.filepath.setReadOnly(True)
        self.filepath.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.filepath.setStyleSheet("border: 1px solid #8B8B7A;border-radius:2px")
        taxon_Vbox = QVBoxLayout()
        taxon_Vbox.addWidget(self.filepath)
        taxon_box = MyGroupBox("&File-paths preview", self)
        taxon_box.setLayout(taxon_Vbox)
        self.progressBar = QProgressBar_beautify()
        self.progressBar.setValue(0)
        self.make_button = Start_Button()
        self.make_button.clicked.connect(self.Work)
        progressBar_hbox = QHBoxLayout()
        progressBar_hbox.addWidget(self.progressBar)
        progressBar_hbox.addWidget(self.make_button)
        progressBar_box = MyGroupBox("&Run and Progress", self)
        progressBar_box.setLayout(progressBar_hbox)
        win_layout_V = QVBoxLayout()
        win_layout_V.setSpacing(20)
        win_layout_V.addWidget(file_box)
        win_layout_V.addWidget(taxon_box)
        win_layout_V.addWidget(progressBar_box)
        self.setLayout(win_layout_V)

    def Work(self):
        """
        运行前检查以及各类参数设置
        :return: 
        """
        file_dir = self.file_edit.text()
        if file_dir == "" or os.path.isfile(file_dir):
            return
        self.make_button.setEnabled(False)
        self.make_button.setText("Run...")
        self.thread_y = RunThread_y(file_dir)
        self.thread_y.trigger.connect(self.set_progressbar_value)
        self.thread_y.trigger_save.connect(self.run_finshed)
        self.thread_y.start()

    def load_file(self):
        """
        加载文件路径到文本框
        :return:
        """
        fir = self.file_edit.text()
        if fir == "" or fir.find(".") != -1:
            PopTipWindow.error("Error!", "Please drag into a folder first！")
            return
        a = get_all_path(self.file_edit.text())
        s = ""
        for line in a:
            s += line + "\n"

        self.filepath.setPlainText(s.strip())

    def set_progressbar_value(self, value):
        """
        进度条属性,如果收到“1314”信号，则弹出error log窗口
        :param value: 进度条数值
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

    def run_finshed(self, savefilepath):
        """
        运行结束时弹窗
        :param savefilepath: 储存文件所在路径
        :return:
        """
        PopTipWindow.success("Successfully produced!", "<b>Saved in: </b>" + savefilepath)

    def closeEvent(self, QCloseEvent):
        """
        重写关闭事件
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
        self.filepath.setPlainText("")
        if platform.system().lower() == "windows":
            try:
                self.thread_y.terminate()
            except:
                pass

            try:
                self.error.close()
            except:
                pass

        self.close()


class RunThread_y(QThread):
    __doc__ = "\n    运行线程\n    "
    trigger = Signal(int)
    trigger_save = Signal(str)

    def __init__(self, file_dir):
        """
        :param file_dir: 输入文件夹的路径

        """
        super(RunThread_y, self).__init__()
        self.file_dir = file_dir
        self.savefile_path = file_dir + "/merge.txt"

    def run(self):
        try:
            try:
                self.trigger.emit(10)
                file_list = get_all_path(self.file_dir)
                file_count = len(file_list)
                offset = int(file_count / 2)
                self.trigger.emit(20)
                with open(self.savefile_path, "w") as out:
                    for n in range(file_count):
                        with open(file_list[n], "r") as f1:
                            for line in f1:
                                out.write(line)

                        out.write("\n")
                        offset = offset + 1
                        proess = offset / file_count * 50 - 1
                        self.trigger.emit(int(proess))

                self.trigger.emit(100)
                self.trigger_save.emit(self.savefile_path)
            except Exception as e:
                try:
                    config.run_error_log = "traceback.format_exc():\n%s" % traceback.format_exc()
                    self.trigger.emit(1314)
                finally:
                    e = None
                    del e

        finally:
            self.exit(0)
