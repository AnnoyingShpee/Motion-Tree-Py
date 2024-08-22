import sys
import ctypes
import traceback
from BuilderPageGUI import BuilderPage
from OutputWindowGUI import OutputWindow
from pathlib import Path
from PySide6.QtCore import Qt, QSize, QProcess, Signal, QObject, QRunnable, Slot, QThreadPool
from PySide6.QtGui import QAction, QIcon, QPixmap, QKeyEvent
from PySide6.QtWidgets import QApplication, QMainWindow, QToolBar, QPushButton, QStatusBar, QFileDialog, QMessageBox, \
    QVBoxLayout, QHBoxLayout, QStackedLayout, QWidget, QStackedWidget, QTabWidget, QGridLayout, QLabel, QLineEdit, QSlider
from MotionTree import MotionTree


class HelpPage(QWidget):
    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout()

