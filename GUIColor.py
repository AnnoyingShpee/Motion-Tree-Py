from PySide6.QtCore import QSize
from PySide6.QtGui import QAction, QIcon, QPalette, QColor
from PySide6.QtWidgets import QMainWindow, QToolBar, QPushButton, QStatusBar, QMessageBox, QVBoxLayout, QHBoxLayout, QWidget


class Color(QWidget):
    def __init__(self, color):
        super(Color, self).__init__()
        self.setAutoFillBackground(True)

        palette = self.palette()
        palette.setColor(QPalette.Window, QColor(color))
        self.setPalette(palette)

