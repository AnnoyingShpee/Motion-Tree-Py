from PySide6.QtCore import QSize
from PySide6.QtGui import QAction, QIcon
from PySide6.QtWidgets import QMainWindow, QToolBar, QPushButton, QStatusBar, QMessageBox, QVBoxLayout, QHBoxLayout, QWidget


class MainWindow(QMainWindow):
    def __init__(self, app):
        super().__init__()
        self.app = app  # Declare an app member
        self.setWindowTitle("Motion Tree Builder")

        # Menubar and menus
        menu_bar = self.menuBar()
        # & sets the first letter as the ALT shortcut (ALT + f)
        file_menu = menu_bar.addMenu("&File")
        quit_action = file_menu.addAction("Quit")
        quit_action.triggered.connect(self.quit_app)

        edit_menu = menu_bar.addMenu("&Edit")
        copy_action = edit_menu.addAction("Copy")
        edit_menu.addAction("Cut")
        edit_menu.addAction("Paste")
        edit_menu.addAction("Undo")
        edit_menu.addAction("Redo")

        window_menu = menu_bar.addMenu("&Window")
        settings_menu = menu_bar.addMenu("&Settings")
        help_menu = menu_bar.addMenu("&Help")

        # The toolbar
        tool_bar = QToolBar("Main Toolbar")
        tool_bar.setIconSize(QSize(16, 16))
        self.addToolBar(tool_bar)

        # Add quit action to toolbar
        tool_bar.addAction(quit_action)

        # Add custom action to toolbar
        action_1 = QAction("Some action", self)
        action_1.setStatusTip("Status message for some action.")
        action_1.triggered.connect(self.toolbar_action_1_click)
        tool_bar.addAction(action_1)

        # Add custom action containing image to toolbar
        action_2 = QAction(QIcon("Dendrogram.png"), "Some other action", self)
        action_2.setStatusTip("Action with image")
        action_2.triggered.connect(self.toolbar_action_2_click)
        action_2.setCheckable(True)
        tool_bar.addAction(action_2)

        action_3 = QAction("Error", self)
        action_3.setStatusTip("About to give an error")
        action_3.triggered.connect(
            lambda error_msg="Error given": self.button_clicked_error(error_msg)
        )
        tool_bar.addAction(action_3)

        # Status bar
        self.setStatusBar(QStatusBar(self))


    def toolbar_action_1_click(self):
        print("Action triggered")
        self.statusBar().showMessage("Action Triggered")

    def toolbar_action_2_click(self):
        print("Image clicked")
        self.statusBar().showMessage("Image Clicked")

    # def button_clicked_error(self, error_msg):
    #     issue = QMessageBox.critical(QWidget, "Error", error_msg, QMessageBox.Ok | QMessageBox.Cancel)
    #     if issue == QMessageBox.Ok:
    #         print("Ok clicked")
    #     elif issue == QMessageBox.Cancel:
    #         print("Cancel clicked")


    def quit_app(self):
        self.app.quit()


class MainWidget(QWidget):
    def __init__(self):
        super().__init__()

        # Buttons
        button = QPushButton("Click")
        button.clicked.connect(self.button_clicked)
        layout = QVBoxLayout()
        layout.addWidget(button)
        self.setLayout(layout)

    def button_clicked(self):
        print("Done")
