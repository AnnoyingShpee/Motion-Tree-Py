import sys
import ctypes
from BuilderPageGUI import BuilderPage, DynDomBuilderPage
from HelpPageGUI import HelpPage
from OutputWindowGUI import OutputWindow
from PySide6.QtCore import Qt, Slot, QThreadPool
from PySide6.QtGui import QKeyEvent
from PySide6.QtWidgets import QApplication, QMainWindow, QStatusBar, QMessageBox, \
    QVBoxLayout, QWidget, QTabWidget
from DataMngr import conn
user32 = ctypes.windll.user32
screensize = user32.GetSystemMetrics(0), user32.GetSystemMetrics(1)
# Used to fit display into correct resolution
user32.SetProcessDPIAware()

gui_app = QApplication(sys.argv)


class MainWindow(QMainWindow):
    def __init__(self, app):
        super().__init__()
        self.app = app  # Declare an app member
        self.setWindowTitle("Motion Tree Builder")
        self.resize(screensize[0] * 0.4, screensize[1] * 0.4)
        self.move(screensize[0] * 0.05, screensize[1] * 0.05)

        self.standard_builder_page = BuilderPage()
        # self.database_page = DatabasePage()
        self.dyndom_builder_page = DynDomBuilderPage()
        self.help_page = HelpPage()

        self.main_widget = QTabWidget()
        self.layout = QVBoxLayout()
        self.main_widget.addTab(self.standard_builder_page, "Standard")
        self.main_widget.addTab(self.dyndom_builder_page, "DynDom")
        self.main_widget.addTab(self.help_page, "Help")

        self.standard_builder_page.build_button.clicked.connect(self.build_button_on_click)
        self.dyndom_builder_page.build_button.clicked.connect(self.build_button_on_click)

        # self.layout.addLayout(image_output_widget.layout, 0, 1)
        # self._layout.addLayout()
        self.setCentralWidget(self.main_widget)

        self.windows = []
        self.max_window_count = 10

        self.thread_pool = QThreadPool()

    def keyPressEvent(self, event: QKeyEvent):
        if event.key() == Qt.Key.Key_Enter or event.key() == Qt.Key.Key_Return:
            self.build_button_on_click()
        elif event.key() == Qt.Key.Key_Escape:
            self.quit_app()
        super().keyPressEvent(event)

    def build_button_on_click(self):
        print(self.main_widget.currentIndex())
        empty_boxes = self.check_protein_text_boxes(self.main_widget.currentIndex())
        if len(empty_boxes) > 0:
            self.create_error_dialog("invalid_text_box", empty_boxes)
            return 1

        # if len(self.windows) >= self.max_window_count:
        #     self.statusBar().showMessage("Too many windows. Please remove some")
        #     self.create_error_dialog("too_many_windows")
        #     return 1
        if self.main_widget.currentIndex() == 0:
            input_path = self.standard_builder_page.file_form_widget.widgets["input_path_text_box"].text()
            output_path = self.standard_builder_page.file_form_widget.widgets["output_path_text_box"].text()
            protein_1_name = self.standard_builder_page.file_form_widget.widgets["protein_1_text_box"].text().lower()
            chain_1 = self.standard_builder_page.file_form_widget.widgets["protein_1_chain_text_box"].text().upper()
            protein_2_name = self.standard_builder_page.file_form_widget.widgets["protein_2_text_box"].text().lower()
            chain_2 = self.standard_builder_page.file_form_widget.widgets["protein_2_chain_text_box"].text().upper()
            spatial_proximity = float(self.standard_builder_page.param_form_widget.widgets["spatial_proximity"].value())
            small_node = int(self.standard_builder_page.param_form_widget.widgets["small_node"].value())
            clust_size = int(self.standard_builder_page.param_form_widget.widgets["clust_size"].value())
            magnitude = int(self.standard_builder_page.param_form_widget.widgets["magnitude"].value())

            print("Building Motion Tree Standard")
            new_window = OutputWindow(self.thread_pool, screensize, input_path, output_path,
                                      protein_1_name, chain_1, protein_2_name, chain_2,
                                      spatial_proximity, small_node, clust_size, magnitude)
            new_window.closed.connect(self.remove_window)
            new_window.show()
            self.windows.append(new_window)
        elif self.main_widget.currentIndex() == 1:
            input_path = self.dyndom_builder_page.file_form_widget.widgets["input_path_text_box"].text()
            output_path = self.dyndom_builder_page.file_form_widget.widgets["output_path_text_box"].text()
            protein_name = self.dyndom_builder_page.file_form_widget.widgets["protein_text_box"].text().upper()
            spatial_proximity = float(
                self.dyndom_builder_page.param_form_widget.widgets["spatial_proximity"].value())
            small_node = int(self.dyndom_builder_page.param_form_widget.widgets["small_node"].value())
            clust_size = int(self.dyndom_builder_page.param_form_widget.widgets["clust_size"].value())
            magnitude = int(self.dyndom_builder_page.param_form_widget.widgets["magnitude"].value())

            print("Building Motion Tree DynDom")
            new_window = OutputWindow(self.thread_pool, screensize, input_path, output_path,
                                      protein_name, None, None, None,
                                      spatial_proximity, small_node, clust_size, magnitude)
            new_window.closed.connect(self.remove_window)
            new_window.show()
            self.windows.append(new_window)

    def check_protein_text_boxes(self, page):
        empty_boxes = []
        if page == 0:
            if len(self.standard_builder_page.file_form_widget.widgets["input_path_text_box"].text().strip()) < 1:
                empty_boxes.append("Input Path")
            if len(self.standard_builder_page.file_form_widget.widgets["output_path_text_box"].text().strip()) < 1:
                empty_boxes.append("Output Path")
            if len(self.standard_builder_page.file_form_widget.widgets["protein_1_text_box"].text().strip()) < 4:
                empty_boxes.append("Protein 1")
            if len(self.standard_builder_page.file_form_widget.widgets["protein_2_text_box"].text().strip()) < 4:
                empty_boxes.append("Protein 2")
            if len(self.standard_builder_page.file_form_widget.widgets["protein_1_chain_text_box"].text().strip()) < 1 or \
                    self.standard_builder_page.file_form_widget.widgets["protein_1_chain_text_box"].text().isdigit():
                empty_boxes.append("Protein 1 Chain")
            if len(self.standard_builder_page.file_form_widget.widgets["protein_2_chain_text_box"].text().strip()) < 1 or \
                    self.standard_builder_page.file_form_widget.widgets["protein_2_chain_text_box"].text().isdigit():
                empty_boxes.append("Protein 2 Chain")
        elif page == 1:
            if len(self.dyndom_builder_page.file_form_widget.widgets["input_path_text_box"].text().strip()) < 1:
                empty_boxes.append("Input Path")
            if len(self.dyndom_builder_page.file_form_widget.widgets["output_path_text_box"].text().strip()) < 1:
                empty_boxes.append("Output Path")
            if len(self.dyndom_builder_page.file_form_widget.widgets["protein_text_box"].text().strip()) < 4:
                empty_boxes.append("Protein")
        return empty_boxes

    def create_error_dialog(self, *args):
        """
        Creates a QMessageBox that informs the user that text boxes have invalid inputs
        :return:
        """
        err_msg = ""
        if args[0] == "invalid_text_box":
            if len(args[1]) > 1:
                err_msg = "Multiple text boxes with invalid inputs!"
            else:
                err_msg = f"{args[1]} text is invalid!"
        elif args[0] == "too_many_windows":
            err_msg = f"Too many windows are opened. Please close a few. Maximum windows allowed to be opened: {self.max_window_count}"
        issue = QMessageBox.critical(self, "Error", err_msg, buttons=QMessageBox.Ok)
        if issue == QMessageBox.Ok:
            print("Ok clicked")

    @Slot(QWidget)
    def remove_window(self, window_to_close):
        self.windows.remove(window_to_close)

    def quit_app(self):
        self.app.quit()
        conn.close()
        sys.exit()

    def closeEvent(self, event) -> None:
        self.quit_app()


window = MainWindow(gui_app)
window.show()
gui_app.exec()

