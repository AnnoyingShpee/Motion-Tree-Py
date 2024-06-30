import sys
import ctypes
from pathlib import Path
from PySide6.QtCore import Qt, QSize, QRunnable, Slot, QThreadPool
from PySide6.QtGui import QAction, QIcon, QPalette, QColor
from PySide6.QtWidgets import QApplication, QMainWindow, QToolBar, QPushButton, QStatusBar, QFileDialog, QMessageBox, \
    QGridLayout, QVBoxLayout, QHBoxLayout, QStackedLayout, QWidget, QGridLayout, QLabel, QLineEdit, QSlider
from MotionTree import MotionTreeTest
from FileMngr import read_file_paths, read_param_file

user32 = ctypes.windll.user32
screensize = user32.GetSystemMetrics(0), user32.GetSystemMetrics(1)
# Used to fit display into correct resolution
user32.SetProcessDPIAware()

gui_app = QApplication(sys.argv)
with open('gui_styles.qss', 'r') as f:
    style = f.read()
    # Set the stylesheet of the application
    gui_app.setStyleSheet(style)


class MainWindow(QMainWindow):
    def __init__(self, app):
        super().__init__()
        self.app = app  # Declare an app member
        self.setWindowTitle("Motion Tree Builder")
        self.resize(screensize[0]*0.8, screensize[1]*0.8)

        self.file_form_widget = FileForm()
        self.param_form_widget = ParametersForm()
        self.image_output_widget = OutputForm()

        self.layout = QGridLayout()
        self.layout.addLayout(self.file_form_widget.layout, 0, 0)
        self.layout.addLayout(self.image_output_widget.layout, 0, 1)
        self.layout.addLayout(self.param_form_widget.layout, 1, 0)
        self.param_form_widget.widgets["build_button"].clicked.connect(self.build_button_on_click)
        self.main_widget = QWidget()
        self.main_widget.setLayout(self.layout)
        # self.layout.addLayout(image_output_widget.layout, 0, 1)
        # self._layout.addLayout()
        self.setCentralWidget(self.main_widget)

        self.files_params_dict = {
            "input_path": self.file_form_widget.widgets["input_path_text_box"].text(),
            "output_path": self.file_form_widget.widgets["output_path_text_box"].text(),
            "protein1": "",
            "chain1id": "",
            "protein2": "",
            "chain2id": "",
            "threshold": self.param_form_widget.widgets["spatial_proximity"].value(),
            "magnitude": self.param_form_widget.widgets["magnitude"].value(),
            "avg_diff": self.param_form_widget.widgets["avg_diff"].value()
        }

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

        # action_3 = QAction("Error", self)
        # action_3.setStatusTip("About to give an error")
        # action_3.triggered.connect(
        #     lambda error_msg="Error given": self.button_clicked_error(error_msg)
        # )
        # tool_bar.addAction(action_3)

        self.param_form_widget.widgets["spatial_proximity"].valueChanged.connect(
            lambda: self.param_slider_on_change("spatial_proximity")
        )
        self.param_form_widget.widgets["magnitude"].valueChanged.connect(
            lambda: self.param_slider_on_change("magnitude")
        )
        self.param_form_widget.widgets["avg_diff"].valueChanged.connect(
            lambda: self.param_slider_on_change("avg_diff")
        )

        # Status bar
        self.setStatusBar(QStatusBar(self))

    def toolbar_action_1_click(self):
        print("Action triggered")
        self.statusBar().showMessage("Action Triggered")

    def toolbar_action_2_click(self):
        print("Image clicked")
        self.statusBar().showMessage("Image Clicked")

    def param_slider_on_change(self, param_name):
        self.param_form_widget.widgets[f"{param_name}_value"].setText(str(self.param_form_widget.widgets[param_name].value()))
        self.files_params_dict[param_name] = self.param_form_widget.widgets[param_name].value()

    def build_button_on_click(self):
        empty_boxes = self.check_protein_text_boxes()
        empty_count = len(empty_boxes)
        if empty_count > 0:
            self.statusBar().showMessage("Invalid text")
            self.invalid_protein_strings_error(empty_count, empty_boxes)
        else:
            self.files_params_dict["protein1"] = self.file_form_widget.widgets["protein_1_text_box"].text()
            self.files_params_dict["protein2"] = self.file_form_widget.widgets["protein_2_text_box"].text()
            self.files_params_dict["chain1id"] = self.file_form_widget.widgets["protein_1_chain_text_box"].text()
            self.files_params_dict["chain2id"] = self.file_form_widget.widgets["protein_1_chain_text_box"].text()
            # print("Text:", self.file_form_widget.widgets["protein_1_text_box"].text())
            print("Building Motion Tree")
            self.statusBar().showMessage("Building Motion Tree")
            files_dict = read_file_paths()
            # Read parameter file to get parameters ( Protein PDB file names, protein chains, window size, domain size, ratio )
            param_dict = read_param_file()
            engine = MotionTreeTest(files_dict, param_dict)
            print(engine.protein_1.main_atoms_coords.shape)
            engine.run()

    def check_protein_text_boxes(self):
        empty_boxes = []
        if len(self.file_form_widget.widgets["input_path_text_box"].text().strip()) < 1:
            empty_boxes.append("Input")
        if len(self.file_form_widget.widgets["output_path_text_box"].text().strip()) < 1:
            empty_boxes.append("Output")
        if len(self.file_form_widget.widgets["protein_1_text_box"].text().strip()) < 4:
            empty_boxes.append("Protein 1")
        if len(self.file_form_widget.widgets["protein_2_text_box"].text().strip()) < 4:
            empty_boxes.append("Protein 2")
        if len(self.file_form_widget.widgets["protein_1_chain_text_box"].text().strip()) < 1 or self.file_form_widget.widgets["protein_1_chain_text_box"].text().isdigit():
            empty_boxes.append("Protein 1 Chain")
        if len(self.file_form_widget.widgets["protein_2_chain_text_box"].text().strip()) < 1 or self.file_form_widget.widgets["protein_2_chain_text_box"].text().isdigit():
            empty_boxes.append("Protein 2 Chain")
        return empty_boxes

    def invalid_protein_strings_error(self, empty_count, empty_boxes):
        if empty_count > 1:
            err_msg = "Multiple boxes with invalid inputs!"
        else:
            err_msg = f"{empty_boxes[0]} text is invalid!"
        issue = QMessageBox.critical(self, "Error", err_msg, buttons=QMessageBox.Ok)
        if issue == QMessageBox.Ok:
            print("Ok clicked")

    def quit_app(self):
        self.app.quit()


class FileForm(QWidget):
    def __init__(self):
        super().__init__()
        self.setObjectName("file-form")
        # Initialise the layout to place the widgets
        self.layout = QGridLayout()

        self._widgets_info = {
            "input_path_label": {
                "row": 0,
                "col": 0
            },
            "input_path_text_box": {
                "row": 1,
                "col": 0
            },
            "input_path_button": {
                "row": 1,
                "col": 1
            },
            "output_path_label": {
                "row": 2,
                "col": 0
            },
            "output_path_text_box": {
                "row": 3,
                "col": 0
            },
            "output_path_button": {
                "row": 3,
                "col": 1
            },
            "protein_1_label": {
                "row": 4,
                "col": 0
            },
            "protein_1_text_box": {
                "placeholder": "4ake",
                "row": 5,
                "col": 0,
                "max_length": 4
            },
            "protein_1_chain_label": {
                "row": 6,
                "col": 0
            },
            "protein_1_chain_text_box": {
                "placeholder": "A",
                "row": 7,
                "col": 0,
                "max_length": 1
            },
            "protein_2_label": {
                "row": 8,
                "col": 0
            },
            "protein_2_text_box": {
                "placeholder": "2eck",
                "row": 9,
                "col": 0,
                "max_length": 4
            },
            "protein_2_chain_label": {
                "row": 10,
                "col": 0
            },
            "protein_2_chain_text_box": {
                "placeholder": "B",
                "row": 11,
                "col": 0,
                "max_length": 1
            }
        }

        self.widgets = {
            "input_path_label": QLabel("PDB files directory. (The directory where the files are located)"),
            "input_path_text_box": QLineEdit("/data/input/pdb"),
            "input_path_button": QPushButton("Choose Folder"),
            "output_path_label": QLabel("Output directory. (The directory where the results are stored)"),
            "output_path_text_box": QLineEdit("/data/output"),
            "output_path_button": QPushButton("Choose Folder"),
            "protein_1_label": QLabel("Protein 1. Give the 4-character code of the protein."),
            "protein_1_text_box": QLineEdit(),
            "protein_1_chain_label": QLabel("The chain of Protein 1."),
            "protein_1_chain_text_box": QLineEdit(),
            "protein_2_label": QLabel("Protein 2. Give the 4-character code of the protein."),
            "protein_2_text_box": QLineEdit(),
            "protein_2_chain_label": QLabel("The chain of Protein 2."),
            "protein_2_chain_text_box": QLineEdit()
        }

        self.widgets["input_path_button"].clicked.connect(
            lambda: self.on_path_button_click("input_path_text_box")
        )

        self.widgets["output_path_button"].clicked.connect(
            lambda: self.on_path_button_click("output_path_text_box")
        )

        # self._widgets["protein_1_text_box"].setProperty("mandatoryField", True)
        self.widgets["protein_1_text_box"].setPlaceholderText(
            self._widgets_info["protein_1_text_box"]["placeholder"]
        )
        self.widgets["protein_1_text_box"].setMaxLength(
            self._widgets_info["protein_1_text_box"]["max_length"]
        )

        self.widgets["protein_1_chain_text_box"].setPlaceholderText(
            self._widgets_info["protein_1_chain_text_box"]["placeholder"]
        )
        self.widgets["protein_1_chain_text_box"].setMaxLength(
            self._widgets_info["protein_1_chain_text_box"]["max_length"]
        )

        # self._widgets["protein_2_text_box"].setProperty("mandatoryField", True)
        self.widgets["protein_2_text_box"].setPlaceholderText(
            self._widgets_info["protein_2_text_box"]["placeholder"]
        )
        self.widgets["protein_2_text_box"].setMaxLength(
            self._widgets_info["protein_2_text_box"]["max_length"]
        )

        self.widgets["protein_2_chain_text_box"].setPlaceholderText(
            self._widgets_info["protein_2_chain_text_box"]["placeholder"]
        )
        self.widgets["protein_2_chain_text_box"].setMaxLength(
            self._widgets_info["protein_2_chain_text_box"]["max_length"]
        )

        for w in self.widgets.keys():
            self.layout.addWidget(self.widgets[w], self._widgets_info[w]["row"], self._widgets_info[w]["col"])

    def on_path_button_click(self, p):
        dlg = str(QFileDialog.getExistingDirectory(self, "Select Folder"))
        cwd = str(Path.cwd())
        cwd = cwd.replace("\\", "/")
        tokens = dlg.split(cwd)
        for input_path in tokens:
            if len(input_path) > 0:
                self.widgets[p].setText(input_path)
                break


class ParametersForm(QWidget):
    def __init__(self):
        super().__init__()
        self.setObjectName("parameter-form")
        # Initialise the layout to place the widgets
        # self.layout = QGridLayout()
        self.layout = QVBoxLayout()

        # self._widgets_info = read_gui_json_file("gui_parameter_form")
        self._widgets_info = {
            "spatial_proximity": {
                "minimum": 4,
                "maximum": 8,
                "default": 7.0,
                "interval": 0.5
            },
            "magnitude": {
                "minimum": 5,
                "maximum": 10,
                "default": 5,
                "interval": 1
            },
            "avg_diff": {
                "minimum": 10,
                "maximum": 50,
                "default": 20,
                "interval": 5
            }
        }

        self.widgets = {
            "spatial_proximity_label": QLabel("The minimum spatial proximity between the closest residues: "),
            "spatial_proximity_min_label": QLabel(str(self._widgets_info["spatial_proximity"]["minimum"])),
            "spatial_proximity": DoubleSlider(Qt.Orientation.Horizontal),
            "spatial_proximity_max_label": QLabel(str(self._widgets_info["spatial_proximity"]["maximum"])),
            "spatial_proximity_value": QLabel(str(self._widgets_info["spatial_proximity"]["default"])),

            "magnitude_label": QLabel("The minimum magnitude for an effective node: "),
            "magnitude_min_label": QLabel(str(self._widgets_info["magnitude"]["minimum"])),
            "magnitude": DoubleSlider(Qt.Orientation.Horizontal),
            "magnitude_max_label": QLabel(str(self._widgets_info["magnitude"]["maximum"])),
            "magnitude_value": QLabel(str(self._widgets_info["magnitude"]["default"])),

            "avg_diff_label": QLabel("The number of residues to average when cluster is large enough: "),
            "avg_diff_min_label": QLabel(str(self._widgets_info["avg_diff"]["minimum"])),
            "avg_diff": DoubleSlider(Qt.Orientation.Horizontal),
            "avg_diff_max_label": QLabel(str(self._widgets_info["avg_diff"]["maximum"])),
            "avg_diff_value": QLabel(str(self._widgets_info["avg_diff"]["default"])),

            "build_button": QPushButton("Build Motion Tree")
        }

        self.widgets["spatial_proximity"].setMinimum(
            self._widgets_info["spatial_proximity"]["minimum"]
        )
        self.widgets["spatial_proximity"].setMaximum(
            self._widgets_info["spatial_proximity"]["maximum"]
        )
        self.widgets["spatial_proximity"].setInterval(
            self._widgets_info["spatial_proximity"]["interval"]
        )
        self.widgets["spatial_proximity"].setValue(
            self._widgets_info["spatial_proximity"]["default"]
        )

        self._spatial_proximity_label_layout = QHBoxLayout()
        self._spatial_proximity_label_layout.addWidget(self.widgets["spatial_proximity_label"])
        self._spatial_proximity_label_layout.addWidget(self.widgets["spatial_proximity_value"])
        self._spatial_proximity_label_widget = QWidget()
        self._spatial_proximity_label_widget.setLayout(self._spatial_proximity_label_layout)

        self._spatial_proximity_slider_layout = QHBoxLayout()
        self._spatial_proximity_slider_layout.addWidget(self.widgets["spatial_proximity_min_label"])
        self._spatial_proximity_slider_layout.addWidget(self.widgets["spatial_proximity"])
        self._spatial_proximity_slider_layout.addWidget(self.widgets["spatial_proximity_max_label"])
        self._spatial_proximity_slider_widget = QWidget()
        self._spatial_proximity_slider_widget.setLayout(self._spatial_proximity_slider_layout)

        self.widgets["magnitude"].setMinimum(
            self._widgets_info["magnitude"]["minimum"]
        )
        self.widgets["magnitude"].setMaximum(
            self._widgets_info["magnitude"]["maximum"]
        )
        self.widgets["magnitude"].setInterval(
            self._widgets_info["magnitude"]["interval"]
        )
        self.widgets["magnitude"].setValue(
            self._widgets_info["magnitude"]["default"]
        )

        self._magnitude_label_layout = QHBoxLayout()
        self._magnitude_label_layout.addWidget(self.widgets["magnitude_label"])
        self._magnitude_label_layout.addWidget(self.widgets["magnitude_value"])
        self._magnitude_label_widget = QWidget()
        self._magnitude_label_widget.setLayout(self._magnitude_label_layout)

        self._magnitude_slider_layout = QHBoxLayout()
        self._magnitude_slider_layout.addWidget(self.widgets["magnitude_min_label"])
        self._magnitude_slider_layout.addWidget(self.widgets["magnitude"])
        self._magnitude_slider_layout.addWidget(self.widgets["magnitude_max_label"])
        self._magnitude_slider_widget = QWidget()
        self._magnitude_slider_widget.setLayout(self._magnitude_slider_layout)

        self.widgets["avg_diff"].setMinimum(
            self._widgets_info["avg_diff"]["minimum"]
        )
        self.widgets["avg_diff"].setMaximum(
            self._widgets_info["avg_diff"]["maximum"]
        )
        self.widgets["avg_diff"].setInterval(
            self._widgets_info["avg_diff"]["interval"]
        )
        self.widgets["avg_diff"].setValue(
            self._widgets_info["avg_diff"]["default"]
        )

        self._avg_diff_label_layout = QHBoxLayout()
        self._avg_diff_label_layout.addWidget(self.widgets["avg_diff_label"])
        self._avg_diff_label_layout.addWidget(self.widgets["avg_diff_value"])
        self._avg_diff_label_widget = QWidget()
        self._avg_diff_label_widget.setLayout(self._avg_diff_label_layout)

        self._avg_diff_slider_layout = QHBoxLayout()
        self._avg_diff_slider_layout.addWidget(self.widgets["avg_diff_min_label"])
        self._avg_diff_slider_layout.addWidget(self.widgets["avg_diff"])
        self._avg_diff_slider_layout.addWidget(self.widgets["avg_diff_max_label"])
        self._avg_diff_slider_widget = QWidget()
        self._avg_diff_slider_widget.setLayout(self._avg_diff_slider_layout)

        self.layout.addWidget(self._spatial_proximity_label_widget)
        self.layout.addWidget(self._spatial_proximity_slider_widget)
        self.layout.addWidget(self._magnitude_label_widget)
        self.layout.addWidget(self._magnitude_slider_widget)
        self.layout.addWidget(self._avg_diff_label_widget)
        self.layout.addWidget(self._avg_diff_slider_widget)
        self.layout.addWidget(self.widgets["build_button"])


class OutputForm(QWidget):
    def __init__(self):
        super().__init__()
        # Initialise the layout to place the widgets
        self.layout = QVBoxLayout()

        self.widgets = {
            "dist_diff_button": QPushButton("Diff Dist Matrix"),
            "dist_diff_mat_img": QLabel("Diff Dist Matrix"),
            "motion_tree_button": QPushButton("Motion Tree"),
            "motion_tree_img": QLabel("Motion Tree"),
            "build_button": QPushButton("Build Motion Tree")
        }

        self._buttons_layout = QHBoxLayout()
        self._buttons_layout.addWidget(self.widgets["dist_diff_button"])
        self._buttons_layout.addWidget(self.widgets["motion_tree_button"])
        self._buttons_layout_widget = QWidget()
        self._buttons_layout_widget.setLayout(self._buttons_layout)

        self._image_output_layout = QStackedLayout()
        self._image_output_layout.addWidget(self.widgets["dist_diff_mat_img"])
        self._image_output_layout.addWidget(self.widgets["motion_tree_img"])
        self._image_output_layout_widget = QWidget()
        self._image_output_layout_widget.setLayout(self._image_output_layout)

        self.layout.addWidget(self._buttons_layout_widget)
        self.layout.addWidget(self._image_output_layout_widget)


class DoubleSlider(QSlider):
    """
    QSlider class to handle float values and intervals as well as custom step values
    https://stackoverflow.com/questions/42820380/use-float-for-qslider
    """
    def __init__(self, *args, **kargs):
        super(DoubleSlider, self).__init__( *args, **kargs)
        self._min = 0
        self._max = 99
        self.interval = 1

    def setValue(self, value):
        index = round((value - self._min) / self.interval)
        return super(DoubleSlider, self).setValue(index)

    def value(self):
        return self.index * self.interval + self._min

    @property
    def index(self):
        return super(DoubleSlider, self).value()

    def setIndex(self, index):
        return super(DoubleSlider, self).setValue(index)

    def setMinimum(self, value):
        self._min = value
        self._range_adjusted()

    def setMaximum(self, value):
        self._max = value
        self._range_adjusted()

    def setInterval(self, value):
        # To avoid division by zero
        if not value:
            raise ValueError('Interval of zero specified')
        self.interval = value
        self._range_adjusted()

    def _range_adjusted(self):
        number_of_steps = int((self._max - self._min) / self.interval)
        super(DoubleSlider, self).setMaximum(number_of_steps)


window = MainWindow(gui_app)
window.show()
gui_app.exec()

