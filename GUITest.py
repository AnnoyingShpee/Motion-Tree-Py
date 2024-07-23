import sys
import ctypes
import traceback

from pathlib import Path
from PySide6.QtCore import Qt, QSize, QProcess, Signal, QObject, QRunnable, Slot, QThreadPool
from PySide6.QtGui import QAction, QIcon, QPixmap, QKeyEvent
from PySide6.QtWidgets import QApplication, QMainWindow, QToolBar, QPushButton, QStatusBar, QFileDialog, QMessageBox, \
    QGridLayout, QVBoxLayout, QHBoxLayout, QStackedLayout, QWidget, QGridLayout, QLabel, QLineEdit, QSlider, QComboBox
from MotionTree import MotionTree
from DataMngr import get_files, get_motion_tree_outputs

user32 = ctypes.windll.user32
screensize = user32.GetSystemMetrics(0), user32.GetSystemMetrics(1)
# Used to fit display into correct resolution
user32.SetProcessDPIAware()

gui_app = QApplication(sys.argv)
with open('gui_styles.qss', 'r') as f:
    style = f.read()
    # Set the stylesheet of the application
    gui_app.setStyleSheet(style)


class WorkerSignals(QObject):
    '''
    Defines the signals available from a running worker thread.

    Supported signals are:
    finished:
        No data
    error:
        tuple (exctype, value, traceback.format_exc() )
    result:
        object data returned from processing, anything
    progress:
        object data indicating progress
    '''
    finished = Signal()
    error = Signal(tuple)
    result = Signal(tuple)
    progress = Signal(object)


class Worker(QRunnable):
    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        # Add the callback to our kwargs
        self.kwargs['progress_callback'] = self.signals.progress

    @Slot()
    def run(self):
        try:
            result = self.fn(*self.args, **self.kwargs)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)  # Return the result of the processing
        finally:
            self.signals.finished.emit()  # Done


class MainWindow(QMainWindow):
    def __init__(self, app):
        super().__init__()
        self.app = app  # Declare an app member
        self.setWindowTitle("Motion Tree Builder")
        self.resize(screensize[0] * 0.4, screensize[1] * 0.4)
        self.move(screensize[0] * 0.1, screensize[1] * 0.1)

        self.file_form_widget = FileForm()
        self.param_form_widget = ParametersForm()

        self.input_layout = QVBoxLayout()
        self.input_layout.addLayout(self.file_form_widget.layout)
        self.input_layout.addLayout(self.param_form_widget.layout)

        self.param_form_widget.widgets["build_button"].clicked.connect(self.build_button_on_click)
        self.main_widget = QWidget()
        self.main_widget.setLayout(self.input_layout)
        # self.layout.addLayout(image_output_widget.layout, 0, 1)
        # self._layout.addLayout()
        self.setCentralWidget(self.main_widget)

        self.input_path = self.file_form_widget.widgets["input_path_text_box"].text()
        self.output_path = self.file_form_widget.widgets["output_path_text_box"].text()
        self.protein_1_name = ""
        self.chain_1 = ""
        self.protein_2_name = ""
        self.chain_2 = ""
        self.spatial_proximity = float(self.param_form_widget.widgets["spatial_proximity"].value())
        self.dissimilarity = int(self.param_form_widget.widgets["dissimilarity"].value())
        self.magnitude = int(self.param_form_widget.widgets["magnitude"].value())

        self.windows = []
        self.max_window_count = 10

        # Menubar and menus
        menu_bar = self.menuBar()
        # & sets the first letter as the ALT shortcut (ALT + f)
        builder_menu = menu_bar.addMenu("&Builder")
        settings_menu = menu_bar.addMenu("&Database")
        help_menu = menu_bar.addMenu("&Help")

        quit_menu = menu_bar.addMenu("Quit")
        quit_menu.triggered.connect(self.quit_app)

        self.param_form_widget.widgets["spatial_proximity"].valueChanged.connect(
            lambda: self.param_slider_on_change("spatial_proximity")
        )
        self.param_form_widget.widgets["magnitude"].valueChanged.connect(
            lambda: self.param_slider_on_change("magnitude")
        )
        self.param_form_widget.widgets["dissimilarity"].valueChanged.connect(
            lambda: self.param_slider_on_change("dissimilarity")
        )

        self.thread_pool = QThreadPool()

        # Status bar
        self.setStatusBar(QStatusBar(self))

    def param_slider_on_change(self, param_name):
        self.param_form_widget.widgets[f"{param_name}_value"].setText(
            str(self.param_form_widget.widgets[param_name].value()))

    def keyPressEvent(self, event: QKeyEvent):
        if event.key() == Qt.Key.Key_Enter or event.key() == Qt.Key.Key_Return:
            self.build_button_on_click()
        elif event.key() == Qt.Key.Key_Escape:
            self.quit_app()
        super().keyPressEvent(event)

    def build_button_on_click(self):
        empty_boxes = self.check_protein_text_boxes()
        if len(empty_boxes) > 0:
            self.statusBar().showMessage("Invalid text")
            self.create_error_dialog("invalid_text_box", empty_boxes)
            return 1

        if len(self.windows) >= self.max_window_count:
            self.statusBar().showMessage("Too many windows. Please remove some")
            self.create_error_dialog("too_many_windows")
            return 1

        self.input_path = self.file_form_widget.widgets["input_path_text_box"].text()
        self.output_path = self.file_form_widget.widgets["output_path_text_box"].text()
        self.protein_1_name = self.file_form_widget.widgets["protein_1_text_box"].text()
        self.chain_1 = self.file_form_widget.widgets["protein_1_chain_text_box"].text()
        self.protein_2_name = self.file_form_widget.widgets["protein_2_text_box"].text()
        self.chain_2 = self.file_form_widget.widgets["protein_2_chain_text_box"].text()
        self.spatial_proximity = float(self.param_form_widget.widgets["spatial_proximity"].value())
        self.dissimilarity = int(self.param_form_widget.widgets["dissimilarity"].value())
        self.magnitude = int(self.param_form_widget.widgets["magnitude"].value())

        print("Building Motion Tree")
        self.statusBar().showMessage("Building Motion Tree")
        new_window = OutputWindow(self.thread_pool, self.input_path, self.output_path,
                                  self.protein_1_name, self.chain_1, self.protein_2_name, self.chain_2,
                                  self.spatial_proximity, self.dissimilarity, self.magnitude)
        new_window.closed.connect(self.remove_window)
        new_window.show()
        self.windows.append(new_window)

        print("starting worker")
        # worker = Worker(self.build_motion_tree)
        # worker.signals.progress.connect(self.display_progress)
        # worker.signals.result.connect(self.display_output)
        # worker.signals.finished.connect(self.thread_complete)
        # worker.signals.error.connect(self.display_error)
        #
        # self.thread_pool.start(worker)

    def check_protein_text_boxes(self):
        empty_boxes = []
        if len(self.file_form_widget.widgets["input_path_text_box"].text().strip()) < 1:
            empty_boxes.append("Input Path")
        if len(self.file_form_widget.widgets["output_path_text_box"].text().strip()) < 1:
            empty_boxes.append("Output Path")
        if len(self.file_form_widget.widgets["protein_1_text_box"].text().strip()) < 4:
            empty_boxes.append("Protein 1")
        if len(self.file_form_widget.widgets["protein_2_text_box"].text().strip()) < 4:
            empty_boxes.append("Protein 2")
        if len(self.file_form_widget.widgets["protein_1_chain_text_box"].text().strip()) < 1 or \
                self.file_form_widget.widgets["protein_1_chain_text_box"].text().isdigit():
            empty_boxes.append("Protein 1 Chain")
        if len(self.file_form_widget.widgets["protein_2_chain_text_box"].text().strip()) < 1 or \
                self.file_form_widget.widgets["protein_2_chain_text_box"].text().isdigit():
            empty_boxes.append("Protein 2 Chain")
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
        sys.exit()


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
            "input_path_text_box": QLineEdit("data/input/pdb"),
            "input_path_button": QPushButton("Choose Folder"),
            "output_path_label": QLabel("Output directory. (The directory where the results are stored)"),
            "output_path_text_box": QLineEdit("data/output"),
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
        cwd = str(Path.cwd()) + "/"
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
        self.layout = QVBoxLayout()

        self._widgets_info = {
            "spatial_proximity": {
                "minimum": 4,
                "maximum": 8,
                "default": 7.0,
                "interval": 0.5
            },
            "dissimilarity": {
                "minimum": 10,
                "maximum": 50,
                "default": 20,
                "interval": 5
            },
            "magnitude": {
                "minimum": 5,
                "maximum": 10,
                "default": 5,
                "interval": 1
            }
        }

        self.widgets = {
            "spatial_proximity_label": QLabel("The minimum spatial proximity between the closest residues: "),
            "spatial_proximity_min_label": QLabel(str(self._widgets_info["spatial_proximity"]["minimum"])),
            "spatial_proximity": DoubleSlider(Qt.Orientation.Horizontal),
            "spatial_proximity_max_label": QLabel(str(self._widgets_info["spatial_proximity"]["maximum"])),
            "spatial_proximity_value": QLabel(str(self._widgets_info["spatial_proximity"]["default"])),

            "dissimilarity_label": QLabel("The number of residues to average when cluster is large enough: "),
            "dissimilarity_min_label": QLabel(str(self._widgets_info["dissimilarity"]["minimum"])),
            "dissimilarity": DoubleSlider(Qt.Orientation.Horizontal),
            "dissimilarity_max_label": QLabel(str(self._widgets_info["dissimilarity"]["maximum"])),
            "dissimilarity_value": QLabel(str(self._widgets_info["dissimilarity"]["default"])),

            "magnitude_label": QLabel("The minimum magnitude for an effective node: "),
            "magnitude_min_label": QLabel(str(self._widgets_info["magnitude"]["minimum"])),
            "magnitude": DoubleSlider(Qt.Orientation.Horizontal),
            "magnitude_max_label": QLabel(str(self._widgets_info["magnitude"]["maximum"])),
            "magnitude_value": QLabel(str(self._widgets_info["magnitude"]["default"])),

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

        self.widgets["dissimilarity"].setMinimum(
            self._widgets_info["dissimilarity"]["minimum"]
        )
        self.widgets["dissimilarity"].setMaximum(
            self._widgets_info["dissimilarity"]["maximum"]
        )
        self.widgets["dissimilarity"].setInterval(
            self._widgets_info["dissimilarity"]["interval"]
        )
        self.widgets["dissimilarity"].setValue(
            self._widgets_info["dissimilarity"]["default"]
        )

        self._dissimilarity_label_layout = QHBoxLayout()
        self._dissimilarity_label_layout.addWidget(self.widgets["dissimilarity_label"])
        self._dissimilarity_label_layout.addWidget(self.widgets["dissimilarity_value"])
        self._dissimilarity_label_widget = QWidget()
        self._dissimilarity_label_widget.setLayout(self._dissimilarity_label_layout)

        self._dissimilarity_slider_layout = QHBoxLayout()
        self._dissimilarity_slider_layout.addWidget(self.widgets["dissimilarity_min_label"])
        self._dissimilarity_slider_layout.addWidget(self.widgets["dissimilarity"])
        self._dissimilarity_slider_layout.addWidget(self.widgets["dissimilarity_max_label"])
        self._dissimilarity_slider_widget = QWidget()
        self._dissimilarity_slider_widget.setLayout(self._dissimilarity_slider_layout)

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

        self.layout.addWidget(self._spatial_proximity_label_widget)
        self.layout.addWidget(self._spatial_proximity_slider_widget)
        self.layout.addWidget(self._dissimilarity_label_widget)
        self.layout.addWidget(self._dissimilarity_slider_widget)
        self.layout.addWidget(self._magnitude_label_widget)
        self.layout.addWidget(self._magnitude_slider_widget)
        self.layout.addWidget(self.widgets["build_button"])


class OutputWindow(QWidget):
    closed = Signal(QWidget)

    def __init__(self, thread_pool, input_path, output_path, protein_1, chain_1, protein_2, chain_2, spat_prox, dis, magnitude):
        super().__init__()

        self.setWindowTitle(f"{protein_1}_{chain_1}_{protein_2}_{chain_2}")

        self.thread_pool = thread_pool

        # Initialise the layout to place the widgets
        self.layout = QVBoxLayout()
        self.input_path = input_path
        self.output_path = output_path
        self.protein_1 = protein_1
        self.chain_1 = chain_1
        self.protein_2 = protein_2
        self.chain_2 = chain_2
        self.spat_prox = spat_prox
        self.diss = dis
        self.magnitude = magnitude

        self.widgets = {
            "diff_dist_button": QPushButton("Difference Dist Matrix"),
            "motion_tree_button": QPushButton("Motion Tree"),
            "diff_dist_mat_img": QLabel(""),
            "motion_tree_img": QLabel("")
        }

        self.widgets["diff_dist_button"].clicked.connect(
            lambda: self.change_image(0)
        )
        self.widgets["motion_tree_button"].clicked.connect(
            lambda: self.change_image(1)
        )

        combi_layout = QVBoxLayout()
        combi_widget = QWidget()
        protein_1_label = QLabel(f"Protein 1 : {self.protein_1} ({self.chain_1})")
        protein_2_label = QLabel(f"Protein 2 : {self.protein_2} ({self.chain_2})")
        sp_label = QLabel(f"Spatial Proximity : {self.spat_prox}")
        dis_label = QLabel(f"Dissimilarity : {self.diss}")
        mag_label = QLabel(f"Magnitude : {magnitude}")
        self.progress_label = QLabel(f"Progress : Running")

        combi_layout.addWidget(protein_1_label)
        combi_layout.addWidget(protein_2_label)
        combi_layout.addWidget(sp_label)
        combi_layout.addWidget(dis_label)
        combi_layout.addWidget(mag_label)
        combi_layout.addWidget(self.progress_label)
        combi_widget.setLayout(combi_layout)

        self._buttons_layout = QHBoxLayout()
        self._buttons_layout.addWidget(self.widgets["diff_dist_button"])
        self._buttons_layout.addWidget(self.widgets["motion_tree_button"])
        self._buttons_widget = QWidget()
        self._buttons_widget.setLayout(self._buttons_layout)

        self.image_output_layout = QStackedLayout()
        self.image_output_layout.addWidget(self.widgets["diff_dist_mat_img"])
        self.image_output_layout.addWidget(self.widgets["motion_tree_img"])
        self.image_output_widget = QWidget()
        self.image_output_widget.setLayout(self.image_output_layout)

        self.layout.addWidget(combi_widget)
        self.layout.addWidget(self._buttons_widget)
        self.layout.addWidget(self.image_output_widget)
        self.setLayout(self.layout)

        worker = Worker(self.build_motion_tree)
        worker.signals.progress.connect(self.display_progress)
        worker.signals.result.connect(self.display_output)
        worker.signals.finished.connect(self.thread_complete)
        worker.signals.error.connect(self.display_error)

        self.thread_pool.start(worker)

    def build_motion_tree(self, progress_callback):
        print("Motion Tree Test Init")
        engine = MotionTree(self.input_path, self.output_path, self.protein_1, self.chain_1, self.protein_2, self.chain_2,
                            self.spat_prox, self.diss, self.magnitude)
        print("Done Tree Class Init")
        progress_callback.emit(f"Initialising {self.protein_1}")
        progress_callback.emit(engine.init_protein(1))
        progress_callback.emit(f"Initialising {self.protein_2}")
        progress_callback.emit(engine.init_protein(2))
        engine.check_sequence_identity()
        progress_callback.emit("Creating Distance Difference Matrix")
        progress_callback.emit(engine.create_distance_difference_matrix())
        progress_callback.emit("Building Motion Tree")
        total_time, nodes, protein_str, param_str = engine.run()
        print("Get time taken")
        progress_callback.emit(
            f"Motion Tree {protein_str} built. \nTime taken: {total_time}s. \nNumber of Effective Nodes: {nodes}")
        print("Checking")
        diff_dist_npy, diff_dist_img, motion_tree_img = get_motion_tree_outputs(
            self.output_path, self.protein_1, self.chain_1, self.protein_2, self.chain_2, self.spat_prox, self.diss, self.magnitude)
        print("Result Callback")
        return diff_dist_npy, diff_dist_img, motion_tree_img, f"{protein_str}_{param_str}"

    def display_progress(self, msg):
        """
        Displays the progress made by the program when the Motion Tree is being built
        :param msg: The message to be displayed in the Status Bar
        :return:
        """
        print(msg)
        self.progress_label.setText(msg)

    def display_output(self, outputs):
        print("paths", outputs)
        diff_dist_img = QPixmap(outputs[1])
        self.widgets["diff_dist_mat_img"].setPixmap(diff_dist_img)
        motion_tree_img = QPixmap(outputs[2])
        self.widgets["motion_tree_img"].setPixmap(motion_tree_img)
        # self.statusBar().showMessage(f"Time taken {str(msg)}")

    def thread_complete(self):
        print("Thread complete")

    def display_error(self, msg):
        print(msg)
        self.progress_label.setText(f"Error : {str(msg[1])}")

    def change_image(self, index):
        self.image_output_layout.setCurrentIndex(index)

    def closeEvent(self, event):
        self.closed.emit(self)
        super().closeEvent(event)

class DoubleSlider(QSlider):
    """
    Custom QSlider class to handle float values and intervals as well as custom step values
    https://stackoverflow.com/questions/42820380/use-float-for-qslider
    """

    def __init__(self, *args, **kargs):
        super(DoubleSlider, self).__init__(*args, **kargs)
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

