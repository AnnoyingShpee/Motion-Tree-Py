from pathlib import Path
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QPushButton,QFileDialog, QVBoxLayout, QHBoxLayout, QWidget, QGridLayout, QLabel, \
    QLineEdit, QSlider


class BuilderPage(QWidget):
    def __init__(self):
        super().__init__()
        self.file_form_widget = FileFormWidget()
        self.param_form_widget = ParametersFormWidget()
        self.build_button = QPushButton("Build Motion Tree")

        self.layout = QVBoxLayout()
        self.layout.addLayout(self.file_form_widget.layout)
        self.layout.addLayout(self.param_form_widget.layout)
        self.layout.addWidget(self.build_button)

        self.setLayout(self.layout)


class FileFormWidget(QWidget):
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


class ParametersFormWidget(QWidget):
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
            "small_node": {
                "minimum": 0,
                "maximum": 10,
                "default": 5,
                "interval": 1
            },
            "clust_size": {
                "minimum": 10,
                "maximum": 50,
                "default": 30,
                "interval": 1
            },
            "magnitude": {
                "minimum": 1,
                "maximum": 30,
                "default": 5,
                "interval": 1
            }
        }

        self.widgets = {
            "spatial_proximity_label": QLabel("The minimum spatial proximity between the closest residues (Default 7): "),
            "spatial_proximity_min_label": QLabel(str(self._widgets_info["spatial_proximity"]["minimum"])),
            "spatial_proximity": DoubleSlider(Qt.Orientation.Horizontal),
            "spatial_proximity_max_label": QLabel(str(self._widgets_info["spatial_proximity"]["maximum"])),
            "spatial_proximity_value": QLabel(str(self._widgets_info["spatial_proximity"]["default"])),

            "small_node_label": QLabel("Number of residues in the smaller domain must be above value for a node to be effective (Default 5): "),
            "small_node_min_label": QLabel(str(self._widgets_info["small_node"]["minimum"])),
            "small_node": DoubleSlider(Qt.Orientation.Horizontal),
            "small_node_max_label": QLabel(str(self._widgets_info["small_node"]["maximum"])),
            "small_node_value": QLabel(str(self._widgets_info["small_node"]["default"])),

            "clust_size_label": QLabel("The minimum number of residues in a node to be considered effective (Default 30): "),
            "clust_size_min_label": QLabel(str(self._widgets_info["clust_size"]["minimum"])),
            "clust_size": DoubleSlider(Qt.Orientation.Horizontal),
            "clust_size_max_label": QLabel(str(self._widgets_info["clust_size"]["maximum"])),
            "clust_size_value": QLabel(str(self._widgets_info["clust_size"]["default"])),

            "magnitude_label": QLabel("The minimum magnitude for an effective node (Default: 5)"),
            "magnitude_min_label": QLabel(str(self._widgets_info["magnitude"]["minimum"])),
            "magnitude": DoubleSlider(Qt.Orientation.Horizontal),
            "magnitude_max_label": QLabel(str(self._widgets_info["magnitude"]["maximum"])),
            "magnitude_value": QLabel(str(self._widgets_info["magnitude"]["default"]))
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

        self.widgets["small_node"].setMinimum(
            self._widgets_info["small_node"]["minimum"]
        )
        self.widgets["small_node"].setMaximum(
            self._widgets_info["small_node"]["maximum"]
        )
        self.widgets["small_node"].setInterval(
            self._widgets_info["small_node"]["interval"]
        )
        self.widgets["small_node"].setValue(
            self._widgets_info["small_node"]["default"]
        )

        self._small_node_label_layout = QHBoxLayout()
        self._small_node_label_layout.addWidget(self.widgets["small_node_label"])
        self._small_node_label_layout.addWidget(self.widgets["small_node_value"])
        self._small_node_label_widget = QWidget()
        self._small_node_label_widget.setLayout(self._small_node_label_layout)

        self._small_node_slider_layout = QHBoxLayout()
        self._small_node_slider_layout.addWidget(self.widgets["small_node_min_label"])
        self._small_node_slider_layout.addWidget(self.widgets["small_node"])
        self._small_node_slider_layout.addWidget(self.widgets["small_node_max_label"])
        self._small_node_slider_widget = QWidget()
        self._small_node_slider_widget.setLayout(self._small_node_slider_layout)

        self.widgets["clust_size"].setMinimum(
            self._widgets_info["clust_size"]["minimum"]
        )
        self.widgets["clust_size"].setMaximum(
            self._widgets_info["clust_size"]["maximum"]
        )
        self.widgets["clust_size"].setInterval(
            self._widgets_info["clust_size"]["interval"]
        )
        self.widgets["clust_size"].setValue(
            self._widgets_info["clust_size"]["default"]
        )

        self._clust_size_label_layout = QHBoxLayout()
        self._clust_size_label_layout.addWidget(self.widgets["clust_size_label"])
        self._clust_size_label_layout.addWidget(self.widgets["clust_size_value"])
        self._clust_size_label_widget = QWidget()
        self._clust_size_label_widget.setLayout(self._clust_size_label_layout)

        self._clust_size_slider_layout = QHBoxLayout()
        self._clust_size_slider_layout.addWidget(self.widgets["clust_size_min_label"])
        self._clust_size_slider_layout.addWidget(self.widgets["clust_size"])
        self._clust_size_slider_layout.addWidget(self.widgets["clust_size_max_label"])
        self._clust_size_slider_widget = QWidget()
        self._clust_size_slider_widget.setLayout(self._clust_size_slider_layout)

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
        self.layout.addWidget(self._small_node_label_widget)
        self.layout.addWidget(self._small_node_slider_widget)
        self.layout.addWidget(self._clust_size_label_widget)
        self.layout.addWidget(self._clust_size_slider_widget)
        self.layout.addWidget(self._magnitude_label_widget)
        self.layout.addWidget(self._magnitude_slider_widget)


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