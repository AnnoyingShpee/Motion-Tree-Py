from PySide6.QtWidgets import QMainWindow, QToolBar, QPushButton, QStatusBar, QMessageBox, QVBoxLayout, QHBoxLayout, \
    QWidget, QLabel, QLineEdit, QFileDialog, QGridLayout


class MainWidget(QWidget):
    def __init__(self):
        super().__init__()
        self.widgets_to_add = []
        self.vertical_layout = QVBoxLayout()
        self.horizontal_layout = QHBoxLayout()

        self.input_path_label = QLabel("PDB files directory. (The directory where the files are located)")
        self.input_path_text_box = QFileDialog()
        # self.input_path_text_box = QLineEdit()
        self.output_path_label = QLabel("Output directory. (The directory where the results are stored)")
        self.output_path_text_box = QLineEdit()
        self.widgets_to_add.extend([self.input_path_label, self.input_path_text_box, self.output_path_label, self.output_path_text_box])

        self.protein_1_label = QLabel("Protein 1. Give the 4-character code of the protein. ")
        self.protein_1_text_box = QLineEdit()
        self.protein_1_text_box.setMaxLength(4)
        self.protein_1_text_box.setPlaceholderText("4ake")
        self.protein_1_chain_label = QLabel("The chain of Protein 1.")
        self.protein_1_chain_text_box = QLineEdit()
        self.protein_1_chain_text_box.setMaxLength(1)
        self.protein_1_chain_text_box.setPlaceholderText("A")
        self.widgets_to_add.extend(
            [self.protein_1_label, self.protein_1_text_box,
             self.protein_1_chain_text_box, self.protein_1_chain_text_box]
        )

        self.protein_2_label = QLabel("Protein 2. Give the 4-character code of the protein. ")
        self.protein_2_text_box = QLineEdit()
        self.protein_2_text_box.setMaxLength(4)
        self.protein_2_text_box.setPlaceholderText("2eck")
        self.protein_2_chain_label = QLabel("The chain of Protein 2.")
        self.protein_2_chain_text_box = QLineEdit()
        self.protein_2_chain_text_box.setMaxLength(1)
        self.protein_2_chain_text_box.setPlaceholderText("A")
        self.widgets_to_add.extend(
            [self.protein_2_label, self.protein_2_text_box,
             self.protein_2_chain_text_box, self.protein_2_chain_text_box]
        )

        self.submit_button = QPushButton("Run")

        for widget in self.widgets_to_add:
            self.vertical_layout.addWidget(widget)

        self.setLayout(self.vertical_layout)

    def submit_error(self, error_msg):
        issue = QMessageBox.critical(self, "Error", error_msg, QMessageBox.Ok | QMessageBox.Cancel)
        if issue == QMessageBox.Ok:
            print("Ok clicked")
        elif issue == QMessageBox.Cancel:
            print("Cancel clicked")

    def add_widget_with_label(self, layout, widget, label_text):
        hbox = QHBoxLayout()
        label = QLabel(label_text)
        hbox.addWidget(label)
        hbox.addWidget(widget)
        layout.addLayout(hbox)

