import sys
import numpy as np
import traceback
from PySide6.QtCore import Signal, QObject, QRunnable, Slot
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import QPushButton, QVBoxLayout, QHBoxLayout, QStackedLayout, QWidget, QLabel
from MotionTree import MotionTree
from DataMngr import conn, get_motion_tree_outputs, check_motion_tree_exists, get_motion_tree, \
    save_results_to_disk, write_info_file, write_to_pdb, write_domains_to_pml, insert_motion_tree, check_nodes_exist, \
    get_nodes, insert_nodes


class OutputWindow(QWidget):
    closed = Signal(QWidget)

    def __init__(self, thread_pool, screensize, input_path, output_path, protein_1, chain_1, protein_2, chain_2, spat_prox,
                 small_node, clust_size, magnitude):
        super().__init__()

        self.setWindowTitle(f"{protein_1}_{chain_1}_{protein_2}_{chain_2}")
        self.move(screensize[0] * 0.05, screensize[1] * 0.05)

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
        self.small_node = small_node
        self.clust_size = clust_size
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
        dis_label = QLabel(f"Cluster Size : {self.clust_size}")
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
        has_motion_tree, has_nodes = -1, -1
        if conn is not None:
            has_motion_tree = check_motion_tree_exists(self.protein_1, self.chain_1, self.protein_2, self.chain_2, self.spat_prox)
            has_nodes = check_nodes_exist(self.protein_1, self.chain_1, self.protein_2, self.chain_2,
                                          self.spat_prox, self.small_node, self.clust_size, self.magnitude)
        engine = MotionTree(self.input_path, self.output_path, self.protein_1, self.chain_1, self.protein_2,
                            self.chain_2, self.spat_prox, self.small_node, self.clust_size, self.magnitude)
        print("Done Tree Class Init")
        progress_callback.emit(f"Initialising {self.protein_1}")
        progress_callback.emit(engine.init_protein(1))
        progress_callback.emit(f"Initialising {self.protein_2}")
        progress_callback.emit(engine.init_protein(2))
        engine.check_sequence_identity()
        engine.check_rmsd()
        while True:
            # If unable to connect to database
            if has_motion_tree == -1 or has_nodes == -1:
                engine.dist_mat_processing()
                progress_callback.emit("Creating Difference Distance Matrix")
                engine.create_distance_difference_matrix()
                progress_callback.emit("Difference Distance Matrix Created")
                progress_callback.emit("Building Motion Tree")
                total_time, num_nodes, protein_str, param_str = engine.run()
                # progress_callback.emit(
                #     f"Motion Tree {protein_str} built. \nTime taken: {total_time}s. \nNumber of Effective Nodes: {nodes}"
                # )
                # print("Checking")
                # print("Result Callback")
            # If database contains motion tree and nodes
            elif has_motion_tree and has_nodes:
                diff_dist_mat, link_mat = get_motion_tree(self.protein_1, self.chain_1, self.protein_2, self.chain_2,
                                                          self.spat_prox)
                nodes = get_nodes(self.protein_1, self.chain_1, self.protein_2, self.chain_2, self.spat_prox,
                                  self.small_node, self.clust_size, self.magnitude)
                if type(diff_dist_mat) == int or type(nodes) == int:
                    continue
                save_results_to_disk(self.output_path, self.protein_1, self.chain_1, self.protein_2, self.chain_2,
                                     self.spat_prox, self.small_node, self.clust_size, self.magnitude, diff_dist_mat, "diff_dist_mat")
                save_results_to_disk(self.output_path, self.protein_1, self.chain_1, self.protein_2, self.chain_2,
                                     self.spat_prox, self.small_node, self.clust_size, self.magnitude, link_mat, "link_mat")
                superimpose_result = write_to_pdb(self.output_path, engine.protein_1, engine.protein_2, self.spat_prox,
                                                  self.small_node, self.clust_size, self.magnitude)
                write_domains_to_pml(self.output_path, engine.protein_1, engine.protein_2, self.spat_prox, self.small_node,
                                     self.clust_size, self.magnitude, nodes)
                write_info_file(self.output_path, engine.protein_1, engine.protein_2, self.spat_prox, self.small_node,
                                self.clust_size, self.magnitude, nodes, superimpose_result)
                protein_str = f"{self.protein_1}_{self.chain_1}_{self.protein_2}_{self.chain_2}"
                param_str = f"sp_{self.spat_prox}_node_{self.small_node}_clust_{self.clust_size}_mag_{self.magnitude}"
                num_nodes = len(nodes)
            # If database does not contain motion tree
            else:
                progress_callback.emit("Processing Distance Matrix")
                engine.dist_mat_processing()
                progress_callback.emit("Creating Difference Distance Matrix")
                diff_dist_mat = np.copy(engine.create_distance_difference_matrix())
                progress_callback.emit("Difference Distance Matrix created")
                progress_callback.emit("Building Motion Tree")
                total_time, num_nodes, protein_str, param_str = engine.run()
                # progress_callback.emit(
                #     f"Motion Tree {protein_str} built. \nTime taken: {total_time}s. \nNumber of Effective Nodes: {nodes}"
                # )
                insert_motion_tree(self.protein_1, self.chain_1, self.protein_2, self.chain_2, self.spat_prox, engine.similarity, total_time, diff_dist_mat, engine.link_mat, has_motion_tree)
                insert_nodes(self.protein_1, self.chain_1, self.protein_2, self.chain_2, self.spat_prox, self.small_node,
                             self.clust_size, self.magnitude, engine.nodes, has_nodes)

            diff_dist_npy, diff_dist_img, motion_tree_img = get_motion_tree_outputs(
                self.output_path, self.protein_1, self.chain_1, self.protein_2, self.chain_2, self.spat_prox,
                self.small_node, self.clust_size, self.magnitude
            )
            progress_callback.emit(
                f"Motion Tree {protein_str} built. \nNumber of Effective Nodes: {num_nodes}"
            )
            return diff_dist_npy, diff_dist_img, motion_tree_img, f"{protein_str}_{param_str}"

    def display_progress(self, msg):
        """
        Displays the progress made by the program when the Motion Tree is being built
        :param msg: The message to be displayed in the Status Bar
        :return:
        """
        print(msg)
        self.progress_label.setStyleSheet("color: 'black';")
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
        self.progress_label.setStyleSheet("color: 'red';")
        self.progress_label.setText(f"Error : {str(msg[1])}")

    def change_image(self, index):
        self.image_output_layout.setCurrentIndex(index)

    def closeEvent(self, event):
        self.closed.emit(self)
        super().closeEvent(event)


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