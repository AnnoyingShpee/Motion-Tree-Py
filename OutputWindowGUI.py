import sys
import traceback
from PySide6.QtCore import Signal, QObject, QRunnable, Slot
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import QPushButton, QVBoxLayout, QHBoxLayout, QStackedLayout, QWidget, QLabel
from MotionTree import MotionTree
from FileMngr import get_motion_tree_outputs, save_results_to_disk, write_info_file, write_to_pdb, write_domains_to_pml, \
    write_to_pdb_dyndom
from DataMngr import conn, check_motion_tree_exists, get_motion_tree, insert_motion_tree, check_nodes_exist, \
    get_nodes, insert_nodes, check_protein_pair_exists, get_protein_pair, insert_protein_pair


class OutputWindow(QWidget):
    closed = Signal(QWidget)

    def __init__(self, thread_pool, screensize, input_path, output_path, protein_1, chain_1=None, protein_2=None,
                 chain_2=None, spat_prox=7, small_node=5, clust_size=30, magnitude=5):
        super().__init__()

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
        if self.protein_2 is None:
            self.setWindowTitle(protein_1)
            protein_1_label = QLabel(f"DynDom Code : {self.protein_1}")
            combi_layout.addWidget(protein_1_label)
        else:
            self.setWindowTitle(f"{protein_1}_{chain_1}_{protein_2}_{chain_2}")
            protein_1_label = QLabel(f"Protein 1 : {self.protein_1} ({self.chain_1})")
            protein_2_label = QLabel(f"Protein 2 : {self.protein_2} ({self.chain_2})")
            combi_layout.addWidget(protein_1_label)
            combi_layout.addWidget(protein_2_label)

        sp_label = QLabel(f"Spatial Proximity : {self.spat_prox}")
        node_label = QLabel(f"Small Domain Size : {self.small_node}")
        clust_label = QLabel(f"Cluster Size : {self.clust_size}")
        mag_label = QLabel(f"Magnitude : {magnitude}")
        self.progress_label = QLabel(f"Progress : Running")

        combi_layout.addWidget(sp_label)
        combi_layout.addWidget(node_label)
        combi_layout.addWidget(clust_label)
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
        has_protein_pair, has_motion_tree, has_nodes = -1, -1, -1
        chain_1 = self.chain_1 if self.chain_1 is not None else "A"
        protein_2 = self.protein_2 if self.protein_2 is not None else self.protein_1
        chain_2 = self.chain_2 if self.chain_2 is not None else "B"
        is_dyndom = False if self.protein_2 is not None else True
        if conn is not None:
            has_protein_pair = check_protein_pair_exists(self.protein_1, chain_1, protein_2, chain_2)
            has_motion_tree = check_motion_tree_exists(self.protein_1, chain_1, protein_2, chain_2, self.spat_prox)
            has_nodes = check_nodes_exist(self.protein_1, chain_1, protein_2, chain_2,
                                          self.spat_prox, self.small_node, self.clust_size, self.magnitude)
        engine = MotionTree(self.input_path, self.output_path, self.protein_1, self.chain_1, self.protein_2,
                            self.chain_2, self.spat_prox, self.small_node, self.clust_size, self.magnitude, is_dyndom)
        print("Done Tree Class Init")
        progress_callback.emit(f"Initialising {self.protein_1}")
        progress_callback.emit(engine.init_protein(1))
        progress_callback.emit(f"Initialising {protein_2}")
        progress_callback.emit(engine.init_protein(2))
        engine.preprocessing()

        while True:
            # If unable to connect to database
            if has_protein_pair == -1 or has_motion_tree == -1 or has_nodes == -1:
                engine.dist_mat_processing()
                progress_callback.emit("Creating Difference Distance Matrix")
                engine.create_distance_difference_matrix()
                progress_callback.emit("Difference Distance Matrix Created")
                progress_callback.emit("Building Motion Tree")
                total_time, num_nodes, protein_str, param_str = engine.run()
            # If database contains motion tree and nodes
            elif has_protein_pair is True and has_motion_tree is True and has_nodes is True:
                rmsd, diff_dist_mat = get_protein_pair(self.protein_1, chain_1, protein_2, chain_2)
                link_mat = get_motion_tree(self.protein_1, chain_1, protein_2, chain_2, self.spat_prox)
                nodes = get_nodes(self.protein_1, chain_1, protein_2, chain_2, self.spat_prox,
                                  self.small_node, self.clust_size, self.magnitude)
                if type(diff_dist_mat) == int or type(link_mat) == int or type(nodes) == int:
                    has_protein_pair, has_motion_tree, has_nodes = -1, -1, -1
                    continue
                save_results_to_disk(self.output_path, self.protein_1, self.chain_1, self.protein_2, self.chain_2,
                                     self.spat_prox, self.small_node, self.clust_size, self.magnitude, diff_dist_mat, "diff_dist_mat")
                save_results_to_disk(self.output_path, self.protein_1, self.chain_1, self.protein_2, self.chain_2,
                                     self.spat_prox, self.small_node, self.clust_size, self.magnitude, link_mat, "dendrogram")
                if is_dyndom:
                    write_to_pdb_dyndom(self.output_path, engine.protein_1, engine.protein_2, self.spat_prox,
                                        self.small_node, self.clust_size, self.magnitude)
                else:
                    write_to_pdb(self.output_path, engine.protein_1, engine.protein_2, self.spat_prox,
                                 self.small_node, self.clust_size, self.magnitude)
                if self.protein_2 is None:
                    is_dyndom = True
                    protein_str = self.protein_1
                else:
                    is_dyndom = False
                    protein_str = f"{self.protein_1}_{self.chain_1}_{self.protein_2}_{self.chain_2}"
                write_domains_to_pml(self.output_path, engine.protein_1, engine.protein_2, self.spat_prox, self.small_node,
                                     self.clust_size, self.magnitude, nodes, is_dyndom)
                write_info_file(self.output_path, engine.protein_1, engine.protein_2, self.spat_prox, self.small_node,
                                self.clust_size, self.magnitude, nodes, rmsd, is_dyndom)
                param_str = f"sp_{self.spat_prox}_node_{self.small_node}_clust_{self.clust_size}_mag_{self.magnitude}"
                num_nodes = len(nodes)
            # If database does not contain motion tree
            else:
                progress_callback.emit("Processing Distance Matrices")
                engine.dist_mat_processing()
                progress_callback.emit("Creating Difference Distance Matrix")

                if has_protein_pair is True:
                    rmsd, diff_dist_mat = get_protein_pair(self.protein_1, chain_1, protein_2, chain_2)
                    engine.create_distance_difference_matrix(diff_dist_mat)
                    is_fail_1 = False
                else:
                    engine.create_distance_difference_matrix(None)
                    is_fail_1 = insert_protein_pair(self.protein_1, chain_1, protein_2, chain_2, engine.rmsd, engine.diff_dist_mat_init)
                progress_callback.emit("Difference Distance Matrix created. Building Motion Tree")
                total_time, num_nodes, protein_str, param_str = engine.run()
                is_fail_2 = insert_motion_tree(self.protein_1, chain_1, protein_2, chain_2,
                                               self.spat_prox, engine.similarity, total_time, engine.link_mat, has_motion_tree)
                is_fail_3 = insert_nodes(self.protein_1, chain_1, protein_2, chain_2, self.spat_prox,
                                         self.small_node, self.clust_size, self.magnitude, engine.nodes, has_nodes)

                if is_fail_1 or is_fail_2 or is_fail_3:
                    has_protein_pair, has_motion_tree, has_nodes = -1, -1, -1
                    continue

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
        self.widgets["diff_dist_mat_img"].setScaledContents(True)
        motion_tree_img = QPixmap(outputs[2])
        self.widgets["motion_tree_img"].setPixmap(motion_tree_img)
        self.widgets["motion_tree_img"].setScaledContents(True)
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