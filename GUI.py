import sys
from PySide6.QtWidgets import QApplication
from GUIWindow import MainWindow
from GUIWidget import MainWidget


app = QApplication(sys.argv)
widget = MainWidget()
window = MainWindow(app, widget)
window.show()
app.exec()

