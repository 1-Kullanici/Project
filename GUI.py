from PyQt6.QtWidgets import QApplication

import GUI_background
import sys # Only needed for access to command line arguments


#######################
# You need one (and only one) QApplication instance per application.
# Pass in sys.argv to allow command line arguments for your app.
# If you know you won't use command line arguments QApplication([]) works too.
# Subclass QMainWindow to customize your application's main window

app = QApplication(sys.argv)

window = GUI_background.MainWindow()
window.setWindowTitle("Numerical Simulation Tool for Diffusion by Ibrahim Furkan Tezcan  |  v1.1-stable  |  Release: 10.Jun.2024")
window.show()

app.exec()
