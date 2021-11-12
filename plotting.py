from PyQt5 import QtGui  # (the example applies equally well to PySide2)
import pyqtgraph as pg

## Always start by initializing Qt (only once per application)
app = QtGui.QApplication([])

## Define a top-level widget to hold everything
w = QtGui.QWidget()

## Create some widgets to be placed inside
btn = QtGui.QPushButton('press me')
text = QtGui.QLineEdit('enter text')
listw = QtGui.QListWidget()
chart = pg.PlotWidget()
chart.plot([1,2,3], [3,2,1])

## Create a grid layout to manage the widgets size and position
layout = QtGui.QGridLayout()
w.setLayout(layout)

## Add widgets to the layout in their proper positions
layout.addWidget(listw, 0, 0)  # list widget goes in bottom-left
layout.addWidget(chart, 0, 1)  # plot goes on right side, spanning 3 rows

## Display the widget as a new window
#w.show()

## Start the Qt event loop
#app.exec_()