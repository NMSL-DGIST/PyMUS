# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_Main.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!
"""

PyMUS: Simulator for virtual experiments on motor unit system

Version 2.0.1

Copyright (C) 2017-Now  Hojeong Kim
Neuromuscular Systems Laboratory
DGIST, Korea

This program is free software: you can redistribute it and/or modify it 
under the terms of the GNU General Public License as published by the Free Software Foundation, 
either version 3 of the License, or any later version.  
This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.  
You should have received a copy of the GNU General Public License along with this program. 
If not, see <http://www.gnu.org/licenses/>.

Please contact us at hojeong.kim03@gmail.com for any inquiries or questions on this program.

"""

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(630, 300)
        MainWindow.setMinimumSize(QtCore.QSize(630, 300))
        MainWindow.setMaximumSize(QtCore.QSize(630, 300))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Calibri"))
        MainWindow.setFont(font)
        MainWindow.setToolButtonStyle(QtCore.Qt.ToolButtonIconOnly)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.textEdit = QtGui.QTextEdit(self.centralwidget)
        self.textEdit.setGeometry(QtCore.QRect(259, 20, 361, 241))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Calibri"))
        font.setPointSize(11)
        font.setBold(False)
        font.setWeight(50)
        self.textEdit.setFont(font)
        self.textEdit.setObjectName(_fromUtf8("textEdit"))
        self.layoutWidget = QtGui.QWidget(self.centralwidget)
        self.layoutWidget.setGeometry(QtCore.QRect(0, 10, 261, 261))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Calibri"))
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.layoutWidget.setFont(font)
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.layoutWidget)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.label = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(11)
        font.setBold(False)
        font.setWeight(50)
        font.setKerning(True)
        self.label.setFont(font)
        self.label.setObjectName(_fromUtf8("label"))
        self.horizontalLayout.addWidget(self.label)
        self.model_comboBox = QtGui.QComboBox(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.model_comboBox.sizePolicy().hasHeightForWidth())
        self.model_comboBox.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Calibri"))
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.model_comboBox.setFont(font)
        self.model_comboBox.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.model_comboBox.setEditable(False)
        self.model_comboBox.setMaxVisibleItems(4)
        self.model_comboBox.setInsertPolicy(QtGui.QComboBox.InsertAtBottom)
        self.model_comboBox.setDuplicatesEnabled(False)
        self.model_comboBox.setFrame(True)
        self.model_comboBox.setObjectName(_fromUtf8("model_comboBox"))
        self.model_comboBox.addItem(_fromUtf8(""))
        self.model_comboBox.addItem(_fromUtf8(""))
        self.model_comboBox.addItem(_fromUtf8(""))
        self.model_comboBox.addItem(_fromUtf8(""))
        self.horizontalLayout.addWidget(self.model_comboBox)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.ParamButton = QtGui.QPushButton(self.layoutWidget)
        self.ParamButton.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.ParamButton.sizePolicy().hasHeightForWidth())
        self.ParamButton.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Calibri"))
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.ParamButton.setFont(font)
        self.ParamButton.setObjectName(_fromUtf8("ParamButton"))
        self.verticalLayout.addWidget(self.ParamButton)
        self.IntegrationSettingButton = QtGui.QPushButton(self.layoutWidget)
        self.IntegrationSettingButton.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.IntegrationSettingButton.sizePolicy().hasHeightForWidth())
        self.IntegrationSettingButton.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Calibri"))
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.IntegrationSettingButton.setFont(font)
        self.IntegrationSettingButton.setObjectName(_fromUtf8("IntegrationSettingButton"))
        self.verticalLayout.addWidget(self.IntegrationSettingButton)
        self.InputSignalButton = QtGui.QPushButton(self.layoutWidget)
        self.InputSignalButton.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.InputSignalButton.sizePolicy().hasHeightForWidth())
        self.InputSignalButton.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Calibri"))
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.InputSignalButton.setFont(font)
        self.InputSignalButton.setObjectName(_fromUtf8("InputSignalButton"))
        self.verticalLayout.addWidget(self.InputSignalButton)
        self.ScopeButton = QtGui.QPushButton(self.layoutWidget)
        self.ScopeButton.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.ScopeButton.sizePolicy().hasHeightForWidth())
        self.ScopeButton.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Calibri"))
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.ScopeButton.setFont(font)
        self.ScopeButton.setAutoDefault(False)
        self.ScopeButton.setDefault(False)
        self.ScopeButton.setFlat(False)
        self.ScopeButton.setObjectName(_fromUtf8("ScopeButton"))
        self.verticalLayout.addWidget(self.ScopeButton)
        self.RunButton = QtGui.QToolButton(self.layoutWidget)
        self.RunButton.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.RunButton.sizePolicy().hasHeightForWidth())
        self.RunButton.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Calibri"))
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.RunButton.setFont(font)
        self.RunButton.setCheckable(True)
        self.RunButton.setChecked(False)
        self.RunButton.setAutoRepeat(False)
        self.RunButton.setObjectName(_fromUtf8("RunButton"))
        self.verticalLayout.addWidget(self.RunButton)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 630, 26))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        self.menuHelp_H = QtGui.QMenu(self.menubar)
        self.menuHelp_H.setObjectName(_fromUtf8("menuHelp_H"))
        MainWindow.setMenuBar(self.menubar)
        self.actionExit = QtGui.QAction(MainWindow)
        self.actionExit.setObjectName(_fromUtf8("actionExit"))
        self.actionAboutThis = QtGui.QAction(MainWindow)
        self.actionAboutThis.setObjectName(_fromUtf8("actionAboutThis"))
        self.actionSave = QtGui.QAction(MainWindow)
        self.actionSave.setObjectName(_fromUtf8("actionSave"))
        self.actionDirOpen = QtGui.QAction(MainWindow)
        self.actionDirOpen.setObjectName(_fromUtf8("actionDirOpen"))
        self.actionModel_Information = QtGui.QAction(MainWindow)
        self.actionModel_Information.setObjectName(_fromUtf8("actionModel_Information"))
        self.menuFile.addAction(self.actionDirOpen)
        self.menuFile.addAction(self.actionSave)
        self.menuFile.addAction(self.actionExit)
        self.menuHelp_H.addAction(self.actionAboutThis)
        self.menuHelp_H.addAction(self.actionModel_Information)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuHelp_H.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QObject.connect(self.actionExit, QtCore.SIGNAL(_fromUtf8("triggered()")), MainWindow.close)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "PyMUS", None))
        self.label.setText(_translate("MainWindow", "  MODEL ", None))
        self.model_comboBox.setItemText(0, _translate("MainWindow", "Please Select !", None))
        self.model_comboBox.setItemText(1, _translate("MainWindow", "Motoneuron", None))
        self.model_comboBox.setItemText(2, _translate("MainWindow", "Muscle Fibers", None))
        self.model_comboBox.setItemText(3, _translate("MainWindow", "Motor Unit", None))
        self.ParamButton.setText(_translate("MainWindow", "Model Parameter Settings", None))
        self.IntegrationSettingButton.setText(_translate("MainWindow", "Simulation Condition Settings", None))
        self.InputSignalButton.setText(_translate("MainWindow", "Input Signal Settings", None))
        self.ScopeButton.setText(_translate("MainWindow", "Output Signal Settings", None))
        self.RunButton.setText(_translate("MainWindow", "Run", None))
        self.menuFile.setTitle(_translate("MainWindow", "File(F)", None))
        self.menuHelp_H.setTitle(_translate("MainWindow", "Help(H)", None))
        self.actionExit.setText(_translate("MainWindow", "Exit", None))
        self.actionAboutThis.setText(_translate("MainWindow", "About", None))
        self.actionSave.setText(_translate("MainWindow", "Save as..", None))
        self.actionDirOpen.setText(_translate("MainWindow", "Open", None))
        self.actionModel_Information.setText(_translate("MainWindow", "Model Information", None))

import icons_rc
