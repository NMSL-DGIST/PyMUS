# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_AboutThis2.ui'
#
# Created: Mon Mar 12 17:55:17 2018
#      by: PyQt4 UI code generator 4.10.4
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

class Ui_AboutThisWindow(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(846, 377)
        Dialog.setMinimumSize(QtCore.QSize(380, 120))
        Dialog.setLayoutDirection(QtCore.Qt.LeftToRight)
        Dialog.setSizeGripEnabled(False)
        Dialog.setModal(False)
        self.ok_pushButton = QtGui.QPushButton(Dialog)
        self.ok_pushButton.setGeometry(QtCore.QRect(730, 330, 91, 31))
        self.ok_pushButton.setObjectName(_fromUtf8("ok_pushButton"))
        self.verticalLayoutWidget = QtGui.QWidget(Dialog)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(160, 20, 671, 291))
        self.verticalLayoutWidget.setObjectName(_fromUtf8("verticalLayoutWidget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setMargin(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label = QtGui.QLabel(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.label.setFrameShape(QtGui.QFrame.NoFrame)
        self.label.setLineWidth(1)
        self.label.setTextFormat(QtCore.Qt.LogText)
        self.label.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label.setWordWrap(False)
        self.label.setMargin(0)
        self.label.setObjectName(_fromUtf8("label"))
        self.verticalLayout.addWidget(self.label)
        self.label_5 = QtGui.QLabel(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_5.sizePolicy().hasHeightForWidth())
        self.label_5.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.label_5.setFont(font)
        self.label_5.setFrameShape(QtGui.QFrame.NoFrame)
        self.label_5.setLineWidth(1)
        self.label_5.setTextFormat(QtCore.Qt.LogText)
        self.label_5.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_5.setWordWrap(False)
        self.label_5.setMargin(0)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.verticalLayout.addWidget(self.label_5)
        self.label_2 = QtGui.QLabel(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.label_2.setFont(font)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.verticalLayout.addWidget(self.label_2)
        self.label_3 = QtGui.QLabel(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_3.sizePolicy().hasHeightForWidth())
        self.label_3.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.label_3.setFont(font)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.verticalLayout.addWidget(self.label_3)
        self.label_4 = QtGui.QLabel(self.verticalLayoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_4.sizePolicy().hasHeightForWidth())
        self.label_4.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.label_4.setFont(font)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.verticalLayout.addWidget(self.label_4)
        self.label_6 = QtGui.QLabel(Dialog)
        self.label_6.setGeometry(QtCore.QRect(20, 20, 121, 131))
        self.label_6.setText(_fromUtf8(""))
        self.label_6.setPixmap(QtGui.QPixmap(_fromUtf8("resources/PyMUS_image.png")))
        self.label_6.setObjectName(_fromUtf8("label_6"))

        self.retranslateUi(Dialog)
        QtCore.QObject.connect(self.ok_pushButton, QtCore.SIGNAL(_fromUtf8("clicked()")), Dialog.close)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "About PyMUS", None))
        self.ok_pushButton.setText(_translate("Dialog", "OK", None))
        self.label.setText(_translate("Dialog", "PyMUS: Simulator for virtual experiments on motor unit system", None))
        self.label_5.setText(_translate("Dialog", "Version 2.0.1", None))
        self.label_2.setText(_translate("Dialog", "Copyright (C) 2017-Now  Hojeong Kim\n"
"Neuromuscular Systems Laboratory\n"
"DGIST, Korea\n"
"", None))
        self.label_3.setText(_translate("Dialog", "This program is free software: \n"
"you can redistribute it and/or modify it under the terms of the GNU General Public License\n"
"as published by the Free Software Foundation, either version 3 of the License, or any later version.\n"
"This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; \n"
"without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. \n"
"See the GNU General Public License for more details.  \n"
"You should have received a copy of the GNU General Public License along with this program. \n"
"If not, see <http://www.gnu.org/licenses/>.\n"
"", None))
        self.label_4.setText(_translate("Dialog", "Please contact us at hojeong.kim03@gmail.com for any inquiries or questions on this program.", None))

