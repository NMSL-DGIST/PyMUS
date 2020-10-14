# -*- coding: utf-8 -*-
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

import matplotlib as mpl
mpl.use('qt4Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
import sys
import os
import warnings
import copy
from Motorunit import *
from PyQt4.QtGui import *
from PyQt4.QtCore import *
from ui_Main import Ui_MainWindow as Ui_M
from ui_InitialValue import Ui_InitialValueWindow as Ui_ID
from ui_Oscilloscope import Ui_OscilloscopeWindow as Ui_OD
from ui_Parameter import Ui_ParameterSettingWindow as Ui_PD
from ui_SignalGenerator import Ui_SignalGeneratorWindow as Ui_SGD
from ui_AboutThis import Ui_AboutThisWindow as Ui_AT

# Main Window class
class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.PSW = ParameterSettingWindow(self) 
        self.ISW = IntegrationSettingWindow(self)
        self.SGW = SignalGeneratorWindow(self)
        self.OW = OscilloscopeWindow(self)
        self.OSC = Oscilloscope()
        self.ATW = AboutThisWindow()
        self.ModelType = ''
        self.uim = Ui_M()
        self.uim.setupUi(self)
        self.uim.ParamButton.setEnabled(False)
        self.uim.IntegrationSettingButton.setEnabled(False)
        self.uim.InputSignalButton.setEnabled(False)
        self.uim.ScopeButton.setEnabled(False)
        self.uim.RunButton.setEnabled(False)
        self.uim.actionSave.setEnabled(False)      
        self.uim.textEdit.setReadOnly(True)
        
        # GUI operation-function
        self.uim.actionDirOpen.triggered.connect(self.openResultDirectory)
        self.uim.actionSave.triggered.connect(self.saveToFile)
        self.uim.actionModel_Information.connect(self.uim.actionModel_Information, SIGNAL("triggered()"), self.openModelInfo)
        self.uim.actionModel_Information.triggered.connect(self.openModelInfo)
        self.uim.actionAboutThis.connect(self.uim.actionAboutThis, SIGNAL("triggered()"), self.openAboutThisWindow)
        self.uim.model_comboBox.connect(self.uim.model_comboBox, SIGNAL('activated(QString)'), self.importModel)
        self.uim.InputSignalButton.clicked.connect(self.openInputSignalWindow)
        self.uim.ParamButton.clicked.connect(self.openParameterWindow)
        self.uim.IntegrationSettingButton.clicked.connect(self.openIntegrationWindow)
        self.uim.ScopeButton.clicked.connect(self.openOscilloscopeWindow)
        self.uim.RunButton.clicked.connect(self.runSimulation)
    
    def importModel(self, text):
        # Question box
        if(self.ModelType != ''):
            result = QMessageBox.question(self, 'Confirm Change...', 'Are you sure you want to change the model ?', QMessageBox.Yes| QMessageBox.No)
            if(result == QMessageBox.No):
                if(self.ModelType == 'Motoneuron'):                
                    self.uim.model_comboBox.setCurrentIndex(1)
                elif(self.ModelType == 'Muscle Fibers'):
                    self.uim.model_comboBox.setCurrentIndex(2)
                elif(self.ModelType == 'Motor Unit'):
                    self.uim.model_comboBox.setCurrentIndex(3)
                return
        
        self.ModelType = text
        # No model
        if(self.ModelType == 'Please Select !'):
            self.ModelType = ''
            self.delete()
            self.uim.ParamButton.setEnabled(False)
            self.uim.IntegrationSettingButton.setEnabled(False)
            self.uim.InputSignalButton.setEnabled(False)
            self.uim.ScopeButton.setEnabled(False)
            self.uim.RunButton.setEnabled(False)
            self.uim.RunButton.setStyleSheet("alternate-background-color: rgb(0, 0, 0)")
            self.uim.actionSave.setEnabled(False)
            self.ATW = AboutThisWindow()
            self.uim.textEdit.clear()
            self.setTextEdit("No Model selected.")
        # MN, MF, NU model
        else:
            self.delete()
            self.uim.ParamButton.setEnabled(True)
            self.uim.IntegrationSettingButton.setEnabled(True)
            self.uim.InputSignalButton.setEnabled(True)
            self.uim.ScopeButton.setEnabled(True)
            self.uim.ScopeButton.setEnabled(True)
            self.uim.RunButton.setEnabled(True)
            self.uim.RunButton.setStyleSheet("background-color: rgb(155, 255, 172)")
            self.uim.actionSave.setEnabled(False)
            self.OSC = Oscilloscope()
            self.ISW = IntegrationSettingWindow(self)        
            self.PSW = ParameterSettingWindow(self)     
            self.OW = OscilloscopeWindow(self) 
            self.SGW = SignalGeneratorWindow(self)
            self.ATW = AboutThisWindow()
            
            # Model creation, default value setting
            if(self.ModelType == 'Motoneuron' or self.ModelType == 'Motor Unit'):
                self.MN = MotoNeuron(1)
                if(self.ModelType == 'Motoneuron'):
                    self.SGW.displayValue()
                    self.MN.setModelParam(self.PSW.mn_setValue, False, False)
                    self.MN.setIntegrationEnv(self.ISW.t_start, self.ISW.t_stop, self.ISW.t_dt, self.ISW.t_pt)
                    self.MN.setInitialValues(self.ISW.mn_setValue)
                    self.MN.setInputSignal(self.SGW.ISG.signalType, self.SGW.ISG.iValue, self.SGW.ISG.pValue, self.SGW.ISG.Is_0, self.SGW.ISG.heav_param)
                    self.MN.setSynConSignal(True, True, True, True, self.SGW.s_SCSG.e_times, self.SGW.s_SCSG.i_times, self.SGW.d_SCSG.e_times, self.SGW.d_SCSG.i_times, self.SGW.s_SCSG.G_e, self.SGW.s_SCSG.G_i, self.SGW.d_SCSG.G_e, self.SGW.d_SCSG.G_i)
                    self.OSC.setModel(self.MN)
            
            if(self.ModelType == 'Muscle Fibers' or self.ModelType == 'Motor Unit'):
                self.MF = MuscleFibers(1)
                if(self.ModelType == 'Muscle Fibers'):
                    self.MF.setModelParam(self.PSW.mf_setValue)
                    self.MF.setIntegrationEnv(self.ISW.t_start, self.ISW.t_stop, self.ISW.t_dt, self.ISW.t_pt)
                    self.MF.setInitialValues(self.ISW.mf_setValue)
                    self.MF.setSpikeSignal(self.SGW.SSG.spike, self.SGW.SSG.spike_idx, self.SGW.SSG.SpikeTimes)
                    self.MF.setXmSignal(self.SGW.XSG.signalType, self.SGW.XSG.times, self.SGW.XSG.xm)
                    self.OSC.setModel(self.MF)
                
                elif(self.ModelType == 'Motor Unit'):
                    self.MU = Motorunit(1, self.MN, self.MF)
                    self.SGW.displayValue()
                    self.MU.setModelParam(self.PSW.mn_setValue, self.PSW.mf_setValue, False, False)
                    self.MU.setIntegrationEnv(self.ISW.t_start, self.ISW.t_stop, self.ISW.t_dt, self.ISW.t_pt)
                    self.MU.setInitialValues(self.ISW.mn_setValue, self.ISW.mf_setValue)
                    self.MU.setInputSignal(self.SGW.ISG.signalType, self.SGW.ISG.iValue, self.SGW.ISG.pValue, self.SGW.ISG.Is_0, self.SGW.ISG.heav_param)
                    self.MU.setSynConSignal(True, True, True, True, self.SGW.s_SCSG.e_times, self.SGW.s_SCSG.i_times, self.SGW.d_SCSG.e_times, self.SGW.d_SCSG.i_times, self.SGW.s_SCSG.G_e, self.SGW.s_SCSG.G_i, self.SGW.d_SCSG.G_e, self.SGW.d_SCSG.G_i)
                    self.MU.setXmSignal(self.SGW.XSG.signalType, self.SGW.XSG.times, self.SGW.XSG.xm)
                    self.OSC.setModel(self.MU)
            self.uim.textEdit.clear()            
            self.setTextEdit("["+self.ModelType+"] Model created.")                    

    def openResultDirectory(self):
        curPath = os.path.dirname( os.path.abspath( sys.argv[0] ) )
        if(sys.platform == 'win32'):
            os.startfile(curPath)
    
    def saveToFile(self):
        curPath = os.path.dirname( os.path.abspath( sys.argv[0] ) )
        t=time.localtime()
        tm_year='{:0>4}'.format(t.tm_year)
        tm_mon='{:0>2}'.format(t.tm_mon)
        tm_mday='{:0>2}'.format(t.tm_mday)
        tm_hour='{:0>2}'.format(t.tm_hour)
        tm_min='{:0>2}'.format(t.tm_min)
        tm_sec='{:0>2}'.format(t.tm_sec)
        timePath='_'+str(tm_year)+str(tm_mon)+str(tm_mday)
        timeName='_'+str(tm_hour)+str(tm_min)+str(tm_sec)

        # file name
        fileName=str(self.OSC.cell.cellType)+timePath+timeName+'.csv'
        filePath = QFileDialog.getSaveFileName(self, "Save Result data in a File. (Must be a csv format)", curPath+'\\'+fileName, "CSV Files (*.csv)")
        if(filePath==''):
            return
            
        # save to file
        if(self.ModelType=="Motoneuron"):            
            df1=pd.DataFrame(self.OSC.mn_ResultArrays)
            cols = ['Time','Is', 'V_soma', '[Ca]_soma','E_Ca_soma','I_Naf_soma','m_Naf_soma','h_Naf_soma','I_Nap_soma','m_Nap_soma','I_Kdr_soma','n_Kdr_soma','I_Kca_soma','I_Can_soma','m_Can_soma','h_Can_soma','I_H_soma','m_H_soma','V_dend','[Ca]_dend','E_Ca_dend','I_Cal_dend','l_Cal_dend','I_Naf_dend','m_Naf_dend','h_Naf_dend','I_Nap_dend','m_Nap_dend','I_Kdr_dend','n_Kdr_dend','I_Kca_dend','I_Can_dend','m_Can_dend','h_Can_dend','I_H_dend','m_H_dend','I_esyn_soma','G_esyn_soma','I_isyn_soma','G_isyn_soma','I_esyn_dend','G_esyn_dend','I_isyn_dend','G_isyn_dend']
        elif(self.ModelType=='Muscle Fibers'):            
            df1=pd.DataFrame(self.OSC.mf_ResultArrays)
            cols = ['Time','F','A','Am','A_tilde','B','CaSP','CaSPB','CaSPT','CaSR','CaSRCS','Cs','R','Spike','T','Vm','XCE','Xm']
        elif(self.ModelType=='Motor Unit'):                
            df1=pd.DataFrame(self.OSC.mu_ResultArrays)
            cols = ['Time','Is','FiringRate','V_soma','[Ca]_soma','E_Ca_soma','I_Naf_soma','m_Naf_soma','h_Naf_soma','I_Nap_soma','m_Nap_soma','I_Kdr_soma','n_Kdr_soma','I_Kca_soma','I_Can_soma','m_Can_soma','h_Can_soma','I_H_soma','m_H_soma','V_dend','[Ca]_dend','E_Ca_dend','I_Cal_dend','l_Cal_dend','I_Naf_dend','m_Naf_dend','h_Naf_dend','I_Nap_dend','m_Nap_dend','I_Kdr_dend','n_Kdr_dend','I_Kca_dend','I_Can_dend','m_Can_dend','h_Can_dend','I_H_dend','m_H_dend','I_esyn_soma','G_esyn_soma','I_isyn_soma','G_isyn_soma','I_esyn_dend','G_esyn_dend','I_isyn_dend','G_isyn_dend','A','Am','A_tilde','B','CaSP','CaSPB','CaSPT','CaSR','CaSRCS','Cs','F','R','Spike','T','Vm','XCE','Xm']
        df = df1.reindex(columns = cols)
        df.to_csv(filePath)
        self.setTextEdit("["+self.ModelType+"] Simulation data saved.")
       
    def openModelInfo(self):
        if(sys.platform == 'win32'):
            os.startfile("model_equations.pdf")
    
    def openAboutThisWindow(self):
        self.ATW.show()
        self.ATW.raise_()
        self.ATW.activateWindow()

    def setTextEdit(self, message):
        cursor = '>>> '
        self.uim.textEdit.append(cursor+message)
    
    def openInputSignalWindow(self):            
        self.SGW.displayValue()
        self.SGW.show()
        self.SGW.raise_()
        self.SGW.activateWindow()

    def runSimulation(self):
        if(self.uim.RunButton.isChecked()):        
            if(self.ModelType == 'Motoneuron'):
                model = self.MN
                scopeList1 = self.OW.mn_ScopeList
                scopeList2 = self.OW.mn_ScopeList
            elif(self.ModelType == 'Muscle Fibers'):
                model = self.MF
                scopeList1 = self.OW.mf_ScopeList
                scopeList2 = self.OW.mf_ScopeList
            elif(self.ModelType == 'Motor Unit'):
                model = self.MU
                scopeList1 = self.OW.mn_ScopeList
                scopeList2 = self.OW.mf_ScopeList
            
            if(self.OW.displayResult == self.OW.uid.radioButton_Individual):
                self.OSC.setScope(list(scopeList1), list(scopeList2), 'Individual')
            elif(self.OW.displayResult == self.OW.uid.radioButton_Combined):
                self.OSC.setScope(list(scopeList1), list(scopeList2), 'Combined')
            
            model.cellState = 'Normal'
            self.uim.model_comboBox.setEnabled(False)
            self.uim.RunButton.setStyleSheet("background-color: rgb(255, 106, 106)")
            self.uim.RunButton.setText("Stop")
            QApplication.processEvents()
            self.setTextEdit("["+self.ModelType+"] Simulation started.")
            
            # run simulation
            try:
                warnings.filterwarnings("error")
                self.OSC.run()
            
            except:
                self.setTextEdit("[Error] Simulation stopped due to integrator error.")
            finally:
                warnings.filterwarnings("default")
            
            # end simulation
            if(model.cellState != 'Stop'):                    
                self.uim.RunButton.setStyleSheet("background-color: rgb(155, 255, 172)")
                self.uim.RunButton.setText("Run")
                self.uim.RunButton.setChecked(False)
                self.setTextEdit("["+self.ModelType+"] Simulation finished.")
            self.setTextEdit("["+self.ModelType+"] Elapsed real time for simulation: "+str(model.simulTime)+" s")
            self.uim.model_comboBox.setEnabled(True)
            self.uim.actionSave.setEnabled(True)
                
        # stop simulation
        else:
            if(self.ModelType == 'Motoneuron'):
                model = self.MN
            elif(self.ModelType == 'Muscle Fibers'):
                model = self.MF
            elif(self.ModelType == 'Motor Unit'):
                model = self.MU
            model.cellState = 'Stop'
            self.uim.model_comboBox.setEnabled(True)
            self.uim.RunButton.setStyleSheet("background-color: rgb(155, 255, 172)")
            self.uim.RunButton.setText('Run')
            self.setTextEdit("["+self.ModelType+"] Simulation finished.")

    def openParameterWindow(self):
        self.PSW.displayValue()
        self.PSW.show()
        self.PSW.raise_()
        self.PSW.activateWindow()

    def openIntegrationWindow(self):
        self.ISW.displayValue()
        self.ISW.show()
        self.ISW.raise_()
        self.ISW.activateWindow()
    
    def openOscilloscopeWindow(self):
        self.OW.displayScope()
        self.OW.show()
        self.OW.raise_()
        self.OW.activateWindow()
   
    def delete(self):    
        # delete instances
        plt.close('all')
        sys.exc_traceback = sys.last_traceback = None
        if(hasattr(self, 'OSC')):
            del self.OSC
        if(hasattr(self, 'ISW')):
            del self.ISW
        if(hasattr(self, 'PSW')):
            del self.PSW
        if(hasattr(self, 'OW')):
            del self.OW
        if(hasattr(self, 'SGW')):
            del self.SGW
        if(hasattr(self, 'ATW')):
            del self.ATW
        if(hasattr(self, 'MN')):
            del self.MN
        if(hasattr(self, 'MF')):
            del self.MF
        if(hasattr(self, 'MU')):
            del self.MU

    def closeEvent(self, event):
        # Question box
        result = QMessageBox.question(self, 'Confirm Exit...', 'Are you sure you want to exit ?', QMessageBox.Yes| QMessageBox.No)
        event.ignore()

        if(result == QMessageBox.Yes):
            self.delete()
            event.accept()

# About PyMUS information Window class
class AboutThisWindow(QDialog):
    def __init__(self):
        super(AboutThisWindow, self).__init__()
        self.uiat = Ui_AT()
        self.uiat.setupUi(self)
        
# Model Parameter Setting Window class       
class ParameterSettingWindow(QDialog):
    def __init__(self, MW):
        super(ParameterSettingWindow, self).__init__()
        self.MW = MW
        self.init = False
        self.const_sEca = False
        self.const_dEca = False
        # subscript, superscript, Greek letter with unicode and unit label
        alpha = u'\u03B1'
        beta = u'\u03B2'
        gamma = u'\u03B3'
        ohm = u'\u2126'
        mu = u'\u03BC'
        dot = u'\u22C5'
        super_2 = u'\u00B2'
        super_2_p = super_2+u'\u207A'
        super_minus_1 = u'\u207B\u00B9' 
        M_ohm = ' (M'+ohm+')'
        mS_cm_2 = ' (mS/cm'+super_2+')'
        mV = ' (mV)'
        mM = ' (mM)'
        ms_1 = ' (ms'+ super_minus_1+')'
        mol_uC_cm_2 = ' (mol/'+mu+'C/cm'+super_2+')'
        mV_ms_1 = ' (mV'+dot+'ms)'+super_minus_1
        sub_i = u'\u1D62'
        sub_0 = u'\u2080'
        sub_1 = u'\u2081'
        sub_2 = u'\u2082'
        sub_3 = u'\u2083'
        sub_4 = u'\u2084'
        sub_5 = u'\u2085'
        sub_6 = u'\u2086'
        M_1 = ' (M'+super_minus_1+')'
        M_ms_1 = ' (M ms'+super_minus_1+')'
        M_1ms_1 = ' (M'+super_minus_1+'ms'+super_minus_1+')'
        ms_1 = ' (ms'+super_minus_1+')'
        ms = ' (ms)'
        s = ' (s)'
        mm = ' (mm)'
        m_s = ' (m/s)'
        mm_1 = ' (mm'+super_minus_1+')' 
        n = ' (N)'
        mm_s_1 = ' (mm s'+super_minus_1+')'
        tau = u'\u03C4' 
        pi = u'\u03C6'
        
        # subscript with HTML
        sub_N = '<sub>N</sub>'
        sub_m = '<sub>m</sub>'
        sub_SD = '<sub>SD</sub>'
        super_DC = '<sup>DC</sup>'
        sub_DS = '<sub>DS</sub>'
        super_AC = '<sup>AC</sup>'
        sub_Na = '<sub>Na</sub>'
        sub_Naf = '<sub>Naf</sub>'
        sub_Nap = '<sub>Nap</sub>'
        sub_K = '<sub>K</sub>'
        sub_Kdr = '<sub>Kdr</sub>'
        sub_KCa = '<sub>K(Ca)</sub>' 
        sub_Can = '<sub>Can</sub>' 
        sub_H = '<sub>H</sub>' 
        sub_esyn = '<sub>esyn</sub>' 
        sub_isyn = '<sub>isyn</sub>' 
        sub_Cal = '<sub>Cal</sub>'     
        sub_Leak_S = '<sub>Leak,S</sub>'
        sub_Leak_D = '<sub>Leak,D</sub>'
        sub_Ca = '<sub>Ca</sub>'
        sub_d = '<sub>d</sub>'
        sub_o = '<sub>o</sub>'
        sub_max = '<sub>max</sub>'
        K_SE = 'K<sub>SE</sub>'
        sub_axon =  '<sub>axon</sub>'
        
        ## initialization
        # set value
        curPath = os.path.dirname( os.path.abspath( sys.argv[0] ) )
        param = pd.read_csv(curPath + '/parameters/MN_Parameters/MN_Parameters.csv')
        setValue = range(len(param.index))
        for i in range(len(param.index)):
            setValue[i] = float(param.Value[i])
        self.mn_idx = [0, 9, 14, 30, 40, 48, 51, 58, 63, 64, 65, 70, 74, 90, 100, 108, 111, 118, 123, 124, 125]
        self.mn_setValue_num = 125
        iter_num = len(self.mn_idx)
        idx = self.mn_idx
        lstSlice_arr = range(len(self.mn_idx)-1)
        for i in range(1, iter_num):
            lstSlice_arr[i-1] = setValue[idx[i-1]:idx[i]]
        self.mn_setValue = lstSlice_arr
        
        param = pd.read_csv(curPath + '/parameters/MF_Parameters/MF_Parameters.csv')
        setValue = range(len(param.index))
        for i in range(len(param.index)):
            setValue[i] = float(param.Value[i])
        self.mf_idx = [0, 15, 20, 28] 
        self.mf_setValue_num = 28
        iter_num = len(self.mf_idx)
        idx = self.mf_idx
        lstSlice_arr = range(len(self.mf_idx)-1)
        for i in range(1, iter_num):
            lstSlice_arr[i-1] = setValue[idx[i-1]:idx[i]]
        self.mf_setValue = lstSlice_arr
        
        # default value
        self.mn_defaultValue = copy.deepcopy(self.mn_setValue) 
        self.mf_defaultValue = copy.deepcopy(self.mf_setValue)
        
        # previous G
        iter_num = len(self.mn_idx)-1
        idx = self.mn_idx
        k = 0
        self.pre_G = range(len(self.mn_idx)-1)
        for i in range(iter_num):        
            self.pre_G[i] = self.mn_setValue[i][0]  
        
        ## GUI
        self.uid = Ui_PD()
        self.uid.setupUi(self)
        self.mn_tableWidget = [self.uid.tableWidget_share,
                               self.uid.tableWidget_Sca,
                               self.uid.tableWidget_Snaf,
                               self.uid.tableWidget_Snap,
                               self.uid.tableWidget_Skdr,
                               self.uid.tableWidget_Skca,
                               self.uid.tableWidget_Scan,
                               self.uid.tableWidget_Sh,
                               self.uid.tableWidget_Sesyn,
                               self.uid.tableWidget_Sisyn,
                               self.uid.tableWidget_Dca,
                               self.uid.tableWidget_Dcal,
                               self.uid.tableWidget_Dnaf,
                               self.uid.tableWidget_Dnap,
                               self.uid.tableWidget_Dkdr,
                               self.uid.tableWidget_Dkca,
                               self.uid.tableWidget_Dcan,
                               self.uid.tableWidget_Dh,
                               self.uid.tableWidget_Desyn,
                               self.uid.tableWidget_Disyn]                       
        self.mf_tableWidget = [self.uid.tableWidget_m1, 
                               self.uid.tableWidget_m2, 
                               self.uid.tableWidget_m3]
        self.uid.mn_lineEdit.setReadOnly(True)
        self.uid.mf_lineEdit.setReadOnly(True)
        self.channel_enable = {'Snaf' : True,
                               'Snap' : True,
                               'Skdr' : True,
                               'Skca' : True,
                               'Scan' : True,
                               'Sh' : True,
                               'Sesyn' : True,
                               'Sisyn' : True,
                               'Dcal' : True,
                               'Dnaf' : True,
                               'Dnap' : True,
                               'Dkdr' : True,
                               'Dkca' : True,
                               'Dcan' : True,
                               'Dh' : True,
                               'Desyn' : True,
                               'Disyn' : True}
        self.ch_enable_arr = [2,2,True,True,True,True,True,True,2,2,2,True,True,True,True,True,True,True,2,2]
        
        # groupbox label
        self.uid.groupBox_Sca.setTitle('d[Ca' + super_2_p + ']' + sub_i + '/dt  ')
        self.uid.label_Snaf.setText('I' + sub_Naf)
        self.uid.label_Snap.setText('I' + sub_Nap)
        self.uid.label_Skdr.setText('I' + sub_Kdr)
        self.uid.label_Skca.setText('I' + sub_KCa)
        self.uid.label_Scan.setText('I' + sub_Can)
        self.uid.label_Sh.setText('I' + sub_H)
        self.uid.label_Sesyn.setText('I' + sub_esyn)
        self.uid.label_Sisyn.setText('I' + sub_isyn)
        self.uid.groupBox_Dca.setTitle('d[Ca' + super_2_p + ']' + sub_i + '/dt  ')
        self.uid.label_Dcal.setText('I' + sub_Cal)
        self.uid.label_Dnaf.setText('I' + sub_Naf)
        self.uid.label_Dnap.setText('I' + sub_Nap)
        self.uid.label_Dkdr.setText('I' + sub_Kdr)
        self.uid.label_Dkca.setText('I' + sub_KCa)
        self.uid.label_Dcan.setText('I' + sub_Can)
        self.uid.label_Dh.setText('I' + sub_H)
        self.uid.label_Desyn.setText('I' + sub_esyn)
        self.uid.label_Disyn.setText('I' + sub_isyn)
        
        # Parameter label
        mn_parameter_1 = QLabel('R' + sub_N + M_ohm)
        mn_parameter_2 = QLabel(tau + sub_m + s)
        mn_parameter_3 = QLabel('VA' + sub_SD + super_DC)
        mn_parameter_4 = QLabel('VA' + sub_DS + super_DC)
        mn_parameter_5 = QLabel('VA' + sub_SD + super_AC)
        mn_parameter_6 = QLabel('p')
        mn_parameter_7 = QLabel('E' + sub_Leak_S + mV)
        mn_parameter_8 = QLabel('E' + sub_Leak_D + mV)
        mu_parameter = QLabel('CV' + sub_axon + m_s)
        mn_parameter_9 = QLabel('f')
        mn_parameter_10 = QLabel('K' + sub_Ca + ms_1)
        mn_parameter_11 = QLabel(alpha + mol_uC_cm_2)
        mn_parameter_12 = QLabel('[Ca' + super_2_p + ']' + sub_o + mM)
        mn_parameter_126 = QLabel('E' + sub_Ca + mV)
        mn_parameter_13 = QLabel('G' + mS_cm_2)
        mn_parameter_14 = QLabel('E' + sub_Na + mV)
        mn_parameter_15 = QLabel(alpha + sub_1 + mV_ms_1) 
        mn_parameter_16 = QLabel(alpha + sub_2 + mV)
        mn_parameter_17 = QLabel(alpha + sub_3 + mV)
        mn_parameter_18 = QLabel(alpha + sub_4)
        mn_parameter_19 = QLabel(beta + sub_1 + mV_ms_1)
        mn_parameter_20 = QLabel(beta + sub_2 + mV)
        mn_parameter_21 = QLabel(beta + sub_3 + mV)
        mn_parameter_22 = QLabel(beta + sub_4)
        mn_parameter_23 = QLabel(gamma + sub_1 + mV)
        mn_parameter_24 = QLabel(gamma + sub_2 + mV)
        mn_parameter_25 = QLabel(gamma + sub_3 + mV)
        mn_parameter_26 = QLabel(gamma + sub_4 + mV)
        mn_parameter_27 = QLabel(gamma + sub_5 + mV)
        mn_parameter_28 = QLabel(gamma + sub_6 + mV)
        mn_parameter_29 = QLabel('G' + mS_cm_2)
        mn_parameter_30 = QLabel('E' + sub_Na + mV)
        mn_parameter_31 = QLabel(alpha + sub_1 + mV_ms_1)
        mn_parameter_32 = QLabel(alpha + sub_2 + mV)
        mn_parameter_33 = QLabel(alpha + sub_3 + mV)
        mn_parameter_34 = QLabel(alpha + sub_4)
        mn_parameter_35 = QLabel(beta + sub_1 + mV_ms_1)
        mn_parameter_36 = QLabel(beta + sub_2 + mV)
        mn_parameter_37 = QLabel(beta + sub_3 + mV)
        mn_parameter_38 = QLabel(beta + sub_4)  
        mn_parameter_39 = QLabel('G' + mS_cm_2)
        mn_parameter_40 = QLabel('E' + sub_K + mV)
        mn_parameter_41 = QLabel(gamma + sub_1 + mV)
        mn_parameter_42 = QLabel(gamma + sub_2 + mV)
        mn_parameter_43 = QLabel(gamma + sub_3 + mV)
        mn_parameter_44 = QLabel(gamma + sub_4 + mV)
        mn_parameter_45 = QLabel(gamma + sub_5 + mV)
        mn_parameter_46 = QLabel(gamma + sub_6 + ms)
        mn_parameter_47 = QLabel('G' + mS_cm_2)
        mn_parameter_48 = QLabel('E' + sub_K + mV)
        mn_parameter_49 = QLabel('K' + sub_d + mM)  
        mn_parameter_50 = QLabel('G' + mS_cm_2)
        mn_parameter_51 = QLabel(gamma + sub_1 + mV)
        mn_parameter_52 = QLabel(gamma + sub_2 + mV)
        mn_parameter_53 = QLabel(gamma + sub_3 + ms)
        mn_parameter_54 = QLabel(gamma + sub_1 + mV)
        mn_parameter_55 = QLabel(gamma + sub_2 + mV)
        mn_parameter_56 = QLabel(gamma + sub_3 + ms)
        mn_parameter_57 = QLabel('G' + mS_cm_2)
        mn_parameter_58 = QLabel('E' + mV)
        mn_parameter_59 = QLabel(gamma + sub_1 + mV)
        mn_parameter_60 = QLabel(gamma + sub_2 + mV)
        mn_parameter_61 = QLabel(gamma + sub_3 + ms)
        mn_parameter_62 = QLabel('E' + sub_esyn + mV)
        mn_parameter_64 = QLabel('E' + sub_isyn + mV)
        mn_parameter_66 = QLabel('f')
        mn_parameter_67 = QLabel('K' + sub_Ca + ms_1)
        mn_parameter_68 = QLabel(alpha + mol_uC_cm_2)
        mn_parameter_69 = QLabel('[Ca' + super_2_p + ']' + sub_o + mM)
        mn_parameter_127 = QLabel('E' + sub_Ca + mV)
        mn_parameter_70 = QLabel('G' + mS_cm_2)
        mn_parameter_71 = QLabel(gamma + sub_1 + mV)
        mn_parameter_72 = QLabel(gamma + sub_2 + mV)
        mn_parameter_73 = QLabel(gamma + sub_3 + ms)
        mn_parameter_74 = QLabel('G' + mS_cm_2)
        mn_parameter_75 = QLabel('E' + sub_Na + mV)
        mn_parameter_76 = QLabel(alpha + sub_1 + mV_ms_1) 
        mn_parameter_77 = QLabel(alpha + sub_2 + mV)
        mn_parameter_78 = QLabel(alpha + sub_3 + mV)
        mn_parameter_79 = QLabel(alpha + sub_4)
        mn_parameter_80 = QLabel(beta + sub_1 + mV_ms_1)
        mn_parameter_81 = QLabel(beta + sub_2 + mV)
        mn_parameter_82 = QLabel(beta + sub_3 + mV)
        mn_parameter_83 = QLabel(beta + sub_4)
        mn_parameter_84 = QLabel(gamma + sub_1 + mV)
        mn_parameter_85 = QLabel(gamma + sub_2 + mV)
        mn_parameter_86 = QLabel(gamma + sub_3 + mV)
        mn_parameter_87 = QLabel(gamma + sub_4 + mV)
        mn_parameter_88 = QLabel(gamma + sub_5 + mV)
        mn_parameter_89 = QLabel(gamma + sub_6 + mV)
        mn_parameter_90 = QLabel('G' + mS_cm_2)
        mn_parameter_91 = QLabel('E' + sub_Na + mV)
        mn_parameter_92 = QLabel(alpha + sub_1 + mV_ms_1)
        mn_parameter_93 = QLabel(alpha + sub_2 + mV)
        mn_parameter_94 = QLabel(alpha + sub_3 + mV)
        mn_parameter_95 = QLabel(alpha + sub_4)
        mn_parameter_96 = QLabel(beta + sub_1 + mV_ms_1)
        mn_parameter_97 = QLabel(beta + sub_2 + mV)
        mn_parameter_98 = QLabel(beta + sub_3 + mV)
        mn_parameter_99 = QLabel(beta + sub_4)
        mn_parameter_100 = QLabel('G' + mS_cm_2)
        mn_parameter_101 = QLabel('E' + sub_K + mV)
        mn_parameter_102 = QLabel(gamma + sub_1 + mV)
        mn_parameter_103 = QLabel(gamma + sub_2 + mV)
        mn_parameter_104 = QLabel(gamma + sub_3 + mV)
        mn_parameter_105 = QLabel(gamma + sub_4 + mV)
        mn_parameter_106 = QLabel(gamma + sub_5 + mV)
        mn_parameter_107 = QLabel(gamma + sub_6 + ms)
        mn_parameter_108 = QLabel('G' + mS_cm_2)
        mn_parameter_109 = QLabel('E' + sub_K + mV)
        mn_parameter_110 = QLabel('K' + sub_d + mM)
        mn_parameter_111 = QLabel('G' + mS_cm_2)
        mn_parameter_112 = QLabel(gamma + sub_1 + mV)
        mn_parameter_113 = QLabel(gamma + sub_2 + mV)
        mn_parameter_114 = QLabel(gamma + sub_3 + ms)
        mn_parameter_115 = QLabel(gamma + sub_1 + mV)
        mn_parameter_116 = QLabel(gamma + sub_2 + mV)
        mn_parameter_117 = QLabel(gamma + sub_3 + ms)
        mn_parameter_118 = QLabel('G' + mS_cm_2)
        mn_parameter_119 = QLabel('E' + mV)
        mn_parameter_120 = QLabel(gamma + sub_1 + mV)
        mn_parameter_121 = QLabel(gamma + sub_2 + mV)
        mn_parameter_122 = QLabel(gamma + sub_3 + ms)
        mn_parameter_123 = QLabel('E' + sub_esyn + mV)
        mn_parameter_125 = QLabel('E' + sub_isyn + mV)
        mf_parameter_1 = QLabel('K1' + M_1ms_1)
        mf_parameter_2 = QLabel('K2' + ms_1)
        mf_parameter_3 = QLabel('K3' + M_1ms_1)
        mf_parameter_4 = QLabel('K4' + ms_1)
        mf_parameter_5 = QLabel('K5' + M_1ms_1)
        mf_parameter_6 = QLabel('K6' + sub_i + ms_1)
        mf_parameter_7 = QLabel('K' + M_1)
        mf_parameter_8 = QLabel('R' + sub_max + ms_1)
        mf_parameter_9 = QLabel('U' + sub_max + M_ms_1)
        mf_parameter_10 = QLabel(tau + sub_1+ms)
        mf_parameter_11 = QLabel(tau + sub_2+ms)
        mf_parameter_12 = QLabel(pi+sub_1+mm_1)
        mf_parameter_13 = QLabel(pi+sub_2)
        mf_parameter_14 = QLabel(pi+sub_3+mm_1)
        mf_parameter_15 = QLabel(pi+sub_4)        
        mf_parameter_16 = QLabel('C1')
        mf_parameter_17 = QLabel('C2')
        mf_parameter_18 = QLabel('C3' + ms)
        mf_parameter_19 = QLabel('C4')
        mf_parameter_20 = QLabel('C5')
        mf_parameter_21 = QLabel(K_SE + mm_1)
        mf_parameter_22 = QLabel('P'+sub_0+n)
        mf_parameter_23 = QLabel('g'+sub_1+mm)
        mf_parameter_24 = QLabel('g'+sub_2+mm)
        mf_parameter_26 = QLabel('a'+sub_0+n)
        mf_parameter_27 = QLabel('b'+sub_0+mm_s_1)
        mf_parameter_28 = QLabel('c'+sub_0+n)
        mf_parameter_29 = QLabel('d'+sub_0+mm_s_1)

        mn_parameter = []
        mn_parameter.append(mn_parameter_1)
        mn_parameter.append(mn_parameter_2)
        mn_parameter.append(mn_parameter_3)
        mn_parameter.append(mn_parameter_4)
        mn_parameter.append(mn_parameter_5)
        mn_parameter.append(mn_parameter_6)
        mn_parameter.append(mn_parameter_7)
        mn_parameter.append(mn_parameter_8)
        mn_parameter.append(mu_parameter)
        mn_parameter.append(mn_parameter_9)
        mn_parameter.append(mn_parameter_10)
        mn_parameter.append(mn_parameter_11)
        mn_parameter.append(mn_parameter_12)
        mn_parameter.append(mn_parameter_126)
        mn_parameter.append(mn_parameter_13)
        mn_parameter.append(mn_parameter_14)
        mn_parameter.append(mn_parameter_15)
        mn_parameter.append(mn_parameter_16)
        mn_parameter.append(mn_parameter_17)
        mn_parameter.append(mn_parameter_18)
        mn_parameter.append(mn_parameter_19)
        mn_parameter.append(mn_parameter_20)      
        mn_parameter.append(mn_parameter_21)
        mn_parameter.append(mn_parameter_22)
        mn_parameter.append(mn_parameter_23)
        mn_parameter.append(mn_parameter_24)
        mn_parameter.append(mn_parameter_25)
        mn_parameter.append(mn_parameter_26)
        mn_parameter.append(mn_parameter_27)
        mn_parameter.append(mn_parameter_28)
        mn_parameter.append(mn_parameter_29)
        mn_parameter.append(mn_parameter_30)    
        mn_parameter.append(mn_parameter_31)
        mn_parameter.append(mn_parameter_32)
        mn_parameter.append(mn_parameter_33)
        mn_parameter.append(mn_parameter_34)
        mn_parameter.append(mn_parameter_35)
        mn_parameter.append(mn_parameter_36)
        mn_parameter.append(mn_parameter_37)
        mn_parameter.append(mn_parameter_38)
        mn_parameter.append(mn_parameter_39)
        mn_parameter.append(mn_parameter_40)    
        mn_parameter.append(mn_parameter_41)
        mn_parameter.append(mn_parameter_42)
        mn_parameter.append(mn_parameter_43)
        mn_parameter.append(mn_parameter_44)
        mn_parameter.append(mn_parameter_45)
        mn_parameter.append(mn_parameter_46)
        mn_parameter.append(mn_parameter_47)
        mn_parameter.append(mn_parameter_48)
        mn_parameter.append(mn_parameter_49)
        mn_parameter.append(mn_parameter_50)    
        mn_parameter.append(mn_parameter_51)
        mn_parameter.append(mn_parameter_52)
        mn_parameter.append(mn_parameter_53)
        mn_parameter.append(mn_parameter_54)
        mn_parameter.append(mn_parameter_55)
        mn_parameter.append(mn_parameter_56)
        mn_parameter.append(mn_parameter_57)
        mn_parameter.append(mn_parameter_58)
        mn_parameter.append(mn_parameter_59)
        mn_parameter.append(mn_parameter_60) 
        mn_parameter.append(mn_parameter_61)
        mn_parameter.append(mn_parameter_62)
        mn_parameter.append(mn_parameter_64)
        mn_parameter.append(mn_parameter_66)
        mn_parameter.append(mn_parameter_67)
        mn_parameter.append(mn_parameter_68)
        mn_parameter.append(mn_parameter_69)
        mn_parameter.append(mn_parameter_127)
        mn_parameter.append(mn_parameter_70)
        mn_parameter.append(mn_parameter_71)
        mn_parameter.append(mn_parameter_72)
        mn_parameter.append(mn_parameter_73)
        mn_parameter.append(mn_parameter_74)
        mn_parameter.append(mn_parameter_75)
        mn_parameter.append(mn_parameter_76)
        mn_parameter.append(mn_parameter_77)
        mn_parameter.append(mn_parameter_78)
        mn_parameter.append(mn_parameter_79)
        mn_parameter.append(mn_parameter_80)      
        mn_parameter.append(mn_parameter_81)
        mn_parameter.append(mn_parameter_82)
        mn_parameter.append(mn_parameter_83)
        mn_parameter.append(mn_parameter_84)
        mn_parameter.append(mn_parameter_85)
        mn_parameter.append(mn_parameter_86)
        mn_parameter.append(mn_parameter_87)
        mn_parameter.append(mn_parameter_88)
        mn_parameter.append(mn_parameter_89)
        mn_parameter.append(mn_parameter_90)    
        mn_parameter.append(mn_parameter_91)
        mn_parameter.append(mn_parameter_92)
        mn_parameter.append(mn_parameter_93)
        mn_parameter.append(mn_parameter_94)
        mn_parameter.append(mn_parameter_95)
        mn_parameter.append(mn_parameter_96)
        mn_parameter.append(mn_parameter_97)
        mn_parameter.append(mn_parameter_98)
        mn_parameter.append(mn_parameter_99)
        mn_parameter.append(mn_parameter_100)    
        mn_parameter.append(mn_parameter_101)
        mn_parameter.append(mn_parameter_102)
        mn_parameter.append(mn_parameter_103)
        mn_parameter.append(mn_parameter_104)
        mn_parameter.append(mn_parameter_105)
        mn_parameter.append(mn_parameter_106)
        mn_parameter.append(mn_parameter_107)
        mn_parameter.append(mn_parameter_108)
        mn_parameter.append(mn_parameter_109)
        mn_parameter.append(mn_parameter_110)    
        mn_parameter.append(mn_parameter_111)
        mn_parameter.append(mn_parameter_112)
        mn_parameter.append(mn_parameter_113)
        mn_parameter.append(mn_parameter_114)
        mn_parameter.append(mn_parameter_115)
        mn_parameter.append(mn_parameter_116)
        mn_parameter.append(mn_parameter_117)
        mn_parameter.append(mn_parameter_118)
        mn_parameter.append(mn_parameter_119)
        mn_parameter.append(mn_parameter_120)
        mn_parameter.append(mn_parameter_121)
        mn_parameter.append(mn_parameter_122)
        mn_parameter.append(mn_parameter_123)
        mn_parameter.append(mn_parameter_125)
        mf_parameter = []
        mf_parameter.append(mf_parameter_1)
        mf_parameter.append(mf_parameter_2)
        mf_parameter.append(mf_parameter_3)
        mf_parameter.append(mf_parameter_4)
        mf_parameter.append(mf_parameter_5)
        mf_parameter.append(mf_parameter_6)
        mf_parameter.append(mf_parameter_7)
        mf_parameter.append(mf_parameter_8)
        mf_parameter.append(mf_parameter_9)
        mf_parameter.append(mf_parameter_10)
        mf_parameter.append(mf_parameter_11)
        mf_parameter.append(mf_parameter_12)
        mf_parameter.append(mf_parameter_13)
        mf_parameter.append(mf_parameter_14)
        mf_parameter.append(mf_parameter_15)
        mf_parameter.append(mf_parameter_16)
        mf_parameter.append(mf_parameter_17)
        mf_parameter.append(mf_parameter_18)
        mf_parameter.append(mf_parameter_19)
        mf_parameter.append(mf_parameter_20)
        mf_parameter.append(mf_parameter_21)
        mf_parameter.append(mf_parameter_22)
        mf_parameter.append(mf_parameter_23)
        mf_parameter.append(mf_parameter_24)
        mf_parameter.append(mf_parameter_26)
        mf_parameter.append(mf_parameter_27)
        mf_parameter.append(mf_parameter_28)
        mf_parameter.append(mf_parameter_29)

        # tableWidget
        # row, col
        col_num = 3
        iter_num = len(self.mn_setValue)
        for i in range(iter_num):
            self.mn_tableWidget[i].setColumnCount(col_num)
            self.mn_tableWidget[i].setRowCount(len(self.mn_setValue[i]))
        
        iter_num = len(self.mf_setValue)
        for i in range(iter_num):
            self.mf_tableWidget[i].setColumnCount(col_num)
            self.mf_tableWidget[i].setRowCount(len(self.mf_setValue[i]))
        
        # label
        iter_num = len(mn_parameter)
        idx = self.mn_idx
        k = 0
        for i in range(iter_num):
            self.mn_tableWidget[k].setCellWidget(i-idx[k], 0, mn_parameter[i])
            if(i == idx[k+1]-1):
                k+=1
        
        iter_num = len(mf_parameter)
        idx = self.mf_idx
        k = 0
        for i in range(iter_num):
            self.mf_tableWidget[k].setCellWidget(i-idx[k], 0, mf_parameter[i])
            if(i == idx[k+1]-1):
                k+=1
        
        # deactivation cell
        iter_num = len(mn_parameter)
        idx = self.mn_idx
        k = 0
        for i in range(iter_num):
            item = QTableWidgetItem()
            item.setFlags(Qt.ItemIsEnabled)
            self.mn_tableWidget[k].setItem(i-idx[k], 2, item)
            if(i == idx[k+1]-1):
                k+=1
        
        iter_num = len(mf_parameter)
        idx = self.mf_idx
        k = 0
        for i in range(iter_num):
            item = QTableWidgetItem()
            item.setFlags(Qt.ItemIsEnabled)
            self.mf_tableWidget[k].setItem(i-idx[k], 2, item)
            if(i == idx[k+1]-1):
                k+=1
        
        # GUI operation-function      
        self.uid.mn_load_pushButton.clicked.connect(self.openFile)
        self.uid.mf_load_pushButton.clicked.connect(self.openFile)
        self.uid.buttonBox.accepted.connect(self.setValue)
        self.uid.applyButton.clicked.connect(self.setValue)
        self.uid.tableWidget_share.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Sca.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Snaf.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Snap.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Skdr.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Skca.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Scan.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Sh.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Sesyn.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Sisyn.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Dca.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Dcal.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Dnaf.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Dnap.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Dkdr.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Dkca.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Dcan.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Dh.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Desyn.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Disyn.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_m1.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_m2.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_m3.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_share.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Sca.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Snaf.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Snap.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Skdr.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Skca.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Scan.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Sh.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Sesyn.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Sisyn.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Dca.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Dcal.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Dnaf.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Dnap.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Dkdr.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Dkca.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Dcan.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Dh.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Desyn.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Disyn.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_m1.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_m2.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_m3.cellChanged.connect(self.checkValue)
        self.uid.checkBox_const_SEca.stateChanged.connect(self.checkValue)
        self.uid.checkBox_const_DEca.stateChanged.connect(self.checkValue)
        
    def keyPressEvent(self, event):
        if(event.key() == 0x01000005 or event.key() == 0x01000004):
            pass
        else:
            # QDialog event
            super(ParameterSettingWindow, self).keyPressEvent(event)

    def checkValue(self, row=3, col=3):
        if(col == 3):
            sub_Ca = '<sub>Ca</sub>'
            mv = ' (mV)'
            row = 4
            l_col = 0
            v_col = 1
            label_1 = QLabel('E' + sub_Ca + mv)
            label_2 = QLabel('E' + sub_Ca + mv)
            act_color = Qt.black
            de_color = Qt.darkGray
            
            # Eca checkBox
            if(self.uid.checkBox_const_SEca.isChecked()):
                self.const_sEca = True
                label_1.setEnabled(True)
                self.MW.SGW.setItemColorEnable(self.uid.tableWidget_Sca, row, v_col, act_color, True)
            else:
                self.const_sEca = False
                label_1.setEnabled(False)
                self.MW.SGW.setItemColorEnable(self.uid.tableWidget_Sca, row, v_col, de_color, False)
            
            if(self.uid.checkBox_const_DEca.isChecked()):
                self.const_dEca = True
                label_2.setEnabled(True)
                self.MW.SGW.setItemColorEnable(self.uid.tableWidget_Dca, row, v_col, act_color, True)
            else:
                self.const_dEca = False
                label_2.setEnabled(False)
                self.MW.SGW.setItemColorEnable(self.uid.tableWidget_Dca, row, v_col, de_color, False)
            self.uid.tableWidget_Sca.setCellWidget(row, l_col, label_1)
            self.uid.tableWidget_Dca.setCellWidget(row, l_col, label_2)
        
        elif(col == 1):
            tableWidget = self.sender()
            
            for i in range(len(self.mn_tableWidget)):
                if(tableWidget == self.mn_tableWidget[i]):
                    item = tableWidget.item(row, col)
                    defValue = self.mn_defaultValue[i][row]
                    setValue = self.mn_setValue[i][row]
                    
                    # opossite reversal potential
                    if(self.init == True and row == 1):
                        if(tableWidget == self.uid.tableWidget_Snaf):
                            op_tablewidget = self.uid.tableWidget_Snap
                        elif(tableWidget == self.uid.tableWidget_Snap):
                            op_tablewidget = self.uid.tableWidget_Snaf
                        elif(tableWidget == self.uid.tableWidget_Skdr):
                            op_tablewidget = self.uid.tableWidget_Skca
                        elif(tableWidget == self.uid.tableWidget_Skca):
                            op_tablewidget = self.uid.tableWidget_Skdr
                        elif(tableWidget == self.uid.tableWidget_Dnaf):
                            op_tablewidget = self.uid.tableWidget_Dnap
                        elif(tableWidget == self.uid.tableWidget_Dnap):
                            op_tablewidget = self.uid.tableWidget_Dnaf
                        elif(tableWidget == self.uid.tableWidget_Dkdr):
                            op_tablewidget = self.uid.tableWidget_Dkca
                        elif(tableWidget == self.uid.tableWidget_Dkca):
                            op_tablewidget = self.uid.tableWidget_Dkdr
                        else:
                            continue
                        thisValue = float(item.text())
                        op_item = op_tablewidget.item(1, 1)
                        opValue = float(op_item.text())
                        if(round(thisValue, 10) != round(opValue, 10)):
                            temp_item = QTableWidgetItem(str(thisValue))
                            op_tablewidget.setItem(1, 1, temp_item)
            
            for i in range(len(self.mf_tableWidget)):
                if(tableWidget == self.mf_tableWidget[i]):
                    item = tableWidget.item(row, col)
                    defValue = self.mf_defaultValue[i][row]
                    setValue = self.mf_setValue[i][row]
                    
                    if(tableWidget == self.uid.tableWidget_m3):
                        # a0, c0 are proportional to P0 rate of change when P0 is changed.
                        if(self.init == True and row == 1):
                            p0 = float(item.text())
                            a0_item = tableWidget.item(4, col)
                            a0 = float(a0_item.text())
                            c0_item = tableWidget.item(6, col)
                            c0 = float(c0_item.text())
                            def_p0 = defValue
                            dp0 = p0/def_p0
                            def_a0 = self.mf_defaultValue[i][4]
                            def_c0 = self.mf_defaultValue[i][6]
                            a0 = def_a0*dp0
                            c0 = def_c0*dp0
                            temp_item = QTableWidgetItem(str(a0))
                            tableWidget.setItem(4, 1, temp_item)
                            temp_item = QTableWidgetItem(str(c0))
                            tableWidget.setItem(6, 1, temp_item)

            # compare with default value
            try:
                if(round(float(item.text()), 10) != round(defValue, 10)):
                    item = QTableWidgetItem()
                    icon = QIcon()
                    icon.addPixmap(QPixmap("./resources/default.png"))
                    item.setIcon(icon)
                    item.setFlags(Qt.ItemIsEnabled)
                    tableWidget.setItem(row, 2, item)
                else:
                    item = QTableWidgetItem()
                    item.setFlags(Qt.ItemIsEnabled)  
                    tableWidget.setItem(row, 2, item) 
            except:
                item = QTableWidgetItem(str(setValue))
                tableWidget.setItem(row, 1, item)
    
    def returnDefaultValue(self, row, col):
        if(col == 2):  
            tableWidget = self.sender()
            
            for i in range(len(self.mn_tableWidget)):
                if(tableWidget == self.mn_tableWidget[i]):
                    # prevent return value of Eca that is deactivated.
                    if(tableWidget==self.uid.tableWidget_Sca):
                        if(row == 4):
                            if(not self.uid.checkBox_const_SEca.isChecked()):
                                continue
                    if(tableWidget==self.uid.tableWidget_Dca):
                        if(row == 4):
                            if(not self.uid.checkBox_const_DEca.isChecked()):
                                continue
                    
                    # return the default value
                    item = QTableWidgetItem(str(self.mn_defaultValue[i][row]))
                    self.mn_tableWidget[i].setItem(row, 1, item)
                    
            for i in range(len(self.mf_tableWidget)):
                if(tableWidget == self.mf_tableWidget[i]):
                    item = QTableWidgetItem(str(self.mf_defaultValue[i][row]))
                    self.mf_tableWidget[i].setItem(row, 1, item)
                
    def saveChEnable(self):
        # save the groupbox enable information
        if(self.uid.groupBox_Snaf.isChecked()):
            self.channel_enable['Snaf'] = True
            self.ch_enable_arr[2] = True
        else:
            self.channel_enable['Snaf'] = False
            self.ch_enable_arr[2] = False
        if(self.uid.groupBox_Snap.isChecked()):
            self.channel_enable['Snap'] = True
            self.ch_enable_arr[3] = True
        else:
            self.channel_enable['Snap'] = False
            self.ch_enable_arr[3] = False
        if(self.uid.groupBox_Skdr.isChecked()):
            self.channel_enable['Skdr'] = True
            self.ch_enable_arr[4] = True
        else:
            self.channel_enable['Skdr'] = False
            self.ch_enable_arr[4] = False
        if(self.uid.groupBox_Skca.isChecked()):
            self.channel_enable['Skca'] = True
            self.ch_enable_arr[5] = True
        else:
            self.channel_enable['Skca'] = False
            self.ch_enable_arr[5] = False
        if(self.uid.groupBox_Scan.isChecked()):
            self.channel_enable['Scan'] = True
            self.ch_enable_arr[6] = True
        else:
            self.channel_enable['Scan'] = False
            self.ch_enable_arr[6] = False
        if(self.uid.groupBox_Sh.isChecked()):
            self.channel_enable['Sh'] = True
            self.ch_enable_arr[7] = True
        else:
            self.channel_enable['Sh'] = False
            self.ch_enable_arr[7] = False
        if(self.uid.groupBox_Sesyn.isChecked()):
            self.channel_enable['Sesyn'] = True
        else:
            self.channel_enable['Sesyn'] = False
        if(self.uid.groupBox_Sisyn.isChecked()):
            self.channel_enable['Sisyn'] = True
        else:
            self.channel_enable['Sisyn'] = False
        if(self.uid.groupBox_Dcal.isChecked()):
            self.channel_enable['Dcal'] = True
            self.ch_enable_arr[11] = True
        else:
            self.channel_enable['Dcal'] = False
            self.ch_enable_arr[11] = False
        if(self.uid.groupBox_Dnaf.isChecked()):
            self.channel_enable['Dnaf'] = True
            self.ch_enable_arr[12] = True
        else:
            self.channel_enable['Dnaf'] = False
            self.ch_enable_arr[12] = False
        if(self.uid.groupBox_Dnap.isChecked()):
            self.channel_enable['Dnap'] = True
            self.ch_enable_arr[13] = True
        else:
            self.channel_enable['Dnap'] = False
            self.ch_enable_arr[13] = False
        if(self.uid.groupBox_Dkdr.isChecked()):
            self.channel_enable['Dkdr'] = True
            self.ch_enable_arr[14] = True
        else:
            self.channel_enable['Dkdr'] = False
            self.ch_enable_arr[14] = False
        if(self.uid.groupBox_Dkca.isChecked()):
            self.channel_enable['Dkca'] = True
            self.ch_enable_arr[15] = True
        else:
            self.channel_enable['Dkca'] = False
            self.ch_enable_arr[15] = False
        if(self.uid.groupBox_Dcan.isChecked()):
            self.channel_enable['Dcan'] = True
            self.ch_enable_arr[16] = True
        else:
            self.channel_enable['Dcan'] = False
            self.ch_enable_arr[16] = False
        if(self.uid.groupBox_Dh.isChecked()):
            self.channel_enable['Dh'] = True
            self.ch_enable_arr[17] = True
        else:
            self.channel_enable['Dh'] = False
            self.ch_enable_arr[17] = False
        if(self.uid.groupBox_Desyn.isChecked()):
            self.channel_enable['Desyn'] = True
        else:
            self.channel_enable['Desyn'] = False
        if(self.uid.groupBox_Disyn.isChecked()):
            self.channel_enable['Disyn'] = True
        else:
            self.channel_enable['Disyn'] = False
            
    def checkChEnable(self):
        # set the groupbox enable based on groupbox enable information
        if(self.channel_enable['Snaf'] == True):
            self.uid.groupBox_Snaf.setChecked(True)
        else:
            self.uid.groupBox_Snaf.setChecked(False)
        if(self.channel_enable['Snap'] == True):
            self.uid.groupBox_Snap.setChecked(True)
        else:
            self.uid.groupBox_Snap.setChecked(False)
        if(self.channel_enable['Skdr'] == True):
            self.uid.groupBox_Skdr.setChecked(True)
        else:
            self.uid.groupBox_Skdr.setChecked(False)
        if(self.channel_enable['Skca'] == True):
            self.uid.groupBox_Skca.setChecked(True)
        else:
            self.uid.groupBox_Skca.setChecked(False)
        if(self.channel_enable['Scan'] == True):
            self.uid.groupBox_Scan.setChecked(True)
        else:
            self.uid.groupBox_Scan.setChecked(False)
        if(self.channel_enable['Sh'] == True):
            self.uid.groupBox_Sh.setChecked(True)
        else:
            self.uid.groupBox_Sh.setChecked(False)
        if(self.channel_enable['Sesyn'] == True):
            self.uid.groupBox_Sesyn.setChecked(True)
        else:
            self.uid.groupBox_Sesyn.setChecked(False)
        if(self.channel_enable['Sisyn'] == True):
            self.uid.groupBox_Sisyn.setChecked(True)
        else:
            self.uid.groupBox_Sisyn.setChecked(False)
        if(self.channel_enable['Dcal'] == True):
            self.uid.groupBox_Dcal.setChecked(True)
        else:
            self.uid.groupBox_Dcal.setChecked(False)
        if(self.channel_enable['Dnaf'] == True):
            self.uid.groupBox_Dnaf.setChecked(True)
        else:
            self.uid.groupBox_Dnaf.setChecked(False)
        if(self.channel_enable['Dnap'] == True):
            self.uid.groupBox_Dnap.setChecked(True)
        else:
            self.uid.groupBox_Dnap.setChecked(False)
        if(self.channel_enable['Dkdr'] == True):
            self.uid.groupBox_Dkdr.setChecked(True)
        else:
            self.uid.groupBox_Dkdr.setChecked(False)
        if(self.channel_enable['Dkca'] == True):
            self.uid.groupBox_Dkca.setChecked(True)
        else:
            self.uid.groupBox_Dkca.setChecked(False)
        if(self.channel_enable['Dcan'] == True):
            self.uid.groupBox_Dcan.setChecked(True)
        else:
            self.uid.groupBox_Dcan.setChecked(False)
        if(self.channel_enable['Dh'] == True):
            self.uid.groupBox_Dh.setChecked(True)
        else:
            self.uid.groupBox_Dh.setChecked(False)
        if(self.channel_enable['Desyn'] == True):
            self.uid.groupBox_Desyn.setChecked(True)
        else:
            self.uid.groupBox_Desyn.setChecked(False)
        if(self.channel_enable['Disyn'] == True):
            self.uid.groupBox_Disyn.setChecked(True)
        else:
            self.uid.groupBox_Disyn.setChecked(False)

    def displayValue(self, value=[], index=2):
        # display set value
        if(self.MW.ModelType == 'Motoneuron' or self.MW.ModelType == 'Motor Unit'):  
            self.checkChEnable()
            if(index != 0):
                param=self.mn_setValue
            else:
                param=value
            
            iter_num = self.mn_setValue_num
            idx = self.mn_idx
            k = 0            
            for i in range(iter_num):
                if(i == idx[k] and self.ch_enable_arr[k] == False):
                    item = QTableWidgetItem(str(self.pre_G[k]))
                else:
                    item = QTableWidgetItem(str(param[k][i-idx[k]]))
                self.mn_tableWidget[k].setItem(i-idx[k], 1, item)
                if(i == idx[k+1]-1):
                    k+=1    
        
        if(self.MW.ModelType == 'Muscle Fibers' or self.MW.ModelType == 'Motor Unit'):
            if(index != 1):
                param=self.mf_setValue
            else:
                param=value            

            iter_num = self.mf_setValue_num
            idx = self.mf_idx
            k = 0            
            for i in range(iter_num):
                item = QTableWidgetItem(str(param[k][i-idx[k]]))
                self.mf_tableWidget[k].setItem(i-idx[k], 1, item)
                if(i == idx[k+1]-1):
                    k+=1
                    
        self.init = True
        
        if(self.MW.ModelType == 'Motoneuron'):
            self.uid.tabWidget.setTabEnabled(0, True)
            self.uid.tabWidget.setTabEnabled(1, False)
            self.checkValue()
            
        elif(self.MW.ModelType == 'Muscle Fibers'):
            self.uid.tabWidget.setTabEnabled(0, False)
            self.uid.tabWidget.setTabEnabled(1, True)
            
        elif(self.MW.ModelType == 'Motor Unit'):
            self.uid.tabWidget.setTabEnabled(0, True)
            self.uid.tabWidget.setTabEnabled(1, True)
            self.uid.tabWidget.setCurrentIndex(0) 
            self.checkValue()

    def openFile(self):
        curPath = os.path.dirname( os.path.abspath( sys.argv[0] ) )
        currentWidget = self.uid.tabWidget.currentWidget()        
        
        try:
            if(currentWidget == self.uid.mn_tab):
                model = 'Motoneuron'                
            elif(currentWidget == self.uid.mf_tab):
                model = 'Muscle Fibers'
            fileName = QFileDialog.getOpenFileName(self,"Open a " + model + " Parameter File (Must be a csv format)", curPath+"/parameters", "CSV Files (*.csv)")
        except:
            self.MW.setTextEdit("[Error] Failed to open file.")
            self.MW.raise_()
            self.MW.activateWindow()
        
        if(fileName != ''):
            self.importData(fileName, currentWidget)
    
    def importData(self, fileName, currentWidget):
        try:
            param = pd.read_csv(unicode(fileName))
            tempValue = range(len(param.index))
            for i in range(len(param.index)):
                tempValue[i] = float(param.Value[i])
            
            if(currentWidget == self.uid.mn_tab):
                if(len(tempValue) != self.mn_setValue_num):
                    raise Exception
                iter_num = len(self.mn_idx)
                idx = self.mn_idx
                lstSlice_arr = range(len(self.mn_idx)-1)
            elif(currentWidget == self.uid.mf_tab):
                if(len(tempValue) != self.mf_setValue_num):
                    raise Exception
                iter_num = len(self.mf_idx)
                idx = self.mf_idx
                lstSlice_arr = range(len(self.mf_idx)-1)
            for i in range(1, iter_num):
                lstSlice_arr[i-1] = tempValue[idx[i-1]:idx[i]] 
                if(math.isnan(tempValue[i])):
                    raise Exception  
            dispValue = lstSlice_arr
        
        except:
            self.MW.setTextEdit("[Error] Failed to load data.")
            self.MW.raise_()
            self.MW.activateWindow()
            return
        
        if(currentWidget == self.uid.mn_tab):
            i = 0
            self.uid.mn_lineEdit.setText(fileName)
        elif(currentWidget == self.uid.mf_tab):
            i = 1
            self.uid.mf_lineEdit.setText(fileName)
        
        self.displayValue(dispValue, i)
        self.uid.tabWidget.setCurrentIndex(i)

    def setValue(self):
        # set the value
        if(self.MW.ModelType == 'Motoneuron' or self.MW.ModelType == 'Motor Unit'):
            self.saveChEnable()
            
            iter_num = self.mn_setValue_num
            idx = self.mn_idx 
            k = 0            
            for i in range(iter_num):
                item = self.mn_tableWidget[k].item(i-idx[k], 1)
                if(i == idx[k] and self.ch_enable_arr[k] == False):
                    # save the previous G value and set the G to zero. 
                    self.pre_G[k] = float(item.text())
                    self.mn_setValue[k][i-idx[k]] = 0.
                else:
                    self.mn_setValue[k][i-idx[k]] = float(item.text()) 
                if(i == idx[k+1]-1):
                    k+=1
                    
            # set the synapic conductance signal based on channel enable information.
            if(self.channel_enable['Sesyn'] == True):
                self.MW.SGW.se_ch = True
            else:
                self.MW.SGW.se_ch = False
            if(self.channel_enable['Sisyn'] == True):
                self.MW.SGW.si_ch = True
            else:
                self.MW.SGW.si_ch = False
            if(self.channel_enable['Desyn'] == True):
                self.MW.SGW.de_ch = True
            else:
                self.MW.SGW.de_ch = False
            if(self.channel_enable['Disyn'] == True):
                self.MW.SGW.di_ch = True
            else:
                self.MW.SGW.di_ch = False
            self.MW.MN.setSynConSignal(self.MW.SGW.se_ch, self.MW.SGW.si_ch, self.MW.SGW.de_ch, self.MW.SGW.di_ch, self.MW.SGW.s_SCSG.e_times, self.MW.SGW.s_SCSG.i_times, self.MW.SGW.d_SCSG.e_times, self.MW.SGW.d_SCSG.i_times, self.MW.SGW.s_SCSG.G_e, self.MW.SGW.s_SCSG.G_i, self.MW.SGW.d_SCSG.G_e, self.MW.SGW.d_SCSG.G_i)

        if(self.MW.ModelType == 'Muscle Fibers' or self.MW.ModelType == 'Motor Unit'):
            iter_num = self.mf_setValue_num
            idx = self.mf_idx
            k = 0            
            for i in range(iter_num):
                item = self.mf_tableWidget[k].item(i-idx[k], 1)
                self.mf_setValue[k][i-idx[k]] = float(item.text())
                if(i == idx[k+1]-1):
                    k+=1

        # Text message
        if(self.MW.ModelType == 'Motoneuron'):
            self.MW.MN.setModelParam(self.mn_setValue, self.const_sEca, self.const_dEca)
            self.MW.setTextEdit("[Motoneuron] Parameter values set.")
        elif(self.MW.ModelType == 'Muscle Fibers'):
            self.MW.MF.setModelParam(self.mf_setValue)
            self.MW.setTextEdit("[Muscle Fibers] Parameter values set.")
        elif(self.MW.ModelType == 'Motor Unit'):
            self.MW.MU.setModelParam(self.mn_setValue, self.mf_setValue, self.const_sEca, self.const_dEca)
            self.MW.setTextEdit("[Motor Unit] Parameter values set.")
            
        # update pop-up window enable 
        self.MW.ISW.setObjectEnabled()
        self.MW.SGW.setObjectEnabled()
        self.MW.OW.setObjectEnabled(False)
            
# Simulation Condition Setting Window class           
class IntegrationSettingWindow(QDialog):
    def __init__(self, MW):
        super(IntegrationSettingWindow, self).__init__()
        self.MW = MW
        
        # subscript, superscript, Greek letter with unicode and unit label
        super_2 = u'\u00B2' 
        super_2_p = super_2+u'\u207A'
        sub_i = u'\u1D62'
        mV = ' (mV)'
        mM = ' (mM)'
        tilde = u'\u0303'
        mm = ' (mm)'
        M = ' (M)'
        
        # subscript with HTML
        sub_Naf = '<sub>Naf</sub>'
        sub_Nap = '<sub>Nap</sub>'
        sub_Kdr = '<sub>Kdr</sub>'
        sub_Can = '<sub>Can</sub>'
        sub_H = '<sub>H</sub>'
        sub_Cal = '<sub>Cal</sub>'
        sub_SR = '<sub>SR</sub>'
        sub_SP = '<sub>SP</sub>'
        sub_CE = '<sub>CE</sub>'
        
        ## initialization
        # set value
        self.t_setTable = [0., 10000., 0.1, 100.]
        self.t_start, self.t_stop, self.t_dt, self.t_pt = self.t_setTable
        self.mn_setValue = [[-70.], 
                            [0.0001], 
                            [0.001, 0.5829], 
                            [0.001],
                            [0.1239],
                            [0.004199, 0.9219],
                            [0.], 
                            [-70.],
                            [0.0001],
                            [0.001],
                            [0.001, 0.5829], 
                            [0.001], 
                            [0.1239], 
                            [0.004199, 0.9219], 
                            [0.]]
        self.mn_idx = [0, 1, 2, 4, 5, 6, 8, 9, 10, 11, 12, 14, 15, 16, 18, 19]
        self.mn_setValue_num = 19
        
        self.mf_setValue = [[0.03, 0.0025, 0., 0.00043, 1.e-7, 0.00007, 0., 0.],
                            [0.],
                            [-8.]] 
        self.mf_idx = [0, 8, 9, 10]     
        self.mf_setValue_num = 10
        
        # default value
        self.t_defaultTable = tuple(self.t_setTable)
        self.def_time_num = 4 
        self.mn_defaultValue = copy.deepcopy(self.mn_setValue) 
        self.mf_defaultValue = copy.deepcopy(self.mf_setValue)
        
        ## GUI
        self.uid = Ui_ID()
        self.uid.setupUi(self)
        self.mn_tableWidget = [self.uid.tableWidget_Sv,
                               self.uid.tableWidget_Sca,
                               self.uid.tableWidget_Snaf,
                               self.uid.tableWidget_Snap,
                               self.uid.tableWidget_Skdr,
                               self.uid.tableWidget_Scan,
                               self.uid.tableWidget_Sh,
                               self.uid.tableWidget_Dv,
                               self.uid.tableWidget_Dca,
                               self.uid.tableWidget_Dcal,
                               self.uid.tableWidget_Dnaf,
                               self.uid.tableWidget_Dnap,
                               self.uid.tableWidget_Dkdr,
                               self.uid.tableWidget_Dcan,
                               self.uid.tableWidget_Dh]
        self.mf_tableWidget = [self.uid.tableWidget_m1, 
                               self.uid.tableWidget_m2, 
                               self.uid.tableWidget_m3]
        
        # groupbox label
        self.uid.groupBox_Sca.setTitle('d[Ca' + super_2_p + ']' + sub_i + '/dt  ')
        self.uid.label_Snaf.setText('I' + sub_Naf)
        self.uid.label_Snap.setText('I' + sub_Nap)
        self.uid.label_Skdr.setText('I' + sub_Kdr)
        self.uid.label_Scan.setText('I' + sub_Can)
        self.uid.label_Sh.setText('I' + sub_H)
        self.uid.groupBox_Dca.setTitle('d[Ca' + super_2_p + ']' + sub_i + '/dt  ')
        self.uid.label_Dcal.setText('I' + sub_Cal)
        self.uid.label_Dnaf.setText('I' + sub_Naf)
        self.uid.label_Dnap.setText('I' + sub_Nap)
        self.uid.label_Dkdr.setText('I' + sub_Kdr)
        self.uid.label_Dcan.setText('I' + sub_Can)
        self.uid.label_Dh.setText('I' + sub_H)
        
        # lineEdit
        self.lineEdit = []
        self.lineEdit.append(self.uid.lineEdit_2)
        self.lineEdit.append(self.uid.lineEdit_3)
        self.lineEdit.append(self.uid.lineEdit_4)
        
        # checkBox
        self.checkBox = []
        self.checkBox.append(self.uid.checkBox_2)
        self.checkBox.append(self.uid.checkBox_3)
        self.checkBox.append(self.uid.checkBox_4)
        self.checkBox[0].setEnabled(False)
        self.checkBox[1].setEnabled(False)
        self.checkBox[2].setEnabled(False)
        
        # ivalue label
        mn_ivalue_1 = QLabel('V' + mV)
        mn_ivalue_2 = QLabel('[Ca' + super_2_p + ']' + sub_i + mM)
        mn_ivalue_3 = QLabel('m')
        mn_ivalue_4 = QLabel('h')
        mn_ivalue_5 = QLabel('m') 
        mn_ivalue_6 = QLabel('n') 
        mn_ivalue_7 = QLabel('m')
        mn_ivalue_8 = QLabel('h')
        mn_ivalue_9 = QLabel('m')
        mn_ivalue_12 = QLabel('V' + mV)
        mn_ivalue_13 = QLabel('[Ca' + super_2_p + ']' + sub_i + mM)
        mn_ivalue_14 = QLabel(' l')
        mn_ivalue_15 = QLabel('m')
        mn_ivalue_16 = QLabel('h')
        mn_ivalue_17 = QLabel('m')
        mn_ivalue_18 = QLabel('n')
        mn_ivalue_19 = QLabel('m')
        mn_ivalue_20 = QLabel('h')
        mn_ivalue_21 = QLabel('m')
        mf_ivalue_1 = QLabel('[CS]' + M)
        mf_ivalue_2 = QLabel('[Ca' + sub_SR + ']' + M)
        mf_ivalue_3 = QLabel('[Ca' + sub_SR + 'CS]' + M)
        mf_ivalue_4 = QLabel('[B]' + M)
        mf_ivalue_5 = QLabel('[Ca' + sub_SP + ']' + M)
        mf_ivalue_6 = QLabel('[T]' + M)
        mf_ivalue_7 = QLabel('[Ca' + sub_SP + 'B]' + M)
        mf_ivalue_8 = QLabel('[Ca' + sub_SP + 'T]' + M)
        mf_ivalue_9 = QLabel('A' + tilde)
        mf_ivalue_10 = QLabel('X' + sub_CE + mm)
        
        mn_ivalue = []
        mn_ivalue.append(mn_ivalue_1)
        mn_ivalue.append(mn_ivalue_2)
        mn_ivalue.append(mn_ivalue_3)
        mn_ivalue.append(mn_ivalue_4)
        mn_ivalue.append(mn_ivalue_5)
        mn_ivalue.append(mn_ivalue_6)
        mn_ivalue.append(mn_ivalue_7)
        mn_ivalue.append(mn_ivalue_8)
        mn_ivalue.append(mn_ivalue_9)
        mn_ivalue.append(mn_ivalue_12)
        mn_ivalue.append(mn_ivalue_13)
        mn_ivalue.append(mn_ivalue_14)
        mn_ivalue.append(mn_ivalue_15)
        mn_ivalue.append(mn_ivalue_16)
        mn_ivalue.append(mn_ivalue_17)
        mn_ivalue.append(mn_ivalue_18)
        mn_ivalue.append(mn_ivalue_19)
        mn_ivalue.append(mn_ivalue_20)      
        mn_ivalue.append(mn_ivalue_21)
        mf_ivalue = []
        mf_ivalue.append(mf_ivalue_1)
        mf_ivalue.append(mf_ivalue_2)
        mf_ivalue.append(mf_ivalue_3)
        mf_ivalue.append(mf_ivalue_4)
        mf_ivalue.append(mf_ivalue_5)
        mf_ivalue.append(mf_ivalue_6)
        mf_ivalue.append(mf_ivalue_7)
        mf_ivalue.append(mf_ivalue_8)
        mf_ivalue.append(mf_ivalue_9)
        mf_ivalue.append(mf_ivalue_10)

        # tableWidget
        # row, column
        col_num = 3
        iter_num = len(self.mn_setValue)
        for i in range(iter_num):
            self.mn_tableWidget[i].setColumnCount(col_num)
            self.mn_tableWidget[i].setRowCount(len(self.mn_setValue[i]))
        
        iter_num = len(self.mf_setValue)
        for i in range(iter_num):
            self.mf_tableWidget[i].setColumnCount(col_num)
            self.mf_tableWidget[i].setRowCount(len(self.mf_setValue[i]))
        
        # label
        iter_num = len(mn_ivalue)
        idx = self.mn_idx 
        k = 0
        for i in range(iter_num):
            self.mn_tableWidget[k].setCellWidget(i-idx[k], 0, mn_ivalue[i])
            if(i == idx[k+1]-1):
                k+=1
        
        iter_num = len(mf_ivalue)
        idx = self.mf_idx
        k = 0
        for i in range(iter_num):
            self.mf_tableWidget[k].setCellWidget(i-idx[k], 0, mf_ivalue[i])
            if(i == idx[k+1]-1):
                k+=1

        # deactivation cell
        iter_num = len(mn_ivalue)
        idx = self.mn_idx 
        k = 0
        for i in range(iter_num):
            item = QTableWidgetItem()
            item.setFlags(Qt.ItemIsEnabled)
            self.mn_tableWidget[k].setItem(i-idx[k], 2, item)
            if(i == idx[k+1]-1):
                k+=1
        
        iter_num = len(mf_ivalue)
        idx = self.mf_idx
        k = 0
        for i in range(iter_num):
            item = QTableWidgetItem()
            item.setFlags(Qt.ItemIsEnabled)
            self.mf_tableWidget[k].setItem(i-idx[k], 2, item)
            if(i == idx[k+1]-1):
                k+=1
            
        # GUI operation-function 
        self.uid.buttonBox.connect(self.uid.buttonBox, SIGNAL("accepted()"), self.setValue)
        self.uid.buttonBox.connect(self.uid.applyButton, SIGNAL("clicked()"), self.setValue)
        self.uid.lineEdit_2.editingFinished.connect(self.checkValue)
        self.uid.lineEdit_3.editingFinished.connect(self.checkValue)
        self.uid.lineEdit_4.editingFinished.connect(self.checkValue)
        self.uid.checkBox_2.stateChanged.connect(self.returnDefaultValue)
        self.uid.checkBox_3.stateChanged.connect(self.returnDefaultValue)
        self.uid.checkBox_4.stateChanged.connect(self.returnDefaultValue)
        self.uid.tableWidget_Sv.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Sca.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Snaf.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Snap.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Skdr.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Scan.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Sh.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Dv.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Dca.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Dcal.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Dnaf.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Dnap.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Dkdr.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Dcan.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Dh.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_m1.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_m2.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_m3.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_Sv.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Sca.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Snaf.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Snap.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Skdr.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Scan.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Sh.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Dv.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Dca.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Dcal.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Dnaf.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Dnap.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Dkdr.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Dcan.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_Dh.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_m1.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_m2.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_m3.cellChanged.connect(self.checkValue)
        
    def keyPressEvent(self, event):
        # operation of Key_Return or Key_Enter
        if(event.key() == 0x01000005 or event.key() == 0x01000004): 
            for i in range(len(self.lineEdit)):            
                self.lineEdit[i].clearFocus()
        else:
            # QDialog event
            super(IntegrationSettingWindow, self).keyPressEvent(event)

    def checkValue(self, row=1, col=3):
        # time tab value change
        if(col == 3):                    
            for i in range(self.def_time_num-1):
                try:                       
                    # set default value if the value is '0' or negative number                        
                    if(np.sign(float(self.lineEdit[i].text())) != 1):
                        self.lineEdit[i].setText(str(self.t_defaultTable[i+1]))
                        
                    # compare with default value
                    if(round(float(self.lineEdit[i].text()), 10) != round(self.t_defaultTable[i+1], 10)):   
                        self.checkBox[i].setCheckState(Qt.Checked)
                        self.checkBox[i].setEnabled(True)
                    else:
                        self.checkBox[i].setCheckState(Qt.Unchecked)
                        self.checkBox[i].setEnabled(False)
                except:
                    self.lineEdit[i].setText(str(self.t_setTable[i+1]))
        
        # tableWidget value change
        elif(col == 1):
            tableWidget = self.sender()
            for i in range(len(self.mn_tableWidget)):
                if(tableWidget == self.mn_tableWidget[i]):
                    item = tableWidget.item(row, col)
                    defValue = self.mn_defaultValue[i][row]
                    setValue = self.mn_setValue[i][row]

            for i in range(len(self.mf_tableWidget)):
                if(tableWidget == self.mf_tableWidget[i]):
                    item = tableWidget.item(row, col)
                    defValue = self.mf_defaultValue[i][row]
                    setValue = self.mf_setValue[i][row]
            
            try:
                # compare with default value
                if(round(float(item.text()), 10) != round(defValue, 10)):
                    item = QTableWidgetItem()
                    icon = QIcon()
                    icon.addPixmap(QPixmap("./resources/default.png"))
                    item.setIcon(icon)
                    item.setFlags(Qt.ItemIsEnabled)
                    tableWidget.setItem(row, 2, item)
                else:
                    item = QTableWidgetItem()
                    item.setFlags(Qt.ItemIsEnabled)  
                    tableWidget.setItem(row, 2, item) 
            except:
                item = QTableWidgetItem(str(setValue))
                tableWidget.setItem(row, 1, item)
                
    def returnDefaultValue(self, row, col=3):
        # time tab
        if(col == 3):
            checkBox = self.sender()
            
            # find checkBox index
            if(not checkBox.isChecked()):
                for i in range(self.def_time_num-1):
                        if(self.checkBox[i] == checkBox):
                            break
                
                # return the default value
                self.lineEdit[i].setText(str(self.t_defaultTable[i+1]))
                self.checkBox[i].setEnabled(False)
                
                # prevent widget focus
                if((i+1) < self.def_time_num-1): 
                    self.lineEdit[i+1].clearFocus()
        
        # tableWidget
        elif(col == 2): 
            tableWidget = self.sender()
            
            # return the default value
            for i in range(len(self.mn_tableWidget)):
                if(tableWidget == self.mn_tableWidget[i]):
                    item = QTableWidgetItem(str(self.mn_defaultValue[i][row]))
                    self.mn_tableWidget[i].setItem(row, 1, item)
            
            for i in range(len(self.mf_tableWidget)):
                if(tableWidget == self.mf_tableWidget[i]):
                    item = QTableWidgetItem(str(self.mf_defaultValue[i][row]))
                    self.mf_tableWidget[i].setItem(row, 1, item)
            
    def setObjectEnabled(self):
        # object enable setting
        if(self.MW.PSW.channel_enable['Snaf'] == True):
            self.uid.groupBox_Snaf.setEnabled(True)
        else:
            self.uid.groupBox_Snaf.setEnabled(False)
        if(self.MW.PSW.channel_enable['Snap'] == True):
            self.uid.groupBox_Snap.setEnabled(True)
        else:
            self.uid.groupBox_Snap.setEnabled(False)
        if(self.MW.PSW.channel_enable['Skdr'] == True):
            self.uid.groupBox_Skdr.setEnabled(True)
        else:
            self.uid.groupBox_Skdr.setEnabled(False)
        if(self.MW.PSW.channel_enable['Scan'] == True):
            self.uid.groupBox_Scan.setEnabled(True)
        else:
            self.uid.groupBox_Scan.setEnabled(False)
        if(self.MW.PSW.channel_enable['Sh'] == True):
            self.uid.groupBox_Sh.setEnabled(True)
        else:
            self.uid.groupBox_Sh.setEnabled(False)
        if(self.MW.PSW.channel_enable['Dcal'] == True):
            self.uid.groupBox_Dcal.setEnabled(True)
        else:
            self.uid.groupBox_Dcal.setEnabled(False)
        if(self.MW.PSW.channel_enable['Dnaf'] == True):
            self.uid.groupBox_Dnaf.setEnabled(True)
        else:
            self.uid.groupBox_Dnaf.setEnabled(False)
        if(self.MW.PSW.channel_enable['Dnap'] == True):
            self.uid.groupBox_Dnap.setEnabled(True)
        else:
            self.uid.groupBox_Dnap.setEnabled(False)
        if(self.MW.PSW.channel_enable['Dkdr'] == True):
            self.uid.groupBox_Dkdr.setEnabled(True)
        else:
            self.uid.groupBox_Dkdr.setEnabled(False)
        if(self.MW.PSW.channel_enable['Dcan'] == True):
            self.uid.groupBox_Dcan.setEnabled(True)
        else:
            self.uid.groupBox_Dcan.setEnabled(False)
        if(self.MW.PSW.channel_enable['Dh'] == True):
            self.uid.groupBox_Dh.setEnabled(True)
        else:
            self.uid.groupBox_Dh.setEnabled(False)
        
    def displayValue(self):
        ## display set value
        # time
        self.uid.lineEdit_2.setText(str(self.t_setTable[1]))
        self.uid.lineEdit_3.setText(str(self.t_setTable[2]))
        self.uid.lineEdit_4.setText(str(self.t_setTable[3]))
        
        # tableWidget
        if(self.MW.ModelType == 'Motoneuron' or self.MW.ModelType == 'Motor Unit'):
            self.setObjectEnabled()
            iter_num = self.mn_setValue_num
            idx = self.mn_idx
            k = 0            
            for i in range(iter_num):
                item = QTableWidgetItem(str(self.mn_setValue[k][i-idx[k]]))
                self.mn_tableWidget[k].setItem(i-idx[k], 1, item)
                if(i == idx[k+1]-1):
                    k+=1
        
        if(self.MW.ModelType == 'Muscle Fibers' or self.MW.ModelType == 'Motor Unit'):
            iter_num = self.mf_setValue_num
            idx = self.mf_idx 
            k = 0            
            for i in range(iter_num):
                item = QTableWidgetItem(str(self.mf_setValue[k][i-idx[k]]))
                self.mf_tableWidget[k].setItem(i-idx[k], 1, item)
                if(i == idx[k+1]-1):
                    k+=1
                    
        # tab enable
        if(self.MW.ModelType == 'Motoneuron'):
            self.uid.tabWidget.setTabEnabled(1, True)
            self.uid.tabWidget.setTabEnabled(2, False)
            self.uid.tabWidget.setCurrentIndex(0)
            
        elif(self.MW.ModelType == 'Muscle Fibers'):
            self.uid.tabWidget.setTabEnabled(1, False)
            self.uid.tabWidget.setTabEnabled(2, True)
            self.uid.tabWidget.setCurrentIndex(0) 
            
        elif(self.MW.ModelType == 'Motor Unit'):
            self.uid.tabWidget.setTabEnabled(1, True)
            self.uid.tabWidget.setTabEnabled(2, True)
            self.uid.tabWidget.setCurrentIndex(0)
            
        self.checkValue()

    def setValue(self):
        ## set the value
        # time
        self.t_stop = float(self.uid.lineEdit_2.text())
        
        # warning of changing dt 
        if(round(self.t_dt, 10) != round(float(self.uid.lineEdit_3.text()), 10)):   
            self.MW.setTextEdit("[Warning] Please regenerate input signals for new sample time.")
        self.t_dt = float(self.uid.lineEdit_3.text())
        self.t_pt = float(self.uid.lineEdit_4.text())
        self.t_setTable = self.t_start, self.t_stop, self.t_dt, self.t_pt

        if(self.MW.ModelType == 'Motoneuron'):
            self.MW.MN.setIntegrationEnv(self.t_start, self.t_stop, self.t_dt, self.t_pt)
        elif(self.MW.ModelType == 'Muscle Fibers'):
            self.MW.MF.setIntegrationEnv(self.t_start, self.t_stop, self.t_dt, self.t_pt)
        elif(self.MW.ModelType == 'Motor Unit'):
            self.MW.MU.setIntegrationEnv(self.t_start, self.t_stop, self.t_dt, self.t_pt)
            
        # tableWidget
        if(self.MW.ModelType == 'Motoneuron' or self.MW.ModelType == 'Motor Unit'):
            iter_num = self.mn_setValue_num
            idx = self.mn_idx
            k = 0            
            for i in range(iter_num):
                item = self.mn_tableWidget[k].item(i-idx[k], 1)
                self.mn_setValue[k][i-idx[k]] = float(item.text())
                if(i == idx[k+1]-1):
                    k+=1
        
        if(self.MW.ModelType == 'Muscle Fibers' or self.MW.ModelType == 'Motor Unit'):
            iter_num = self.mf_setValue_num 
            idx = self.mf_idx
            k = 0            
            for i in range(iter_num):
                item = self.mf_tableWidget[k].item(i-idx[k], 1)
                self.mf_setValue[k][i-idx[k]] = float(item.text())
                if(i == idx[k+1]-1):
                    k+=1
        
        # Text message        
        if(self.MW.ModelType == 'Motoneuron'):
            self.MW.MN.setInitialValues(self.mn_setValue)
            self.MW.setTextEdit("[Motoneuron] Simulation conditions set.")
        elif(self.MW.ModelType == 'Muscle Fibers'):
            self.MW.MF.setInitialValues(self.mf_setValue)
            self.MW.setTextEdit("[Muscle Fibers] Simulation conditions set.")
        elif(self.MW.ModelType == 'Motor Unit'):
            self.MW.MU.setInitialValues(self.mn_setValue, self.mf_setValue)
            self.MW.setTextEdit("[Motor Unit] Simulation conditions set.")
        
# Input Signals Window class
class SignalGeneratorWindow(QDialog):
    def __init__(self, MW):
        super(SignalGeneratorWindow, self).__init__()
        self.MW = MW
        #figure
        self.is_fig = None
        self.syn_fig = None
        self.sp_fig = None
        self.xm_fig = None
        # synapse channel enable
        self.se_ch=True
        self.si_ch=True
        self.de_ch=True
        self.di_ch=True
        self.syn_ch = [self.se_ch, self.si_ch, self.de_ch, self.di_ch]
        self.expIs_time = []
        self.expIs = []
        self.exp_syn_time = [0,0,0,0]
        self.exp_syn = [0,0,0,0]
        self.ExpSpike = []
        self.expXm_time = []
        self.expXm = []
        # simulation time
        self.t_start, self.t_final, self.t_dt, self.t_pt = self.MW.ISW.t_setTable
        period = self.t_final/2
        # generated signal type
        self.gen_isType = 'Step'
        self.gen_synType = ['Step','Step','Step','Step']
        self.gen_spikeType = 'User'
        self.gen_xmType = 'isometric'
        
        ## initialization
        # set value
        self.setTable = [0., 20., 5000., 100., 1500., 10., 0., -8., 0., -16., 600., 800.]
        self.isF_setValue = [6.795, 0.5, 1000., 1300., -0.3, 3000., 3300., 0.37, 5000., 9000., 0.45, 6500., 6800., -0.37, 9000., 9300., 17.]
        sesynF_setValue = [0.029, 0.02, 1000., 1300., -0.028, 3000., 3300., 0.013, 5000., 9000., 0.03, 6500., 6800., -0.1, 9000., 9300., 1., 0.5, 0.03]
        sisynF_setValue = [0.029, 0.02, 1000., 1300., -0.028, 3000., 3300., 0.013, 5000., 9000., 0.03, 6500., 6800., -0.1, 9000., 9300., 1., 2., 0.06]
        desynF_setValue = [0.029, 0.02, 1000., 1300., -0.028, 3000., 3300., 0.013, 5000., 9000., 0.03, 6500., 6800., -0.1, 9000., 9300., 1., 0.5, 0.03]
        disynF_setValue = [0.029, 0.02, 1000., 1300., -0.028, 3000., 3300., 0.013, 5000., 9000., 0.03, 6500., 6800., -0.1, 9000., 9300., 1., 2., 0.06]
        self.synF_setValue = [sesynF_setValue,sisynF_setValue,desynF_setValue,disynF_setValue]
        self.synF_idx = [0,19,38,57,76,77]
        self.synF_num = 76
        sesynT_setValue = [0., 0.1, period, 0.5, 0.03]
        sisynT_setValue = [0.1, 0,  period, 2.,  0.06]
        desynT_setValue = [0., 0.1, period, 0.5, 0.03]
        disynT_setValue = [0., 0.1,  period, 2.,  0.06]
        self.synT_setValue = [sesynT_setValue, sisynT_setValue, desynT_setValue, disynT_setValue]
        self.synT_idx = [0,5,10,15,20,21]
        self.synT_num = 20
        syn_idx = [0,19,24,43,48,67,72,91,96,97]
        syn_num = 96
        
        # default value
        self.default_num = 12
        self.defaultTable = tuple(self.setTable)
        self.isF_defValue = tuple(self.isF_setValue)
        self.synF_defValue = copy.deepcopy(self.synF_setValue)
        self.synT_defValue = copy.deepcopy(self.synT_setValue) 
        
        # generated signals
        self.gen_isF_heav = self.isF_setValue[:]
        self.gen_is_iValue = 0.
        self.gen_is_pValue = 20.
        self.gen_is_period = period
        gen_sesynF_heav = sesynF_setValue[:17]
        gen_sisynF_heav = sisynF_setValue[:17]
        gen_desynF_heav = desynF_setValue[:17]
        gen_disynF_heav = disynF_setValue[:17]
        self.gen_spike_t1 = 100.
        self.gen_spike_t2 = 1500.
        self.gen_spike_hz = 10.
        self.gen_spike_scale = 0.
        self.gen_xm_value = -8.
        self.gen_xm_value1 = 0.
        self.gen_xm_value2 = -16.
        self.gen_xm_t1 = 600.
        self.gen_xm_t2 = 800.
        self.gen_synF_heav=[gen_sesynF_heav, gen_sisynF_heav, gen_desynF_heav, gen_disynF_heav]
        self.gen_synF_tau=[0.5, 2., 0.5, 2.]
        self.gen_synF_std_max=[0.03, 0.06, 0.03, 0.06]
        self.gen_synT_iValue=[0., 0., 0., 0.]
        self.gen_synT_pValue=[0.12, 0.12, 0.12, 0.12] 
        self.gen_synT_period=[period, period, period, period]
        self.gen_synT_tau=[0.5, 2., 0.5, 2.]
        self.gen_synT_std_max =[0.03, 0.06, 0.03, 0.06]        

        # Generator class
        self.ISG = InputSignalGenerator('Step', self.t_final, self.t_dt, self.gen_isF_heav)
        self.ISG.genSignal()
        self.s_SCSG = SynConSignalGenerator('Excitatory', 'Step', self.t_final, self.t_dt, gen_sesynF_heav)
        self.s_SCSG.genSignal()
        self.s_SCSG.setValue('Inhibitory', 'Step', self.t_final, self.t_dt, gen_sisynF_heav)
        self.s_SCSG.genSignal()
        self.d_SCSG = SynConSignalGenerator('Excitatory', 'Step', self.t_final, self.t_dt, gen_desynF_heav)
        self.d_SCSG.genSignal()
        self.d_SCSG.setValue('Inhibitory', 'Step', self.t_final, self.t_dt, gen_disynF_heav)
        self.d_SCSG.genSignal()
        self.SSG = SpikeSignalGenerator('User', self.t_final, self.t_dt, self.setTable[3], self.setTable[4], self.setTable[5], self.setTable[6])
        self.SSG.genSignal()
        self.XSG = XmSignalGenerator('Isometric', self.t_final, self.t_dt, self.setTable[7])
        self.XSG.genSignal()
        
        # subscript, superscript, Greek letter with unicode and unit label
        super_2 = u'\u00B2'
        mS_cm_2 = ' (mS/cm'+super_2+')'
        nA = ' (nA)'
        ms = ' (ms)'
        sub_0 = u'\u2080'
        sub_1 = u'\u2081'
        sub_2 = u'\u2082'
        sub_3 = u'\u2083'
        sub_4 = u'\u2084'
        sub_5 = u'\u2085'
        tau = u'\u03C4'
        
        # subscript with HTML
        sub_max = '<sub>max</sub>'
        sub_p = '<sub>p</sub>'
        sub_on = '<sub>on</sub>'
        sub_off = '<sub>off</sub>'
        
        ## GUI
        self.uid = Ui_SGD()
        self.uid.setupUi(self)
        self.uid.sp_lineEdit.setReadOnly(True)
        self.uid.xm_lineEdit.setReadOnly(True)
        
        #groupBox
        self.groupBox_syn = [self.uid.groupBox_sesyn,self.uid.groupBox_sisyn,self.uid.groupBox_desyn,self.uid.groupBox_disyn]
        
        # tableWidget
        self.syn_tableWidget_F = [self.uid.tableWidget_sesynF,self.uid.tableWidget_sisynF,self.uid.tableWidget_desynF,self.uid.tableWidget_disynF]
        self.syn_tableWidget_T = [self.uid.tableWidget_sesynT,self.uid.tableWidget_sisynT,self.uid.tableWidget_desynT,self.uid.tableWidget_disynT]
        
        # lineEdit
        self.syn_lineEdit = [self.uid.sesyn_lineEdit,self.uid.sisyn_lineEdit,self.uid.desyn_lineEdit,self.uid.disyn_lineEdit]
        
        # pushButton
        self.syn_load_pushButton = [self.uid.sesyn_load_pushButton,self.uid.sisyn_load_pushButton,self.uid.desyn_load_pushButton,self.uid.disyn_load_pushButton]
        
        # Isoma radio button
        self.isFixBtn = self.uid.radioButton
        self.isTriBtn = self.uid.radioButton_2
        self.isImtBtn = self.uid.radioButton_17
        
        # Isyn noise checkbox
        self.sesyn_noise = self.uid.sesyn_noise
        self.sisyn_noise = self.uid.sisyn_noise
        self.desyn_noise = self.uid.desyn_noise
        self.disyn_noise = self.uid.disyn_noise
        self.syn_noise = [self.sesyn_noise,self.sisyn_noise,self.desyn_noise,self.disyn_noise]
        
        # Isyn radio button
        syn_SeFixBtn = self.uid.radioButton_9 
        syn_SeTriBtn = self.uid.radioButton_10 
        syn_SeImtBtn = self.uid.radioButton_18
        syn_SiFixBtn = self.uid.radioButton_11
        syn_SiTriBtn = self.uid.radioButton_12
        syn_SiImtBtn = self.uid.radioButton_19
        syn_DeFixBtn = self.uid.radioButton_13 
        syn_DeTriBtn = self.uid.radioButton_14 
        syn_DeImtBtn = self.uid.radioButton_20 
        syn_DiFixBtn = self.uid.radioButton_15 
        syn_DiTriBtn = self.uid.radioButton_16 
        syn_DiImtBtn = self.uid.radioButton_21
        self.syn_FixBtn = [syn_SeFixBtn,syn_SiFixBtn,syn_DeFixBtn,syn_DiFixBtn]
        self.syn_TriBtn = [syn_SeTriBtn,syn_SiTriBtn,syn_DeTriBtn,syn_DiTriBtn]
        self.syn_ImtBtn = [syn_SeImtBtn,syn_SiImtBtn,syn_DeImtBtn,syn_DiImtBtn]
        
        # Iaxon radio button
        self.spUserBtn = self.uid.radioButton_3
        self.spExpBtn = self.uid.radioButton_4
        
        # Xm radio button
        self.xmIsomeBtn = self.uid.radioButton_5 
        self.xmIsokiBtn = self.uid.radioButton_6 
        self.xmDyBtn = self.uid.radioButton_7 
        self.xmExpBtn = self.uid.radioButton_8 
        
        # parameter label 
        isF_parameter_1 = QLabel('I' + sub_0 + nA)
        isF_parameter_2 = QLabel('I' + sub_p + sub_1 +nA)
        isF_parameter_3 = QLabel('T' + sub_on + ' ' + sub_p + sub_1 + ms)
        isF_parameter_4 = QLabel('T' + sub_off + ' '+ sub_p + sub_1 + ms)
        isF_parameter_5 = QLabel('I' + sub_p + sub_2 + nA)
        isF_parameter_6 = QLabel('T' + sub_on + ' ' + sub_p + sub_2 + ms)
        isF_parameter_7 = QLabel('T' + sub_off + ' '+ sub_p + sub_2 + ms)
        isF_parameter_8 = QLabel('I' + sub_p + sub_3 + nA)
        isF_parameter_9 = QLabel('T' + sub_on + ' ' + sub_p + sub_3 + ms)
        isF_parameter_10 = QLabel('T' + sub_off + ' ' + sub_p + sub_3 + ms)
        isF_parameter_11 = QLabel('I' + sub_p + sub_4 + nA)
        isF_parameter_12 = QLabel('T' + sub_on + ' ' + sub_p + sub_4 + ms)
        isF_parameter_13 = QLabel('T' + sub_off + ' ' + sub_p + sub_4 + ms)
        isF_parameter_14 = QLabel('I' + sub_p + sub_5 + nA)
        isF_parameter_15 = QLabel('T' + sub_on + ' ' + sub_p + sub_5 + ms)
        isF_parameter_16 = QLabel('T' + sub_off + ' ' + sub_p + sub_5 + ms)
        isF_parameter_17 = QLabel('scale factor')
        sesynF_parameter_1 = QLabel('G' + sub_0 + mS_cm_2)
        sesynF_parameter_2 = QLabel('G' + sub_p + sub_1 + mS_cm_2)
        sesynF_parameter_3 = QLabel('T' + sub_on + ' ' + sub_p + sub_1 + ms)
        sesynF_parameter_4 = QLabel('T' + sub_off + ' '+ sub_p + sub_1 + ms)
        sesynF_parameter_5 = QLabel('G' + sub_p + sub_2 + mS_cm_2)
        sesynF_parameter_6 = QLabel('T' + sub_on + ' ' + sub_p + sub_2 + ms)
        sesynF_parameter_7 = QLabel('T' + sub_off + ' '+ sub_p + sub_2 + ms)
        sesynF_parameter_8 = QLabel('G' + sub_p + sub_3 + mS_cm_2)
        sesynF_parameter_9 = QLabel('T' + sub_on + ' ' + sub_p + sub_3 + ms)
        sesynF_parameter_10 = QLabel('T' + sub_off + ' ' + sub_p + sub_3 + ms)
        sesynF_parameter_11 = QLabel('G' + sub_p + sub_4 + mS_cm_2)
        sesynF_parameter_12 = QLabel('T' + sub_on + ' ' + sub_p + sub_4 + ms)
        sesynF_parameter_13 = QLabel('T' + sub_off + ' ' + sub_p + sub_4 + ms)
        sesynF_parameter_14 = QLabel('G' + sub_p + sub_5 + mS_cm_2)
        sesynF_parameter_15 = QLabel('T' + sub_on + ' ' + sub_p + sub_5 + ms)
        sesynF_parameter_16 = QLabel('T' + sub_off + ' ' + sub_p + sub_5 + ms)
        sesynF_parameter_17 = QLabel('scale factor')
        sesynF_parameter_18 = QLabel(tau + ms)
        sesynF_parameter_19 = QLabel('std' + sub_max)
        sesynF_parameter_18.setEnabled(False)
        sesynF_parameter_19.setEnabled(False)
        sesynT_parameter_1 = QLabel('G' + sub_0 + mS_cm_2)
        sesynT_parameter_2 = QLabel('G' + sub_p + mS_cm_2)
        sesynT_parameter_3 = QLabel('Duration' + ms)
        sesynT_parameter_4 = QLabel(tau + ms)
        sesynT_parameter_5 = QLabel('std' + sub_max)
        sesynT_parameter_4.setEnabled(False)
        sesynT_parameter_5.setEnabled(False)
        sisynF_parameter_1 = QLabel('G' + sub_0 + mS_cm_2)
        sisynF_parameter_2 = QLabel('G' + sub_p + sub_1 +mS_cm_2)
        sisynF_parameter_3 = QLabel('T' + sub_on + ' ' + sub_p + sub_1 + ms)
        sisynF_parameter_4 = QLabel('T' + sub_off + ' '+ sub_p + sub_1 + ms)
        sisynF_parameter_5 = QLabel('G' + sub_p + sub_2 + mS_cm_2)
        sisynF_parameter_6 = QLabel('T' + sub_on + ' ' + sub_p + sub_2 + ms)
        sisynF_parameter_7 = QLabel('T' + sub_off + ' '+ sub_p + sub_2 + ms)
        sisynF_parameter_8 = QLabel('G' + sub_p + sub_3 + mS_cm_2)
        sisynF_parameter_9 = QLabel('T' + sub_on + ' ' + sub_p + sub_3 + ms)
        sisynF_parameter_10 = QLabel('T' + sub_off + ' ' + sub_p + sub_3 + ms)
        sisynF_parameter_11 = QLabel('G' + sub_p + sub_4 + mS_cm_2)
        sisynF_parameter_12 = QLabel('T' + sub_on + ' ' + sub_p + sub_4 + ms)
        sisynF_parameter_13 = QLabel('T' + sub_off + ' ' + sub_p + sub_4 + ms)
        sisynF_parameter_14 = QLabel('G' + sub_p + sub_5 + mS_cm_2)
        sisynF_parameter_15 = QLabel('T' + sub_on + ' ' + sub_p + sub_5 + ms)
        sisynF_parameter_16 = QLabel('T' + sub_off + ' ' + sub_p + sub_5 + ms)
        sisynF_parameter_17 = QLabel('scale factor')
        sisynF_parameter_18 = QLabel(tau + ms)
        sisynF_parameter_19 = QLabel('std' + sub_max)
        sisynF_parameter_18.setEnabled(False)
        sisynF_parameter_19.setEnabled(False)
        sisynT_parameter_1 = QLabel('G' + sub_0 + mS_cm_2)
        sisynT_parameter_2 = QLabel('G' + sub_p + mS_cm_2)
        sisynT_parameter_3 = QLabel('Duration' + ms)
        sisynT_parameter_4 = QLabel(tau + ms)
        sisynT_parameter_5 = QLabel('std' + sub_max)
        sisynT_parameter_4.setEnabled(False)
        sisynT_parameter_5.setEnabled(False)
        desynF_parameter_1 = QLabel('G' + sub_0 + mS_cm_2)
        desynF_parameter_2 = QLabel('G' + sub_p + sub_1 +mS_cm_2)
        desynF_parameter_3 = QLabel('T' + sub_on + ' ' + sub_p + sub_1 + ms)
        desynF_parameter_4 = QLabel('T' + sub_off + ' '+ sub_p + sub_1 + ms)
        desynF_parameter_5 = QLabel('G' + sub_p + sub_2 + mS_cm_2)
        desynF_parameter_6 = QLabel('T' + sub_on + ' ' + sub_p + sub_2 + ms)
        desynF_parameter_7 = QLabel('T' + sub_off + ' '+ sub_p + sub_2 + ms)
        desynF_parameter_8 = QLabel('G' + sub_p + sub_3 + mS_cm_2)
        desynF_parameter_9 = QLabel('T' + sub_on + ' ' + sub_p + sub_3 + ms)
        desynF_parameter_10 = QLabel('T' + sub_off + ' ' + sub_p + sub_3 + ms)
        desynF_parameter_11 = QLabel('G' + sub_p + sub_4 + mS_cm_2)
        desynF_parameter_12 = QLabel('T' + sub_on + ' ' + sub_p + sub_4 + ms)
        desynF_parameter_13 = QLabel('T' + sub_off + ' ' + sub_p + sub_4 + ms)
        desynF_parameter_14 = QLabel('G' + sub_p + sub_5 + mS_cm_2)
        desynF_parameter_15 = QLabel('T' + sub_on + ' ' + sub_p + sub_5 + ms)
        desynF_parameter_16 = QLabel('T' + sub_off + ' ' + sub_p + sub_5 + ms)
        desynF_parameter_17 = QLabel('scale factor')
        desynF_parameter_18 = QLabel(tau + ms)
        desynF_parameter_19 = QLabel('std' + sub_max)
        desynF_parameter_18.setEnabled(False)
        desynF_parameter_19.setEnabled(False)
        desynT_parameter_1 = QLabel('G' + sub_0 + mS_cm_2)
        desynT_parameter_2 = QLabel('G' + sub_p + mS_cm_2)
        desynT_parameter_3 = QLabel('Duration' + ms)
        desynT_parameter_4 = QLabel(tau + ms)
        desynT_parameter_5 = QLabel('std' + sub_max)
        desynT_parameter_4 = QLabel(tau + ms)
        desynT_parameter_5 = QLabel('std' + sub_max)
        desynT_parameter_4.setEnabled(False)
        desynT_parameter_5.setEnabled(False)
        disynF_parameter_1 = QLabel('G' + sub_0 + mS_cm_2)
        disynF_parameter_2 = QLabel('G' + sub_p + sub_1 +mS_cm_2)
        disynF_parameter_3 = QLabel('T' + sub_on + ' ' + sub_p + sub_1 + ms)
        disynF_parameter_4 = QLabel('T' + sub_off + ' '+ sub_p + sub_1 + ms)
        disynF_parameter_5 = QLabel('G' + sub_p + sub_2 + mS_cm_2)
        disynF_parameter_6 = QLabel('T' + sub_on + ' ' + sub_p + sub_2 + ms)
        disynF_parameter_7 = QLabel('T' + sub_off + ' '+ sub_p + sub_2 + ms)
        disynF_parameter_8 = QLabel('G' + sub_p + sub_3 + mS_cm_2)
        disynF_parameter_9 = QLabel('T' + sub_on + ' ' + sub_p + sub_3 + ms)
        disynF_parameter_10 = QLabel('T' + sub_off + ' ' + sub_p + sub_3 + ms)
        disynF_parameter_11 = QLabel('G' + sub_p + sub_4 + mS_cm_2)
        disynF_parameter_12 = QLabel('T' + sub_on + ' ' + sub_p + sub_4 + ms)
        disynF_parameter_13 = QLabel('T' + sub_off + ' ' + sub_p + sub_4 + ms)
        disynF_parameter_14 = QLabel('G' + sub_p + sub_5 + mS_cm_2)
        disynF_parameter_15 = QLabel('T' + sub_on + ' ' + sub_p + sub_5 + ms)
        disynF_parameter_16 = QLabel('T' + sub_off + ' ' + sub_p + sub_5 + ms)
        disynF_parameter_17 = QLabel('scale factor')
        disynF_parameter_18 = QLabel(tau + ms)
        disynF_parameter_19 = QLabel('std' + sub_max)
        disynF_parameter_18.setEnabled(False)
        disynF_parameter_19.setEnabled(False)
        disynT_parameter_1 = QLabel('G' + sub_0 + mS_cm_2)
        disynT_parameter_2 = QLabel('G' + sub_p + mS_cm_2)
        disynT_parameter_3 = QLabel('Duration' + ms)
        disynT_parameter_4 = QLabel(tau + ms)
        disynT_parameter_5 = QLabel('std' + sub_max)
        disynT_parameter_4.setEnabled(False)
        disynT_parameter_5.setEnabled(False)
        
        isF_parameter = []
        isF_parameter.append(isF_parameter_1)
        isF_parameter.append(isF_parameter_2)
        isF_parameter.append(isF_parameter_3)
        isF_parameter.append(isF_parameter_4)
        isF_parameter.append(isF_parameter_5)
        isF_parameter.append(isF_parameter_6)
        isF_parameter.append(isF_parameter_7)
        isF_parameter.append(isF_parameter_8)
        isF_parameter.append(isF_parameter_9)
        isF_parameter.append(isF_parameter_10)
        isF_parameter.append(isF_parameter_11)
        isF_parameter.append(isF_parameter_12)
        isF_parameter.append(isF_parameter_13)
        isF_parameter.append(isF_parameter_14)
        isF_parameter.append(isF_parameter_15)
        isF_parameter.append(isF_parameter_16)
        isF_parameter.append(isF_parameter_17)
        sesynF_parameter = []
        sesynF_parameter.append(sesynF_parameter_1)
        sesynF_parameter.append(sesynF_parameter_2)
        sesynF_parameter.append(sesynF_parameter_3)
        sesynF_parameter.append(sesynF_parameter_4)
        sesynF_parameter.append(sesynF_parameter_5)
        sesynF_parameter.append(sesynF_parameter_6)
        sesynF_parameter.append(sesynF_parameter_7)
        sesynF_parameter.append(sesynF_parameter_8)
        sesynF_parameter.append(sesynF_parameter_9)
        sesynF_parameter.append(sesynF_parameter_10)
        sesynF_parameter.append(sesynF_parameter_11)
        sesynF_parameter.append(sesynF_parameter_12)
        sesynF_parameter.append(sesynF_parameter_13)
        sesynF_parameter.append(sesynF_parameter_14)
        sesynF_parameter.append(sesynF_parameter_15)
        sesynF_parameter.append(sesynF_parameter_16)
        sesynF_parameter.append(sesynF_parameter_17)
        sesynF_parameter.append(sesynF_parameter_18)
        sesynF_parameter.append(sesynF_parameter_19)
        sesynT_parameter = []
        sesynT_parameter.append(sesynT_parameter_1)
        sesynT_parameter.append(sesynT_parameter_2)
        sesynT_parameter.append(sesynT_parameter_3)
        sesynT_parameter.append(sesynT_parameter_4)
        sesynT_parameter.append(sesynT_parameter_5)
        sisynF_parameter = []
        sisynF_parameter.append(sisynF_parameter_1)
        sisynF_parameter.append(sisynF_parameter_2)
        sisynF_parameter.append(sisynF_parameter_3)
        sisynF_parameter.append(sisynF_parameter_4)
        sisynF_parameter.append(sisynF_parameter_5)
        sisynF_parameter.append(sisynF_parameter_6)
        sisynF_parameter.append(sisynF_parameter_7)
        sisynF_parameter.append(sisynF_parameter_8)
        sisynF_parameter.append(sisynF_parameter_9)
        sisynF_parameter.append(sisynF_parameter_10)
        sisynF_parameter.append(sisynF_parameter_11)
        sisynF_parameter.append(sisynF_parameter_12)
        sisynF_parameter.append(sisynF_parameter_13)
        sisynF_parameter.append(sisynF_parameter_14)
        sisynF_parameter.append(sisynF_parameter_15)
        sisynF_parameter.append(sisynF_parameter_16)
        sisynF_parameter.append(sisynF_parameter_17)
        sisynF_parameter.append(sisynF_parameter_18)
        sisynF_parameter.append(sisynF_parameter_19)
        sisynT_parameter = []
        sisynT_parameter.append(sisynT_parameter_1)
        sisynT_parameter.append(sisynT_parameter_2)
        sisynT_parameter.append(sisynT_parameter_3)
        sisynT_parameter.append(sisynT_parameter_4)
        sisynT_parameter.append(sisynT_parameter_5)
        desynF_parameter = []
        desynF_parameter.append(desynF_parameter_1)
        desynF_parameter.append(desynF_parameter_2)
        desynF_parameter.append(desynF_parameter_3)
        desynF_parameter.append(desynF_parameter_4)
        desynF_parameter.append(desynF_parameter_5)
        desynF_parameter.append(desynF_parameter_6)
        desynF_parameter.append(desynF_parameter_7)
        desynF_parameter.append(desynF_parameter_8)
        desynF_parameter.append(desynF_parameter_9)
        desynF_parameter.append(desynF_parameter_10)
        desynF_parameter.append(desynF_parameter_11)
        desynF_parameter.append(desynF_parameter_12)
        desynF_parameter.append(desynF_parameter_13)
        desynF_parameter.append(desynF_parameter_14)
        desynF_parameter.append(desynF_parameter_15)
        desynF_parameter.append(desynF_parameter_16)
        desynF_parameter.append(desynF_parameter_17)
        desynF_parameter.append(desynF_parameter_18)
        desynF_parameter.append(desynF_parameter_19)
        desynT_parameter = []
        desynT_parameter.append(desynT_parameter_1)
        desynT_parameter.append(desynT_parameter_2)
        desynT_parameter.append(desynT_parameter_3)
        desynT_parameter.append(desynT_parameter_4)
        desynT_parameter.append(desynT_parameter_5)
        disynF_parameter = []
        disynF_parameter.append(disynF_parameter_1)
        disynF_parameter.append(disynF_parameter_2)
        disynF_parameter.append(disynF_parameter_3)
        disynF_parameter.append(disynF_parameter_4)
        disynF_parameter.append(disynF_parameter_5)
        disynF_parameter.append(disynF_parameter_6)
        disynF_parameter.append(disynF_parameter_7)
        disynF_parameter.append(disynF_parameter_8)
        disynF_parameter.append(disynF_parameter_9)
        disynF_parameter.append(disynF_parameter_10)
        disynF_parameter.append(disynF_parameter_11)
        disynF_parameter.append(disynF_parameter_12)
        disynF_parameter.append(disynF_parameter_13)
        disynF_parameter.append(disynF_parameter_14)
        disynF_parameter.append(disynF_parameter_15)
        disynF_parameter.append(disynF_parameter_16)
        disynF_parameter.append(disynF_parameter_17)
        disynF_parameter.append(disynF_parameter_18)
        disynF_parameter.append(disynF_parameter_19)
        disynT_parameter = []
        disynT_parameter.append(disynT_parameter_1)
        disynT_parameter.append(disynT_parameter_2)
        disynT_parameter.append(disynT_parameter_3)
        disynT_parameter.append(disynT_parameter_4)
        disynT_parameter.append(disynT_parameter_5)

        ## tableWidget
        # row, column
        self.isF_idx = len(isF_parameter)
        self.uid.tableWidget_isF.setColumnCount(3)
        self.uid.tableWidget_isF.setRowCount(self.isF_idx)
        
        col_num = 3
        iter_num = len(self.syn_tableWidget_F)
        for i in range(iter_num):
            self.syn_tableWidget_F[i].setColumnCount(col_num)
            self.syn_tableWidget_F[i].setRowCount(len(self.synF_setValue[i]))
        
        iter_num = len(self.syn_tableWidget_T)
        for i in range(iter_num):
            self.syn_tableWidget_T[i].setColumnCount(col_num)
            self.syn_tableWidget_T[i].setRowCount(len(self.synT_setValue[i]))
        
        # column width
        iter_num = len(self.syn_tableWidget_F)
        for i in range(iter_num):
            tableWidget_header = self.syn_tableWidget_F[i].horizontalHeader()
            tableWidget_header.setDefaultSectionSize(85) 
        
        iter_num = len(self.syn_tableWidget_T)
        for i in range(iter_num):
            tableWidget_header = self.syn_tableWidget_T[i].horizontalHeader()
            tableWidget_header.setDefaultSectionSize(90)    
        
        # parameter lebel
        for i in range(self.isF_idx):
            self.uid.tableWidget_isF.setCellWidget(i, 0, isF_parameter[i])
        synF_idx = len(sesynF_parameter)
        synT_idx = len(sesynT_parameter)
        
        for i in range(synF_idx):
            self.uid.tableWidget_sesynF.setCellWidget(i, 0, sesynF_parameter[i])
            self.uid.tableWidget_sisynF.setCellWidget(i, 0, sisynF_parameter[i])
            self.uid.tableWidget_desynF.setCellWidget(i, 0, desynF_parameter[i])
            self.uid.tableWidget_disynF.setCellWidget(i, 0, disynF_parameter[i])
        
        for i in range(synT_idx):
            self.uid.tableWidget_sesynT.setCellWidget(i, 0, sesynT_parameter[i])
            self.uid.tableWidget_sisynT.setCellWidget(i, 0, sisynT_parameter[i])
            self.uid.tableWidget_desynT.setCellWidget(i, 0, desynT_parameter[i])
            self.uid.tableWidget_disynT.setCellWidget(i, 0, disynT_parameter[i])
        
        # deactivation cell
        for i in range(self.isF_idx):
            item = QTableWidgetItem()
            item.setFlags(Qt.ItemIsEnabled)
            self.uid.tableWidget_isF.setItem(i, 2, item)
        
        iter_num = syn_num 
        idx = syn_idx
        syn_tableWidget = np.asarray(zip(self.syn_tableWidget_F + self.syn_tableWidget_T)).flatten()
        k = 0
        for i in range(iter_num):
            item = QTableWidgetItem()
            item.setFlags(Qt.ItemIsEnabled)
            syn_tableWidget[k].setItem(i-idx[k], 2, item)
            if(i == idx[k+1]-1):
                k+=1
        
        # lineEdit
        self.lineEdit = []
        self.lineEdit.append(self.uid.lineEdit_1)        
        self.lineEdit.append(self.uid.lineEdit_2)
        self.lineEdit.append(self.uid.lineEdit_3)
        self.lineEdit.append(self.uid.lineEdit_4)
        self.lineEdit.append(self.uid.lineEdit_5)
        self.lineEdit.append(self.uid.lineEdit_6)
        self.lineEdit.append(self.uid.lineEdit_7)
        self.lineEdit.append(self.uid.lineEdit_8)
        self.lineEdit.append(self.uid.lineEdit_9)
        self.lineEdit.append(self.uid.lineEdit_10)
        self.lineEdit.append(self.uid.lineEdit_11)
        self.lineEdit.append(self.uid.lineEdit_12)
        
        # checkbox
        self.checkBox = []
        self.checkBox.append(self.uid.checkBox_1)
        self.checkBox.append(self.uid.checkBox_2)
        self.checkBox.append(self.uid.checkBox_3)
        self.checkBox.append(self.uid.checkBox_4)
        self.checkBox.append(self.uid.checkBox_5)
        self.checkBox.append(self.uid.checkBox_6)
        self.checkBox.append(self.uid.checkBox_7)
        self.checkBox.append(self.uid.checkBox_8)
        self.checkBox.append(self.uid.checkBox_9)
        self.checkBox.append(self.uid.checkBox_10)
        self.checkBox.append(self.uid.checkBox_11)
        self.checkBox.append(self.uid.checkBox_12)
        
        # initial value for velocity of Xm
        self.uid.lineEdit_8.setText(str(self.setTable[7]))
        self.uid.lineEdit_9.setText(str(self.setTable[8]))
        self.uid.lineEdit_10.setText(str(self.setTable[9]))
        self.uid.lineEdit_11.setText(str(self.setTable[10]))
        self.uid.lineEdit_12.setText(str(self.setTable[11]))
        
        # set button
        self.isSetButton = self.isFixBtn
        syn_SeSetButton = syn_SeFixBtn
        syn_SiSetButton = syn_SiFixBtn
        syn_DeSetButton = syn_DeFixBtn
        syn_DiSetButton = syn_DiFixBtn     
        self.syn_setButton = [syn_SeSetButton, syn_SiSetButton, syn_DeSetButton, syn_DiSetButton]
        self.spSetButton = self.spUserBtn
        self.xmSetButton = self.xmIsomeBtn
        self.sesyn_setNoise = False
        self.sisyn_setNoise = False
        self.desyn_setNoise = False
        self.disyn_setNoise = False
        
        # set button activation
        self.isSetButton.setChecked(True)
        syn_SeSetButton.setChecked(True)
        syn_SiSetButton.setChecked(True)
        syn_DeSetButton.setChecked(True)
        syn_DiSetButton.setChecked(True)
        self.spSetButton.setChecked(True)
        self.xmSetButton.setChecked(True)
        self.sesyn_noise.setChecked(self.sesyn_setNoise)
        self.sisyn_noise.setChecked(self.sisyn_setNoise)
        self.desyn_noise.setChecked(self.desyn_setNoise)
        self.disyn_noise.setChecked(self.disyn_setNoise)
        
        # initial deactivation
        self.uid.applyButton.setEnabled(False)
        self.uid.okButton.setEnabled(False)
        self.uid.label_1.setEnabled(False)
        self.uid.label_2.setEnabled(False)
        self.uid.label_3.setEnabled(False)
        self.lineEdit[0].setEnabled(False)
        self.lineEdit[1].setEnabled(False)
        self.lineEdit[2].setEnabled(False)
        self.checkBox[0].setEnabled(False)
        self.checkBox[1].setEnabled(False)
        self.checkBox[2].setEnabled(False)
        self.uid.is_lineEdit.setEnabled(False)
        self.uid.is_load_pushButton.setEnabled(False)
        self.uid.tableWidget_sesynT.setEnabled(False)
        self.uid.sesyn_lineEdit.setEnabled(False)
        self.uid.sesyn_load_pushButton.setEnabled(False)
        self.uid.tableWidget_sisynT.setEnabled(False)
        self.uid.sisyn_lineEdit.setEnabled(False)
        self.uid.sisyn_load_pushButton.setEnabled(False)
        self.uid.tableWidget_desynT.setEnabled(False)
        self.uid.desyn_lineEdit.setEnabled(False)
        self.uid.desyn_load_pushButton.setEnabled(False)
        self.uid.tableWidget_disynT.setEnabled(False)
        self.uid.disyn_lineEdit.setEnabled(False)
        self.uid.disyn_load_pushButton.setEnabled(False)
        self.uid.sp_lineEdit.setEnabled(False)
        self.uid.sp_load_pushButton.setEnabled(False)
        self.uid.label_10.setEnabled(False)
        self.uid.label_11.setEnabled(False)
        self.uid.label_12.setEnabled(False)
        self.uid.label_13.setEnabled(False)
        self.uid.label_14.setEnabled(False)
        self.uid.velocity_label.setEnabled(False)
        self.checkBox[8].setEnabled(False)
        self.checkBox[9].setEnabled(False)
        self.checkBox[10].setEnabled(False)
        self.checkBox[11].setEnabled(False)
        self.lineEdit[8].setEnabled(False)
        self.lineEdit[9].setEnabled(False)
        self.lineEdit[10].setEnabled(False)
        self.lineEdit[11].setEnabled(False)
        self.uid.xm_lineEdit.setEnabled(False)
        self.uid.xm_load_pushButton.setEnabled(False)
        
        # GUI operation-function 
        self.uid.buttonBox.connect(self.uid.buttonBox, SIGNAL("accepted()"), self.setValue)
        self.uid.buttonBox.connect(self.uid.applyButton, SIGNAL("clicked()"), self.setValue)
        self.uid.is_load_pushButton.clicked.connect(self.openFile)
        self.uid.sesyn_load_pushButton.clicked.connect(self.openFile)
        self.uid.sisyn_load_pushButton.clicked.connect(self.openFile)
        self.uid.desyn_load_pushButton.clicked.connect(self.openFile)
        self.uid.disyn_load_pushButton.clicked.connect(self.openFile)
        self.uid.sp_load_pushButton.clicked.connect(self.openFile)
        self.uid.xm_load_pushButton.clicked.connect(self.openFile)
        self.uid.generate_pushButton.clicked.connect(self.genSignal)
        self.isFixBtn.clicked.connect(self.setObjectEnabled)
        self.isTriBtn.clicked.connect(self.setObjectEnabled)
        self.isImtBtn.clicked.connect(self.setObjectEnabled)
        syn_SeFixBtn.clicked.connect(self.setObjectEnabled)
        syn_SeTriBtn.clicked.connect(self.setObjectEnabled)
        syn_SeImtBtn.clicked.connect(self.setObjectEnabled)
        syn_SiFixBtn.clicked.connect(self.setObjectEnabled)
        syn_SiTriBtn.clicked.connect(self.setObjectEnabled)
        syn_SiImtBtn.clicked.connect(self.setObjectEnabled)
        syn_DeFixBtn.clicked.connect(self.setObjectEnabled)
        syn_DeTriBtn.clicked.connect(self.setObjectEnabled)
        syn_DeImtBtn.clicked.connect(self.setObjectEnabled)
        syn_DiFixBtn.clicked.connect(self.setObjectEnabled)
        syn_DiTriBtn.clicked.connect(self.setObjectEnabled)
        syn_DiImtBtn.clicked.connect(self.setObjectEnabled)
        self.spUserBtn.clicked.connect(self.setObjectEnabled)
        self.spExpBtn.clicked.connect(self.setObjectEnabled)
        self.xmIsomeBtn.clicked.connect(self.setObjectEnabled)
        self.xmIsokiBtn.clicked.connect(self.setObjectEnabled)
        self.xmDyBtn.clicked.connect(self.setObjectEnabled)
        self.xmExpBtn.clicked.connect(self.setObjectEnabled)
        self.uid.lineEdit_1.editingFinished.connect(self.checkValue)
        self.uid.lineEdit_2.editingFinished.connect(self.checkValue)
        self.uid.lineEdit_3.editingFinished.connect(self.checkValue)
        self.uid.lineEdit_4.editingFinished.connect(self.checkValue)
        self.uid.lineEdit_5.editingFinished.connect(self.checkValue)
        self.uid.lineEdit_6.editingFinished.connect(self.checkValue)
        self.uid.lineEdit_7.editingFinished.connect(self.checkValue)
        self.uid.lineEdit_8.editingFinished.connect(self.checkValue)
        self.uid.lineEdit_9.editingFinished.connect(self.checkValue)
        self.uid.lineEdit_10.editingFinished.connect(self.checkValue)
        self.uid.lineEdit_11.editingFinished.connect(self.checkValue)
        self.uid.lineEdit_12.editingFinished.connect(self.checkValue)
        self.uid.tableWidget_isF.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_sesynF.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_sisynF.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_desynF.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_disynF.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_sesynT.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_sisynT.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_desynT.cellChanged.connect(self.checkValue)
        self.uid.tableWidget_disynT.cellChanged.connect(self.checkValue)
        self.uid.lineEdit_9.editingFinished.connect(self.calXmVelocity)
        self.uid.lineEdit_10.editingFinished.connect(self.calXmVelocity)
        self.uid.lineEdit_11.editingFinished.connect(self.calXmVelocity)
        self.uid.lineEdit_12.editingFinished.connect(self.calXmVelocity)
        self.uid.checkBox_1.stateChanged.connect(self.returnDefaultValue)
        self.uid.checkBox_2.stateChanged.connect(self.returnDefaultValue)
        self.uid.checkBox_3.stateChanged.connect(self.returnDefaultValue)
        self.uid.checkBox_4.stateChanged.connect(self.returnDefaultValue)
        self.uid.checkBox_5.stateChanged.connect(self.returnDefaultValue)
        self.uid.checkBox_6.stateChanged.connect(self.returnDefaultValue)
        self.uid.checkBox_7.stateChanged.connect(self.returnDefaultValue)
        self.uid.checkBox_8.stateChanged.connect(self.returnDefaultValue)
        self.uid.checkBox_9.stateChanged.connect(self.returnDefaultValue)
        self.uid.checkBox_10.stateChanged.connect(self.returnDefaultValue)
        self.uid.checkBox_11.stateChanged.connect(self.returnDefaultValue)
        self.uid.checkBox_12.stateChanged.connect(self.returnDefaultValue)
        self.uid.tableWidget_isF.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_sesynF.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_sisynF.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_desynF.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_disynF.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_sesynT.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_sisynT.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_desynT.cellClicked.connect(self.returnDefaultValue)
        self.uid.tableWidget_disynT.cellClicked.connect(self.returnDefaultValue)
        self.sesyn_noise.stateChanged.connect(self.setObjectEnabled)
        self.sisyn_noise.stateChanged.connect(self.setObjectEnabled)
        self.desyn_noise.stateChanged.connect(self.setObjectEnabled)
        self.disyn_noise.stateChanged.connect(self.setObjectEnabled)

    def keyPressEvent(self, event):
        # operation of Key_Return or Key_Enter
        if(event.key() == 0x01000005 or event.key() == 0x01000004):
            for i in range(len(self.lineEdit)):            
                self.lineEdit[i].clearFocus()
        else:
            # QDialog event
            super(SignalGeneratorWindow, self).keyPressEvent(event)

    def checkValue(self, row=0, col=3):
        # tableWidget
        if(col == 1):
            tableWidget = self.sender()
                        
            if(tableWidget == self.uid.tableWidget_isF):
                item = self.uid.tableWidget_isF.item(row, col)
                defValue = self.isF_defValue[row]
                setValue = self.isF_setValue[row]
            
            for i in range(len(self.syn_tableWidget_F)):
                if(tableWidget == self.syn_tableWidget_F[i]):
                    item = tableWidget.item(row, col)
                    defValue = self.synF_defValue[i][row]
                    setValue = self.synF_setValue[i][row]
            
            for i in range(len(self.syn_tableWidget_T)):
                if(tableWidget == self.syn_tableWidget_T[i]):
                    item = tableWidget.item(row, col)
                    defValue = self.synT_defValue[i][row]
                    setValue = self.synT_setValue[i][row]
            
            try:                
                # compare with default value
                if(round(float(item.text()), 10) != round(defValue, 10)):
                    item = QTableWidgetItem()
                    icon = QIcon()
                    icon.addPixmap(QPixmap("./resources/default.png"))
                    item.setIcon(icon)
                    item.setFlags(Qt.ItemIsEnabled)
                    tableWidget.setItem(row, 2, item)
                else:
                    item = QTableWidgetItem()
                    item.setFlags(Qt.ItemIsEnabled)  
                    tableWidget.setItem(row, 2, item)
            
            except:
                item = QTableWidgetItem(str(setValue))
                tableWidget.setItem(row, 1, item)
        
        # lineEdit
        elif(col == 3):
            # compare with default value
            for i in range(self.default_num):
                try:
                    if(round(float(self.lineEdit[i].text()), 10) != round(self.defaultTable[i], 10)):
                        self.checkBox[i].setCheckState(Qt.Checked)
                        self.checkBox[i].setEnabled(True)
                        self.setObjectEnabled() 
                    else:
                        self.checkBox[i].setCheckState(Qt.Unchecked)
                        self.checkBox[i].setEnabled(False)
                
                except:
                    self.lineEdit[i].setText(str(self.setTable[i]))
    
    def returnDefaultValue(self, row=0, col=3):
        # tableWidget
        if(col == 2):
            tableWidget = self.sender()
            
            if(tableWidget == self.uid.tableWidget_isF):
                defValue = self.isF_defValue
                item = QTableWidgetItem(str(defValue[row]))
                tableWidget.setItem(row, 1, item)
            
            for i in range(len(self.syn_tableWidget_F)):
                if(tableWidget == self.syn_tableWidget_F[i]):
                    # prevent return value of SynCon noise parameters that is deactivated.
                    if(row == 17 or row == 18):
                        if(not self.syn_noise[i].isChecked()):
                            continue
                    # return the default value        
                    item = QTableWidgetItem(str(self.synF_defValue[i][row]))
                    self.syn_tableWidget_F[i].setItem(row, 1, item)
            
            for i in range(len(self.syn_tableWidget_T)):
                if(tableWidget == self.syn_tableWidget_T[i]):
                    # prevent return value of SynCon noise parameters that is deactivated.
                    if(row == 3 or row == 4):
                        if(not self.syn_noise[i].isChecked()):
                            continue                    
                    # return the default value   
                    item = QTableWidgetItem(str(self.synT_defValue[i][row]))
                    self.syn_tableWidget_T[i].setItem(row, 1, item)
        
        # lineEdit
        elif(col == 3):
            checkBox = self.sender()
            
            if(not checkBox.isChecked()):
                # find checkBox index
                for i in range(self.default_num):
                    if(self.checkBox[i] == checkBox):
                        break
                
                # return the default value    
                self.lineEdit[i].setText(str(self.defaultTable[i]))
                self.checkBox[i].setEnabled(False)
                # prevent widget focus
                if((i+1) < self.default_num):
                    self.lineEdit[i+1].clearFocus()            
            
    def setItemColorEnable(self, tableWidget, row, col, color, enable):
        # change tau, std_max label color according to noise checkbox
        # to inform activation or deactivation of value.
        item_v = tableWidget.item(row, col)
        item = QTableWidgetItem(item_v.text())
        item.setTextColor(color)
        
        if(enable == True):
            item.setFlags(item.flags() | Qt.ItemIsEnabled)
        else:
            item.setFlags(item.flags() & Qt.ItemIsEnabled)
        
        tableWidget.setItem(row, col, item)
    
    def setObjectEnabled(self):
        ## object enable setting
        # Isoma
        if(self.MW.ModelType == 'Motoneuron' or self.MW.ModelType == 'Motor Unit'):
            if(self.isFixBtn.isChecked()):
                self.uid.tableWidget_isF.setEnabled(True)
                self.uid.label_1.setEnabled(False)
                self.uid.label_2.setEnabled(False)
                self.uid.label_3.setEnabled(False)
                self.lineEdit[0].setEnabled(False)
                self.lineEdit[1].setEnabled(False)
                self.lineEdit[2].setEnabled(False)
                self.checkBox[0].setEnabled(False)
                self.checkBox[1].setEnabled(False)
                self.checkBox[2].setEnabled(False)
                self.uid.is_lineEdit.setEnabled(False)
                self.uid.is_load_pushButton.setEnabled(False)
            
            elif(self.isTriBtn.isChecked()):
                self.uid.label_1.setEnabled(True)
                self.uid.label_2.setEnabled(True)
                self.uid.label_3.setEnabled(True)
                if(self.checkBox[0].isChecked()):
                    self.checkBox[0].setEnabled(True)
                if(self.checkBox[1].isChecked()):
                    self.checkBox[1].setEnabled(True)
                if(self.checkBox[2].isChecked()):
                    self.checkBox[2].setEnabled(True)
                self.lineEdit[0].setEnabled(True)
                self.lineEdit[1].setEnabled(True)
                self.lineEdit[2].setEnabled(True)
                self.uid.tableWidget_isF.setEnabled(False)
                self.uid.is_lineEdit.setEnabled(False)
                self.uid.is_load_pushButton.setEnabled(False)
            
            elif(self.isImtBtn.isChecked()):
                self.uid.is_lineEdit.setEnabled(True)
                self.uid.is_load_pushButton.setEnabled(True)
                self.uid.label_1.setEnabled(False)
                self.uid.label_2.setEnabled(False)
                self.uid.label_3.setEnabled(False)
                self.lineEdit[0].setEnabled(False)
                self.lineEdit[1].setEnabled(False)
                self.lineEdit[2].setEnabled(False)
                self.checkBox[0].setEnabled(False)
                self.checkBox[1].setEnabled(False)
                self.checkBox[2].setEnabled(False)
                self.uid.tableWidget_isF.setEnabled(False)
            
            # Isyn
            tau = u'\u03C4'
            ms = ' (ms)'
            sub_max = '<sub>max</sub>'
            rowF = 17
            rowT = 3
            v_col = 1
            l_col = 0
            channel_enable = range(4)
            if(self.MW.PSW.channel_enable['Sesyn'] == True):
                channel_enable[0] = True
            else:
                channel_enable[0] = False
            if(self.MW.PSW.channel_enable['Sisyn'] == True):
                channel_enable[1] = True
            else:
                channel_enable[1] = False
            if(self.MW.PSW.channel_enable['Desyn'] == True):
                channel_enable[2] = True
            else:
                channel_enable[2] = False
            if(self.MW.PSW.channel_enable['Disyn'] == True):
                channel_enable[3] = True
            else:
                channel_enable[3] = False
            
            for i in range(len(channel_enable)):
                if(channel_enable[i] == True):
                    self.syn_ch[i] = True
                    self.groupBox_syn[i].setEnabled(True)
                    label_1 = QLabel(tau + ms)
                    label_2 = QLabel('std' + sub_max)
                    label_3 = QLabel(tau + ms)
                    label_4 = QLabel('std' + sub_max)
                    
                    if(self.syn_FixBtn[i].isChecked()):
                        self.syn_noise[i].setEnabled(True)
                        self.syn_tableWidget_F[i].setEnabled(True)
                        self.syn_tableWidget_T[i].setEnabled(False)
                        self.syn_lineEdit[i].setEnabled(False)
                        self.syn_load_pushButton[i].setEnabled(False)
                        act_colorF = Qt.black
                        act_colorT = Qt.darkGray
                        de_colorF = Qt.darkGray
                        de_colorT = Qt.darkGray
                    
                    elif(self.syn_TriBtn[i].isChecked()):
                        self.syn_noise[i].setEnabled(True)
                        self.syn_tableWidget_F[i].setEnabled(False)
                        self.syn_tableWidget_T[i].setEnabled(True)
                        self.syn_lineEdit[i].setEnabled(False)
                        self.syn_load_pushButton[i].setEnabled(False)
                        act_colorF = Qt.darkGray
                        act_colorT = Qt.black
                        de_colorF = Qt.darkGray
                        de_colorT = Qt.darkGray
                    
                    elif(self.syn_ImtBtn[i].isChecked()):
                        self.syn_noise[i].setEnabled(False)
                        self.syn_tableWidget_F[i].setEnabled(False)
                        self.syn_tableWidget_T[i].setEnabled(False)
                        self.syn_lineEdit[i].setEnabled(True)
                        self.syn_load_pushButton[i].setEnabled(True)
                        act_colorF = Qt.darkGray
                        act_colorT = Qt.darkGray
                        de_colorF = Qt.darkGray
                        de_colorT = Qt.darkGray
                    
                    if(self.syn_noise[i].isChecked()):
                        label_1.setEnabled(True)
                        label_2.setEnabled(True)
                        label_3.setEnabled(True)
                        label_4.setEnabled(True)
                        self.setItemColorEnable(self.syn_tableWidget_F[i], rowF,   v_col, act_colorF, True)
                        self.setItemColorEnable(self.syn_tableWidget_F[i], rowF+1, v_col, act_colorF, True)
                        self.setItemColorEnable(self.syn_tableWidget_T[i], rowT,   v_col, act_colorT, True)
                        self.setItemColorEnable(self.syn_tableWidget_T[i], rowT+1, v_col, act_colorT, True)
                    
                    else:
                        label_1.setEnabled(False)
                        label_2.setEnabled(False)
                        label_3.setEnabled(False)
                        label_4.setEnabled(False)
                        self.setItemColorEnable(self.syn_tableWidget_F[i], rowF,   v_col, de_colorF, False)
                        self.setItemColorEnable(self.syn_tableWidget_F[i], rowF+1, v_col, de_colorF, False)
                        self.setItemColorEnable(self.syn_tableWidget_T[i], rowT,   v_col, de_colorT, False)
                        self.setItemColorEnable(self.syn_tableWidget_T[i], rowT+1, v_col, de_colorT, False)
                    self.syn_tableWidget_F[i].setCellWidget(rowF,   l_col, label_1)
                    self.syn_tableWidget_F[i].setCellWidget(rowF+1, l_col, label_2)
                    self.syn_tableWidget_T[i].setCellWidget(rowT,   l_col, label_3)
                    self.syn_tableWidget_T[i].setCellWidget(rowT+1, l_col, label_4)
                
                else:
                    self.syn_ch[i] = False
                    self.groupBox_syn[i].setEnabled(False)
                    self.syn_noise[i].setEnabled(False)
                    self.setItemColorEnable(self.syn_tableWidget_F[i], rowF,   v_col, Qt.darkGray, False)
                    self.setItemColorEnable(self.syn_tableWidget_F[i], rowF+1, v_col, Qt.darkGray, False)
                    self.setItemColorEnable(self.syn_tableWidget_T[i], rowT,   v_col, Qt.darkGray, False)
                    self.setItemColorEnable(self.syn_tableWidget_T[i], rowT+1, v_col, Qt.darkGray, False)

        if(self.MW.ModelType == 'Muscle Fibers'):
            if(self.spUserBtn.isChecked()):
                self.uid.label_4.setEnabled(True)
                self.uid.label_5.setEnabled(True)
                self.uid.label_6.setEnabled(True)
                self.uid.label_7.setEnabled(True)
                if(self.checkBox[3].isChecked()):
                    self.checkBox[3].setEnabled(True)
                if(self.checkBox[4].isChecked()):  
                    self.checkBox[4].setEnabled(True)
                if(self.checkBox[5].isChecked()):
                    self.checkBox[5].setEnabled(True)
                if(self.checkBox[6].isChecked()):
                    self.checkBox[6].setEnabled(True)
                self.lineEdit[3].setEnabled(True)
                self.lineEdit[4].setEnabled(True)
                self.lineEdit[5].setEnabled(True)
                self.lineEdit[6].setEnabled(True)
                self.uid.sp_lineEdit.setEnabled(False)
                self.uid.sp_load_pushButton.setEnabled(False)
            
            elif(self.spExpBtn.isChecked()):
                self.uid.sp_lineEdit.setEnabled(True)
                self.uid.sp_load_pushButton.setEnabled(True)
                self.uid.label_4.setEnabled(False)
                self.uid.label_5.setEnabled(False)
                self.uid.label_6.setEnabled(False)
                self.uid.label_7.setEnabled(False)
                self.checkBox[3].setEnabled(False)
                self.checkBox[4].setEnabled(False)
                self.checkBox[5].setEnabled(False)
                self.checkBox[6].setEnabled(False)
                self.lineEdit[3].setEnabled(False)
                self.lineEdit[4].setEnabled(False)
                self.lineEdit[5].setEnabled(False)
                self.lineEdit[6].setEnabled(False)
                
        if(self.MW.ModelType == 'Muscle Fibers' or self.MW.ModelType == 'Motor Unit'):
            if(self.xmIsomeBtn.isChecked()):
                self.uid.label_9.setEnabled(True)
                self.lineEdit[7].setEnabled(True)
                if(self.checkBox[7].isChecked()):
                    self.checkBox[7].setEnabled(True)
                self.uid.label_10.setEnabled(False)
                self.uid.label_11.setEnabled(False)
                self.uid.label_12.setEnabled(False)
                self.uid.label_13.setEnabled(False)
                self.uid.label_14.setEnabled(False)
                self.uid.velocity_label.setEnabled(False)
                self.checkBox[8].setEnabled(False)
                self.checkBox[9].setEnabled(False)
                self.checkBox[10].setEnabled(False)
                self.checkBox[11].setEnabled(False)
                self.lineEdit[8].setEnabled(False)
                self.lineEdit[9].setEnabled(False)
                self.lineEdit[10].setEnabled(False)
                self.lineEdit[11].setEnabled(False)
                self.uid.xm_lineEdit.setEnabled(False)
                self.uid.xm_load_pushButton.setEnabled(False)  
            
            elif(self.xmIsokiBtn.isChecked()):
                self.uid.label_9.setEnabled(False)
                self.lineEdit[7].setEnabled(False)
                self.checkBox[7].setEnabled(False)
                self.uid.label_10.setEnabled(True)
                self.uid.label_11.setEnabled(True)
                self.uid.label_12.setEnabled(True)
                self.uid.label_13.setEnabled(True)
                self.uid.label_14.setEnabled(True)
                self.uid.velocity_label.setEnabled(True)
                if(self.checkBox[8].isChecked()):
                    self.checkBox[8].setEnabled(True)
                if(self.checkBox[9].isChecked()):
                    self.checkBox[9].setEnabled(True)
                if(self.checkBox[10].isChecked()):
                    self.checkBox[10].setEnabled(True)
                if(self.checkBox[11].isChecked()):
                    self.checkBox[11].setEnabled(True)
                self.lineEdit[8].setEnabled(True)
                self.lineEdit[9].setEnabled(True)
                self.lineEdit[10].setEnabled(True)                
                self.lineEdit[11].setEnabled(True)
                self.uid.xm_lineEdit.setEnabled(False)
                self.uid.xm_load_pushButton.setEnabled(False)
            
            elif(self.xmDyBtn.isChecked()):
                self.uid.label_9.setEnabled(False)
                self.lineEdit[7].setEnabled(False)
                self.checkBox[7].setEnabled(False)
                self.uid.label_10.setEnabled(False)
                self.uid.label_11.setEnabled(False)
                self.uid.label_12.setEnabled(False)
                self.uid.label_13.setEnabled(False)
                self.uid.label_14.setEnabled(False)
                self.uid.velocity_label.setEnabled(False)
                self.checkBox[8].setEnabled(False)
                self.checkBox[9].setEnabled(False)
                self.checkBox[10].setEnabled(False)
                self.checkBox[11].setEnabled(False)
                self.lineEdit[8].setEnabled(False)
                self.lineEdit[9].setEnabled(False)
                self.lineEdit[10].setEnabled(False)
                self.lineEdit[11].setEnabled(False)
                self.uid.xm_lineEdit.setEnabled(False)
                self.uid.xm_load_pushButton.setEnabled(False)
            
            elif(self.xmExpBtn.isChecked()):
                self.uid.label_9.setEnabled(False)
                self.lineEdit[7].setEnabled(False)
                self.checkBox[7].setEnabled(False)
                self.uid.label_10.setEnabled(False)
                self.uid.label_11.setEnabled(False)
                self.uid.label_12.setEnabled(False)
                self.uid.label_13.setEnabled(False)
                self.uid.label_14.setEnabled(False)
                self.uid.velocity_label.setEnabled(False)
                self.checkBox[8].setEnabled(False)
                self.checkBox[9].setEnabled(False)
                self.checkBox[10].setEnabled(False)
                self.checkBox[11].setEnabled(False)
                self.lineEdit[8].setEnabled(False)
                self.lineEdit[9].setEnabled(False)
                self.lineEdit[10].setEnabled(False)
                self.lineEdit[11].setEnabled(False)
                self.uid.xm_lineEdit.setEnabled(True)
                self.uid.xm_load_pushButton.setEnabled(True)
            
    def displayValue(self):
        # display set value
        if(self.MW.ModelType == 'Motoneuron' or self.MW.ModelType == 'Motor Unit'):
            for i in range(self.isF_idx):
                item = QTableWidgetItem(str(self.isF_setValue[i]))
                self.uid.tableWidget_isF.setItem(i, 1, item)
            self.uid.lineEdit_1.setText(str(self.setTable[0]))
            self.uid.lineEdit_2.setText(str(self.setTable[1]))
            self.uid.lineEdit_3.setText(str(self.setTable[2]))
            
            iter_num = self.synF_num 
            idx = self.synF_idx 
            k = 0            
            for i in range(iter_num):
                item = QTableWidgetItem(str(self.synF_setValue[k][i-idx[k]]))
                self.syn_tableWidget_F[k].setItem(i-idx[k], 1, item)
                if(i == idx[k+1]-1):
                    k+=1
            
            iter_num = self.synT_num 
            idx = self.synT_idx
            k = 0            
            for i in range(iter_num):
                item = QTableWidgetItem(str(self.synT_setValue[k][i-idx[k]]))
                self.syn_tableWidget_T[k].setItem(i-idx[k], 1, item)
                if(i == idx[k+1]-1):
                    k+=1
        
        if(self.MW.ModelType == 'Muscle Fibers' or self.MW.ModelType == 'Motor Unit'):
            if(self.MW.ModelType == 'Muscle Fibers'):
                self.uid.lineEdit_4.setText(str(self.setTable[3]))
                self.uid.lineEdit_5.setText(str(self.setTable[4]))
                self.uid.lineEdit_6.setText(str(self.setTable[5]))
                self.uid.lineEdit_7.setText(str(self.setTable[6]))
            self.uid.lineEdit_8.setText(str(self.setTable[7]))
            self.uid.lineEdit_9.setText(str(self.setTable[8]))
            self.uid.lineEdit_10.setText(str(self.setTable[9]))
            self.uid.lineEdit_11.setText(str(self.setTable[10]))
            self.uid.lineEdit_12.setText(str(self.setTable[11]))
            self.calXmVelocity()
            
        # tab activation or deactivation
        if(self.MW.ModelType == 'Motoneuron'): 
            self.uid.tabWidget.setTabEnabled(0, True)
            self.uid.tabWidget.setTabEnabled(1, True) 
            self.uid.tabWidget.setTabEnabled(2, False)
            self.uid.tabWidget.setTabEnabled(3, False)
            self.uid.tabWidget.setCurrentIndex(0) 
        elif(self.MW.ModelType == 'Muscle Fibers'):
            self.uid.tabWidget.setTabEnabled(0, False)
            self.uid.tabWidget.setTabEnabled(1, False) 
            self.uid.tabWidget.setTabEnabled(2, True)
            self.uid.tabWidget.setTabEnabled(3, True)
            self.uid.tabWidget.setCurrentIndex(2)
        elif(self.MW.ModelType == 'Motor Unit'):
            self.uid.tabWidget.setTabEnabled(0, True)
            self.uid.tabWidget.setTabEnabled(1, True) 
            self.uid.tabWidget.setTabEnabled(2, False)
            self.uid.tabWidget.setTabEnabled(3, True)
            self.uid.tabWidget.setCurrentIndex(0)
            
        # radio button set
        self.isSetButton.setChecked(True)
        for i in range(len(self.syn_setButton)):
            self.syn_setButton[i].setChecked(True)
        self.spSetButton.setChecked(True)
        self.xmSetButton.setChecked(True)
        self.sesyn_noise.setChecked(self.sesyn_setNoise)
        self.sisyn_noise.setChecked(self.sisyn_setNoise)
        self.desyn_noise.setChecked(self.desyn_setNoise)
        self.disyn_noise.setChecked(self.disyn_setNoise)
        self.setObjectEnabled()
        self.checkValue()
        
        # ButtonBox enable
        self.uid.applyButton.setEnabled(False)
        self.uid.okButton.setEnabled(False)

    def calXmVelocity(self):
        # calculate slope from isokinetic Xm according to changing the value
        try:
            t1 = float(self.uid.lineEdit_11.text())
            t1 = self.setTable[10]
            t2 = float(self.uid.lineEdit_12.text())
            t2 = self.setTable[11]
            x1 = float(self.uid.lineEdit_9.text())
            x1 = self.setTable[8]
            x2 = float(self.uid.lineEdit_10.text())
            x2 = self.setTable[9]
            ms = 0.001
            v = (x2-x1)/((t2-t1)*ms)
            self.uid.velocity_label.setText(str(v))    
        
        except:
            self.uid.velocity_label.setText('')
        
    def setValue(self):
        ## set the value
        if(self.MW.ModelType == 'Motoneuron' or self.MW.ModelType == 'Motor Unit'):
            # Isoma
            if(self.gen_isType == 'Step'):
                self.isF_setValue = self.gen_isF_heav[:]
                self.isSetButton = self.isFixBtn
            elif(self.gen_isType == 'Ramp'):   
                self.setTable[0] = self.gen_is_iValue
                self.setTable[1] = self.gen_is_pValue
                self.setTable[2] = self.gen_is_period
                self.isSetButton = self.isTriBtn
            elif(self.gen_isType == 'Import'):
                self.isSetButton = self.isImtBtn
            
            # Isyn
            # noise checkbox set
            if(self.sesyn_noise.isChecked()):
                self.sesyn_setNoise = True
            else:
                self.sesyn_setNoise = False
            if(self.sisyn_noise.isChecked()):
                self.sisyn_setNoise = True
            else:
                self.sisyn_setNoise = False
            if(self.desyn_noise.isChecked()):
                self.desyn_setNoise = True
            else:
                self.desyn_setNoise = False
            if(self.disyn_noise.isChecked()):
                self.disyn_setNoise = True
            else:
                self.disyn_setNoise = False
            
            for i in range(len(self.gen_synType)):
                if(self.gen_synType[i] == 'Step'):
                    new_arr = [self.gen_synF_tau[i], self.gen_synF_std_max[i]]
                    com_arr = self.gen_synF_heav[i] + new_arr
                    self.synF_setValue[i] = com_arr[:]
                    self.syn_setButton[i] = self.syn_FixBtn[i]
                
                elif(self.gen_synType[i] == 'Ramp'):
                    self.synT_setValue[i] = [self.gen_synT_iValue[i],self.gen_synT_pValue[i],self.gen_synT_period[i],self.gen_synT_tau[i],self.gen_synT_std_max[i]]
                    self.syn_setButton[i] = self.syn_TriBtn[i]
                
                elif(self.gen_synType[i] == 'Import'):
                    self.syn_setButton[i] = self.syn_ImtBtn[i]

        if(self.MW.ModelType == 'Muscle Fibers' or self.MW.ModelType == 'Motor Unit'): 
            # Iaxon
            if(self.MW.ModelType == 'Muscle Fibers'):
                if(self.gen_spikeType == 'User'):
                    self.setTable[3] = self.gen_spike_t1
                    self.setTable[4] = self.gen_spike_t2
                    self.setTable[5] = self.gen_spike_hz
                    self.setTable[6] = self.gen_spike_scale
                    self.spSetButton = self.spUserBtn
                
                elif(self.gen_spikeType == 'Exp'):
                    self.spSetButton = self.spExpBtn
            # Xm
            if(self.gen_xmType == 'Isometric'):
                self.setTable[7] = self.gen_xm_value
                self.xmSetButton = self.xmIsomeBtn
            
            elif(self.gen_xmType == 'Isokinetic'):
                self.xmSetButton = self.xmIsokiBtn
                self.setTable[8] = self.gen_xm_value1
                self.setTable[9] = self.gen_xm_value2
                self.setTable[10] = self.gen_xm_t1
                self.setTable[11] = self.gen_xm_t2
            
            elif(self.gen_xmType == 'Exp'):
                self.xmSetButton = self.xmExpBtn
            
            elif(self.gen_xmType == 'Dynamic'):
                self.xmSetButton = self.xmDyBtn                
        
        self.t_start, self.t_final, self.t_dt, self.t_pt = self.MW.ISW.t_setTable
        
        # Text message
        if(self.MW.ModelType == 'Motoneuron'):
            self.MW.MN.setInputSignal(self.ISG.signalType, self.ISG.iValue, self.ISG.pValue, self.ISG.Is_0, self.ISG.heav_param, self.ISG.period, self.ISG.times, self.ISG.Is)
            self.MW.MN.setSynConSignal(self.se_ch, self.si_ch, self.de_ch, self.di_ch, self.s_SCSG.e_times, self.s_SCSG.i_times, self.d_SCSG.e_times, self.d_SCSG.i_times, self.s_SCSG.G_e, self.s_SCSG.G_i, self.d_SCSG.G_e, self.d_SCSG.G_i)
            self.MW.setTextEdit("[Motoneuron] Input signals set.")
        elif(self.MW.ModelType == 'Muscle Fibers'):
            self.MW.MF.setSpikeSignal(self.SSG.spike, self.SSG.spike_idx, self.SSG.SpikeTimes)                
            self.MW.MF.setXmSignal(self.XSG.signalType, self.XSG.times, self.XSG.xm)
            self.MW.setTextEdit("[Muscle Fibers] Input signals set.")
        elif(self.MW.ModelType == 'Motor Unit'):
            self.MW.MU.setInputSignal(self.ISG.signalType, self.ISG.iValue, self.ISG.pValue, self.ISG.Is_0, self.ISG.heav_param, self.ISG.period, self.ISG.times, self.ISG.Is)
            self.MW.MU.setSynConSignal(self.se_ch, self.si_ch, self.de_ch, self.di_ch, self.s_SCSG.e_times, self.s_SCSG.i_times, self.d_SCSG.e_times, self.d_SCSG.i_times, self.s_SCSG.G_e, self.s_SCSG.G_i, self.d_SCSG.G_e, self.d_SCSG.G_i)
            self.MW.MU.setXmSignal(self.XSG.signalType, self.XSG.times, self.XSG.xm)
            self.MW.setTextEdit("[Motor Unit] Input signals set.")

        # save tab index
        currentWidget = self.uid.tabWidget.currentWidget()
        if(currentWidget == self.uid.is_tab):
            i = 0
        elif(currentWidget == self.uid.syn_tab):
            i = 1
        elif(currentWidget == self.uid.sp_tab):
            i = 2
        elif(currentWidget == self.uid.xm_tab):
            i = 3

        self.displayValue()
        self.uid.tabWidget.setCurrentIndex(i)

    def openFile(self):
        curPath = os.path.dirname( os.path.abspath( sys.argv[0] ) )
        pushButton = self.sender()
        
        try:
            if(pushButton == self.uid.is_load_pushButton):
                signal = 'Isoma'
            elif(pushButton == self.uid.sesyn_load_pushButton):
                signal = 'G_esyn_soma'
            elif(pushButton == self.uid.sisyn_load_pushButton):
                signal = 'G_isyn_soma'
            elif(pushButton == self.uid.desyn_load_pushButton):
                signal = 'G_esyn_dend'
            elif(pushButton == self.uid.disyn_load_pushButton):
                signal = 'G_isyn_dend'
            elif(pushButton == self.uid.sp_load_pushButton):
                signal = 'Iaxon'
            elif(pushButton == self.uid.xm_load_pushButton):
                signal = 'Xm'
            fileName = QFileDialog.getOpenFileName(self, "Open a " + signal + " data File (Must be a csv format)", curPath+"/parameters", "CSV Files (*.csv)")
        
        except:
            self.MW.setTextEdit("[Error] Failed to open " + signal + "data file.")
            self.MW.raise_()
            self.MW.activateWindow()
        
        if(fileName != ''):
            self.importData(fileName, pushButton, signal)

    def importData(self, fileName, pushButton, signal):
        currentWidget = self.uid.tabWidget.currentWidget()
        
        try:
            data = pd.read_csv(unicode(fileName))
            import_time = np.zeros(len(data.index))
            import_data = np.zeros(len(data.index))
            
            # Isoma
            if(pushButton == self.uid.is_load_pushButton):                       
                for i in range(len(data.index)):
                    import_time[i] = float(data.time[i])
                    import_data[i] = float(data.Is[i])
                    if(math.isnan(import_time[i]) and math.isnan(import_data[i])):
                        raise Exception
                        
                self.expIs_time = import_time
                self.expIs = import_data
                self.uid.is_lineEdit.setText(fileName)
            
            # Isyn
            elif(currentWidget == self.uid.syn_tab):
                for i in range(len(self.syn_load_pushButton)):
                    if(pushButton == self.syn_load_pushButton[i]):
                        if(i==0):
                            G = data.G_esyn_soma
                        elif(i==1):
                            G = data.G_isyn_soma
                        elif(i==2):
                            G = data.G_esyn_dend
                        elif(i==3):
                            G = data.G_isyn_dend
                        
                        for k in range(len(data.index)):
                            import_time[k] = float(data.time[k])
                            import_data[k] = float(G[k])
                            if(math.isnan(import_time[k]) and math.isnan(import_data[k])):
                                raise Exception
                        self.exp_syn_time[i] = import_time
                        self.exp_syn[i] = import_data
                        self.syn_lineEdit[i].setText(fileName)
            
            # Iaxon
            elif(pushButton == self.uid.sp_load_pushButton):
                for i in range(len(data.index)):
                    import_data[i] = float(data.time[i])
                self.ExpSpike = import_data
                self.uid.sp_lineEdit.setText(fileName)
            
            # Xm
            elif(pushButton == self.uid.xm_load_pushButton):                        
                for i in range(len(data.index)):
                    import_time[i] = float(data.time[i])
                    import_data[i] = float(data.Xm[i])
                    if(math.isnan(import_time[i]) and math.isnan(import_data[i])):
                        raise Exception
                
                self.expXm_time = import_time
                self.expXm = import_data
                self.uid.xm_lineEdit.setText(fileName) 
        
        except:
            self.MW.setTextEdit("[Error] Failed to load " + signal + " data.")
            self.MW.raise_()
            self.MW.activateWindow()
            return

    def genSignal(self):
        # generate signal
        self.t_start, self.t_final, self.t_dt, self.t_pt = self.MW.ISW.t_setTable
        num_steps = int(np.floor((self.t_final-self.t_start)/self.t_dt)+1) 
        times = np.linspace(self.t_start, self.t_final, num_steps)
        currentWidget = self.uid.tabWidget.currentWidget()
        
        # Isoma
        if(currentWidget == self.uid.is_tab):
            if(self.isFixBtn.isChecked()):
                signalType = 'Step'
                for i in range(len(self.gen_isF_heav)):
                    item = self.uid.tableWidget_isF.item(i, 1)
                    self.gen_isF_heav[i] = float(item.text())
                self.ISG.setValue(signalType, self.t_final, self.t_dt, heav_param=self.gen_isF_heav)
            
            elif(self.isTriBtn.isChecked()): 
                signalType = 'Ramp'
                self.gen_is_iValue = float(self.uid.lineEdit_1.text())
                self.gen_is_pValue = float(self.uid.lineEdit_2.text())
                self.gen_is_period = float(self.uid.lineEdit_3.text())
                self.ISG.setValue(signalType, self.t_final, self.t_dt, iValue=self.gen_is_iValue, pValue=self.gen_is_pValue, period=self.gen_is_period)
            
            elif(self.isImtBtn.isChecked()): 
                signalType = 'Import'
                self.ISG.setValue(signalType, self.t_final, self.t_dt, 0, 0, 0, 0, self.expIs_time, self.expIs)
            
            try:
                self.ISG.genSignal()
                self.gen_isType = signalType
            
            except:
                self.MW.setTextEdit("[Error] Failed to generate Isoma signal.")
                self.MW.raise_()
                self.MW.activateWindow()
                # ButtonBox enable
                self.uid.applyButton.setEnabled(False)
                self.uid.okButton.setEnabled(False)
                return
            
            # x, y array
            x_array = self.ISG.times
            y_array = self.ISG.Is

        # Isyn
        elif(currentWidget == self.uid.syn_tab):  
            for i in range(len(self.groupBox_syn)):
                if(self.groupBox_syn[i].isEnabled() == True):
                    if(i==0):
                        synType = 'Excitatory'
                        signal = 'G_esyn_soma'
                    elif(i==1):
                        synType = 'Inhibitory'
                        signal = 'G_isyn_soma'
                    elif(i==2):
                        synType = 'Excitatory'
                        signal = 'G_esyn_dend'
                    elif(i==3):
                        synType = 'Inhibitory'
                        signal = 'G_isyn_dend'
                    
                    if(self.syn_FixBtn[i].isChecked()):
                        signalType = 'Step'
                        for k in range(len(self.gen_synF_heav[0])):
                            item = self.syn_tableWidget_F[i].item(k, 1)
                            self.gen_synF_heav[i][k] = float(item.text())
                        heav = self.gen_synF_heav[i]
                        iValue = 0.
                        pValue = 0.
                        period = 0.
                        exp_time = 0.
                        exp_G = 0.
                        
                        # noise O
                        if(self.syn_noise[i].isChecked()):
                            noise = True
                            item = self.syn_tableWidget_F[i].item(k+1, 1)
                            self.gen_synF_tau[i] = float(item.text())
                            item = self.syn_tableWidget_F[i].item(k+2, 1)
                            self.gen_synF_std_max[i] = float(item.text())
                            tau=self.gen_synF_tau[i]
                            std_max=self.gen_synF_std_max[i]
                        
                        # noise X
                        else:
                            noise = False
                            tau=0.
                            std_max=0.                
                    
                    elif(self.syn_TriBtn[i].isChecked()): 
                        signalType = 'Ramp'   
                        item = self.syn_tableWidget_T[i].item(0, 1)
                        self.gen_synT_iValue[i] = float(item.text())
                        item = self.syn_tableWidget_T[i].item(1, 1)
                        self.gen_synT_pValue[i] = float(item.text())
                        item = self.syn_tableWidget_T[i].item(2, 1)
                        self.gen_synT_period[i] = float(item.text())
                        heav = []
                        iValue = self.gen_synT_iValue[i]
                        pValue = self.gen_synT_pValue[i]
                        period = self.gen_synT_period[i]
                        exp_time = 0.
                        exp_G = 0.
                        
                        if(self.syn_noise[i].isChecked()):
                            noise = True
                            item = self.syn_tableWidget_T[i].item(3, 1)
                            self.gen_synT_tau[i] = float(item.text())
                            item = self.syn_tableWidget_T[i].item(4, 1)
                            self.gen_synT_std_max[i] = float(item.text())
                            tau=self.gen_synT_tau[i]
                            std_max=self.gen_synT_std_max[i]
                        
                        else:
                            noise = False
                            tau=0.
                            std_max=0.
                    
                    elif(self.syn_ImtBtn[i].isChecked()):
                        signalType = 'Import'
                        heav = []
                        iValue = 0.
                        pValue = 0.
                        period = 0.
                        noise = False
                        tau=0.
                        std_max=0.
                        exp_time = self.exp_syn_time[i]
                        exp_G = self.exp_syn[i]
                    
                    try:
                        if(i < 2):
                            self.s_SCSG.setValue(synType, signalType, self.t_final, self.t_dt, heav, iValue, pValue, period, tau, std_max, noise, exp_time, exp_G)
                            self.s_SCSG.genSignal()
                        
                        else:
                            self.d_SCSG.setValue(synType, signalType, self.t_final, self.t_dt, heav, iValue, pValue, period, tau, std_max, noise, exp_time, exp_G)
                            self.d_SCSG.genSignal()
                        self.gen_synType[i] = signalType
                    
                    except:
                        self.MW.setTextEdit("[Error] Failed to generate " + signal + " signal.")
                        self.MW.raise_()
                        self.MW.activateWindow()
                        # ButtonBox enable
                        self.uid.applyButton.setEnabled(False)
                        self.uid.okButton.setEnabled(False)
                        return
                
                # no synapse channel
                else:
                    pass
                
            # x, y array
            Se_x_array = self.s_SCSG.e_times
            Si_x_array = self.s_SCSG.i_times
            Se_y_array = self.s_SCSG.G_e
            Si_y_array = self.s_SCSG.G_i
            De_x_array = self.d_SCSG.e_times
            Di_x_array = self.d_SCSG.i_times
            De_y_array = self.d_SCSG.G_e
            Di_y_array = self.d_SCSG.G_i
        
        # Iaxon
        elif(currentWidget == self.uid.sp_tab):            
            if(self.spUserBtn.isChecked()):
                signalType = 'User'
                self.gen_spike_t1 = float(self.uid.lineEdit_4.text())                       
                self.gen_spike_t2 = float(self.uid.lineEdit_5.text())
                self.gen_spike_hz = float(self.uid.lineEdit_6.text())
                self.gen_spike_scale = float(self.uid.lineEdit_7.text())
                self.SSG.setValue(signalType, self.t_final, self.t_dt, self.gen_spike_t1, self.gen_spike_t2, self.gen_spike_hz, self.gen_spike_scale)
            
            elif(self.spExpBtn.isChecked()):
                signalType = 'Exp'
                self.SSG.setValue(signalType, self.t_final, self.t_dt, 0,0,0,0, self.ExpSpike)
            
            try:
                self.SSG.genSignal()
                self.gen_spikeType = signalType
                self.gen_spike = True
            
            except:
                self.MW.setTextEdit("[Error] Failed to generate Iaxon signal.")
                self.MW.raise_()
                self.MW.activateWindow()
                # ButtonBox enable
                self.uid.applyButton.setEnabled(False)
                self.uid.okButton.setEnabled(False)
                return
            
            # x, y array
            x_array = times
            y_array = self.SSG.spike
                
        # Xm
        elif(currentWidget == self.uid.xm_tab):  
            if(self.xmIsomeBtn.isChecked()):
                signalType = 'Isometric'
                self.gen_xm_value = float(self.uid.lineEdit_8.text())
                self.XSG.setValue(signalType, self.t_final, self.t_dt, self.gen_xm_value)
            
            elif(self.xmIsokiBtn.isChecked()):
                signalType = 'Isokinetic'
                self.gen_xm_value1 = float(self.uid.lineEdit_9.text())                       
                self.gen_xm_value2 = float(self.uid.lineEdit_10.text())                      
                self.gen_xm_t1 = float(self.uid.lineEdit_11.text())                       
                self.gen_xm_t2 = float(self.uid.lineEdit_12.text())
                self.XSG.setValue(signalType, self.t_final, self.t_dt, self.gen_xm_value1, self.gen_xm_value2, self.gen_xm_t1, self.gen_xm_t2)
            
            elif(self.xmExpBtn.isChecked()):
                signalType = 'Exp'
                self.XSG.setValue(signalType, self.t_final, self.t_dt, 0,0,0,0, self.expXm_time, self.expXm)
            
            elif(self.xmDyBtn.isChecked()):
                signalType = 'Dynamic'
                self.XSG.setValue(signalType, self.t_final, self.t_dt)
            

            self.XSG.genSignal()
            self.gen_xmType = signalType
            self.gen_xm = True

            
            # x, y array     
            x_array = self.XSG.times
            y_array = self.XSG.xm

        ## plotting the signal   
        plt.rc('figure', figsize=(6, 5))
        ModelType = self.MW.ModelType
        
        # Isoma
        if(currentWidget == self.uid.is_tab):
            # figure creation
            if(self.is_fig == None):
                self.is_fig, self.is_ax=plt.subplots(nrows=1, ncols=1)
                self.is_fig.canvas.set_window_title('Intracellular Stimulation')
                self.is_fig.canvas.mpl_connect('close_event', self.FigureCloseEvent)
            self.is_ax.clear() 
            self.is_ax.set_xlabel('Time (ms)') 
            self.is_ax.set_ylabel('Is (nA)')
            self.is_ax.grid()
            self.is_ax.plot(x_array, y_array, 'r') 
            self.is_fig.canvas.draw()
            self.is_fig.show()
            self.is_fig.canvas.manager.window.raise_()
            self.MW.setTextEdit("["+ModelType+"] Somatic input signal (Isoma) generated.")
    
        # Isyn  
        elif(currentWidget == self.uid.syn_tab):
            # figure creation
            if(self.syn_fig == None):
                self.syn_fig = plt.figure()
                self.syn_fig.canvas.set_window_title('Synaptic Stimulation')
                self.syn_fig.canvas.mpl_connect('close_event', self.FigureCloseEvent)
            
            else:
                self.syn_fig.clear()
            
            s_logic = self.uid.groupBox_sesyn.isEnabled() | self.uid.groupBox_sisyn.isEnabled()
            d_logic = self.uid.groupBox_desyn.isEnabled() | self.uid.groupBox_disyn.isEnabled()
            
            if(s_logic & d_logic): # 2 axes
                self.s_ax = self.syn_fig.add_subplot(2, 1, 1)
                self.d_ax = self.syn_fig.add_subplot(2, 1, 2)
            
            elif(s_logic ^ d_logic):# 1 axis
                if(s_logic == True):
                    self.s_ax = self.syn_fig.add_subplot(1, 1, 1)
                
                else:
                    self.d_ax = self.syn_fig.add_subplot(1, 1, 1)
            
            elif(s_logic & d_logic == False):
                self.MW.setTextEdit("["+ModelType+"] No synaptic input included.")
                # no data display
                self.syn_fig.canvas.draw()
                self.syn_fig.show()
                self.MW.raise_()
                self.MW.activateWindow()
                return
            
            if(self.uid.groupBox_sesyn.isEnabled() == True):
                self.s_ax.plot(Se_x_array, Se_y_array, 'r', label='G_esyn_soma')
            if(self.uid.groupBox_sisyn.isEnabled() == True):
                self.s_ax.plot(Si_x_array, Si_y_array, 'b', label='G_isyn_soma')
            if(self.uid.groupBox_desyn.isEnabled() == True):
                self.d_ax.plot(De_x_array, De_y_array, 'r', label='G_esyn_dend')
            if(self.uid.groupBox_disyn.isEnabled() == True):
                self.d_ax.plot(Di_x_array, Di_y_array, 'b', label='G_isyn_dend')
            
            if(s_logic == True):
                self.s_ax.set_xlabel('Time (ms)')
                self.s_ax.set_ylabel('G (mS/cm^2)')
                self.s_ax.grid()
                self.s_ax.legend(loc='best')
            
            if(d_logic == True):
                self.d_ax.set_xlabel('Time (ms)') 
                self.d_ax.set_ylabel('G (mS/cm^2)')
                self.d_ax.grid()
                self.d_ax.legend(loc='best')
            self.syn_fig.canvas.draw()
            self.syn_fig.show()
            self.syn_fig.canvas.manager.window.raise_()
            self.MW.setTextEdit("["+ModelType+"] Synaptic input signal (Isyn) generated.")
        
        # Iaxon
        elif(currentWidget == self.uid.sp_tab):
            ### figure creation
            if(self.sp_fig == None):
                self.sp_fig, self.sp_ax=plt.subplots(nrows=1, ncols=1)
                self.sp_fig.canvas.set_window_title('Axonal Stimulation')  
                self.sp_fig.canvas.mpl_connect('close_event', self.FigureCloseEvent)
            self.sp_ax.clear() 
            self.sp_ax.set_xlabel('Time (ms)') 
            self.sp_ax.set_ylabel('Iaxon')
            self.sp_ax.grid()
            self.sp_ax.plot(x_array, y_array, 'r') 
            self.sp_fig.canvas.draw()
            self.sp_fig.show()
            self.sp_fig.canvas.manager.window.raise_()
            self.MW.setTextEdit("["+ModelType+"] Axonal input signal (Iaxon) generated.")
               
        # Xm
        elif(currentWidget == self.uid.xm_tab):
            ### figure creation
            if(self.xm_fig == None):
                self.xm_fig, self.xm_ax=plt.subplots(nrows=1, ncols=1)
                self.xm_fig.canvas.set_window_title('Muscle Length')
                self.xm_fig.canvas.mpl_connect('close_event', self.FigureCloseEvent)
            self.xm_ax.clear() 
            self.xm_ax.set_xlabel('Time (ms)') 
            self.xm_ax.set_ylabel('Xm (mm)')
            self.xm_ax.grid()
            self.xm_ax.plot(x_array, y_array, 'r') 
            self.xm_fig.canvas.draw()
            self.xm_fig.show()
            self.xm_fig.canvas.manager.window.raise_()
            self.MW.setTextEdit("["+ModelType+"] Muscle length signal (Xm) generated.")
        
        # ButtonBox enable
        self.uid.applyButton.setEnabled(True)
        self.uid.okButton.setEnabled(True)
        
    def FigureCloseEvent(self, event):
        figure = event.canvas.figure
        
        if(figure == self.is_fig):
            self.is_fig = None
        elif(figure == self.syn_fig):
            self.syn_fig = None
        elif(figure == self.sp_fig):
            self.sp_fig = None
        elif(figure == self.xm_fig):
            self.xm_fig = None
        
# Output Signals Window class         
class OscilloscopeWindow(QDialog):
    def __init__(self, MW):
        super(OscilloscopeWindow, self).__init__()
        self.MW = MW
        
        # subscript, superscript, Greek letter with unicode and unit label
        super_2 = u'\u00B2'
        super_2_p = super_2 + u'\u207A'
        mV = ' (mV)'
        mM = ' (mM)'
        cm_2 = 'cm'+super_2
        mA_cm_2 = ' (mA/' + cm_2 + ')'
        mS_cm_2 = ' (mS/' + cm_2 + ')'
        M = ' (M)'
        sub_i = u'\u1D62'
        tilde = u'\u0303'
        
        # subscript with HTML
        sub_Naf = '<sub>Naf</sub>' 
        sub_Nap = '<sub>Nap</sub>'
        sub_Kdr = '<sub>Kdr</sub>' 
        sub_KCa = '<sub>K(Ca)</sub>' 
        sub_Can = '<sub>Can</sub>' 
        sub_H = '<sub>H</sub>' 
        sub_esyn = '<sub>esyn</sub>' 
        sub_isyn = '<sub>isyn</sub>'
        sub_Cal = '<sub>Cal</sub>'
        sub_Ca = '<sub>Ca</sub>' 
        sub_SR = '<sub>SR</sub>'
        sub_SP = '<sub>SP</sub>' 
        sub_CE = '<sub>CE</sub>'
        
        ## initialization
        # set value
        self.mn_itemList = ['Is', 'V_soma', 'Firing_rate', '[Ca]_soma', 'E_Ca_soma', 'I_Naf_soma', 'm_Naf_soma', 'h_Naf_soma', 'I_Nap_soma', 'm_Nap_soma', 'I_Kdr_soma', 'n_Kdr_soma', 'I_Kca_soma', 'I_Can_soma', 'm_Can_soma', 'h_Can_soma', 'I_H_soma', 'm_H_soma', 'I_esyn_soma', 'G_esyn_soma', 'I_isyn_soma', 'G_isyn_soma', 'V_dend', '[Ca]_dend', 'E_Ca_dend', 'I_Cal_dend', 'l_Cal_dend', 'I_Naf_dend', 'm_Naf_dend', 'h_Naf_dend', 'I_Nap_dend', 'm_Nap_dend', 'I_Kdr_dend', 'n_Kdr_dend', 'I_Kca_dend', 'I_Can_dend', 'm_Can_dend', 'h_Can_dend', 'I_H_dend', 'm_H_dend', 'I_esyn_dend', 'G_esyn_dend', 'I_isyn_dend', 'G_isyn_dend']
        self.mf_itemList = ['F', 'Spike', 'R', 'Xm', 'Vm', 'Am', 'Cs', 'CaSR', 'CaSRCS', 'B', 'CaSP', 'T' ,'CaSPB', 'CaSPT', 'A_tilde', 'XCE', 'A']
        self.mn_list_num = len(self.mn_itemList)
        self.mf_list_num = len(self.mf_itemList)
        self.mn_ScopeList = set()
        self.mf_ScopeList = set()
        
        # check box state
        self.mn_checkbox_state = []
        self.mf_checkbox_state = []
        for i in range(self.mn_list_num):
            self.mn_checkbox_state.append(False) # unchecked
        for i in range(self.mf_list_num):
            self.mf_checkbox_state.append(False)
        
        ## GUI
        self.uid = Ui_OD()
        self.uid.setupUi(self)
        
        # groupbox Label
        self.uid.checkBox_Sv.setText('V' + mV)
        self.uid.groupBox_Sca.setTitle('d[Ca' + super_2_p + ']' + sub_i + '/dt  ')
        self.uid.checkBox_Sca.setText('[Ca' + super_2_p + ']' + sub_i + mM)
        self.uid.label_Svca.setText('E' + sub_Ca + mV)
        self.uid.checkBox_SnaI.setText('I' + mA_cm_2)
        self.uid.label_Snaf.setText('I' + sub_Naf)
        self.uid.checkBox_SnapI.setText('I' + mA_cm_2)
        self.uid.label_Snap.setText('I' + sub_Nap)
        self.uid.checkBox_SkdrI.setText('I' + mA_cm_2)
        self.uid.label_Skdr.setText('I' + sub_Kdr)
        self.uid.checkBox_SkcaI.setText('I' + mA_cm_2)
        self.uid.label_Skca.setText('I' + sub_KCa)
        self.uid.checkBox_ScaI.setText('I' + mA_cm_2)
        self.uid.label_Scan.setText('I' + sub_Can)
        self.uid.checkBox_ShI.setText('I' + mA_cm_2)
        self.uid.label_Sh.setText('I' + sub_H)
        self.uid.checkBox_SesynI.setText('I ' + mA_cm_2)
        self.uid.label_Sesyn.setText('I' + sub_esyn)
        self.uid.checkBox_SisynI.setText('I ' + mA_cm_2)
        self.uid.label_Sisyn.setText('I' + sub_isyn)
        self.uid.checkBox_Sgesyn.setText('G' + mS_cm_2)
        self.uid.checkBox_Sgisyn.setText('G' + mS_cm_2)
        self.uid.checkBox_Dv.setText('V' + mV)
        self.uid.groupBox_Dca.setTitle('d[Ca' + super_2_p + ']' + sub_i + '/dt  ')
        self.uid.checkBox_Dca.setText('[Ca' + super_2_p + ']' + sub_i + mM)
        self.uid.label_Dvca.setText('E' + sub_Ca + mV)
        self.uid.checkBox_DcalI.setText('I' + mA_cm_2)
        self.uid.label_Dcal.setText('I' + sub_Cal)
        self.uid.checkBox_DnaI.setText('I' + mA_cm_2)
        self.uid.label_Dnaf.setText('I' + sub_Naf)
        self.uid.checkBox_DnapI.setText('I' + mA_cm_2)
        self.uid.label_Dnap.setText('I' + sub_Nap)
        self.uid.checkBox_DkdrI.setText('I' + mA_cm_2)
        self.uid.label_Dkdr.setText('I' + sub_Kdr)
        self.uid.checkBox_DkcaI.setText('I' + mA_cm_2)
        self.uid.label_Dkca.setText('I' + sub_KCa)
        self.uid.checkBox_DcaI.setText('I' + mA_cm_2)
        self.uid.label_Dcan.setText('I' + sub_Can)
        self.uid.checkBox_DhI.setText('I' + mA_cm_2)
        self.uid.label_Dh.setText('I' + sub_H)
        self.uid.checkBox_DesynI.setText('I ' + mA_cm_2)
        self.uid.label_Desyn.setText('I' + sub_esyn)
        self.uid.checkBox_DisynI.setText('I ' + mA_cm_2)
        self.uid.label_Disyn.setText('I' + sub_isyn)
        self.uid.checkBox_Dgesyn.setText('G' + mS_cm_2)
        self.uid.checkBox_Dgisyn.setText('G' + mS_cm_2)
        self.uid.checkBox_CS.setText('[CS]' + M)
        self.uid.label_CaSR.setText('[Ca' + sub_SR + ']' + M)
        self.uid.label_CaSRCS.setText('[Ca' + sub_SR + 'CS]' + M)
        self.uid.checkBox_B.setText('[B]' + M)
        self.uid.label_CaSP.setText('[Ca' + sub_SP + ']' + M)
        self.uid.checkBox_T.setText('[T]' + M)
        self.uid.label_CaSPB.setText('[Ca' + sub_SP + 'B]' + M)
        self.uid.label_CaSPT.setText('[Ca' + sub_SP + 'T]' + M)
        self.uid.checkBox_Xm.setText('Xm (mm)')
        self.uid.checkBox_Vm.setText('Vm (mm/s)')
        self.uid.checkBox_Am.setText('Am (mm/s' + super_2 + ')')
        self.uid.checkBox_F.setText('F (N)')
        self.uid.label_A.setText('A' + tilde)
        self.uid.label_XCE.setText('X' + sub_CE + ' (mm)')
        
        # set radio button
        self.displayResult = self.uid.radioButton_Individual
        self.displayResult.setChecked(True)
        
        # checkbox_list
        self.mn_checkbox_list = []
        self.mn_checkbox_list.append(self.uid.checkBox_SIs)
        self.mn_checkbox_list.append(self.uid.checkBox_Sv)
        self.mn_checkbox_list.append(self.uid.checkBox_SR)
        self.mn_checkbox_list.append(self.uid.checkBox_Sca)
        self.mn_checkbox_list.append(self.uid.checkBox_Svca)
        self.mn_checkbox_list.append(self.uid.checkBox_SnaI)
        self.mn_checkbox_list.append(self.uid.checkBox_Snam)
        self.mn_checkbox_list.append(self.uid.checkBox_Snah)
        self.mn_checkbox_list.append(self.uid.checkBox_SnapI)
        self.mn_checkbox_list.append(self.uid.checkBox_Snapm)
        self.mn_checkbox_list.append(self.uid.checkBox_SkdrI)
        self.mn_checkbox_list.append(self.uid.checkBox_Skdr)
        self.mn_checkbox_list.append(self.uid.checkBox_SkcaI)
        self.mn_checkbox_list.append(self.uid.checkBox_ScaI)
        self.mn_checkbox_list.append(self.uid.checkBox_Scam)
        self.mn_checkbox_list.append(self.uid.checkBox_Scah)
        self.mn_checkbox_list.append(self.uid.checkBox_ShI)
        self.mn_checkbox_list.append(self.uid.checkBox_Shm)
        self.mn_checkbox_list.append(self.uid.checkBox_SesynI)
        self.mn_checkbox_list.append(self.uid.checkBox_Sgesyn)
        self.mn_checkbox_list.append(self.uid.checkBox_SisynI)
        self.mn_checkbox_list.append(self.uid.checkBox_Sgisyn)
        self.mn_checkbox_list.append(self.uid.checkBox_Dv)
        self.mn_checkbox_list.append(self.uid.checkBox_Dca)
        self.mn_checkbox_list.append(self.uid.checkBox_Dvca)
        self.mn_checkbox_list.append(self.uid.checkBox_DcalI)
        self.mn_checkbox_list.append(self.uid.checkBox_Dcal)
        self.mn_checkbox_list.append(self.uid.checkBox_DnaI)
        self.mn_checkbox_list.append(self.uid.checkBox_Dnam)
        self.mn_checkbox_list.append(self.uid.checkBox_Dnah)
        self.mn_checkbox_list.append(self.uid.checkBox_DnapI)
        self.mn_checkbox_list.append(self.uid.checkBox_Dnapm)
        self.mn_checkbox_list.append(self.uid.checkBox_DkdrI)
        self.mn_checkbox_list.append(self.uid.checkBox_Dkdr)
        self.mn_checkbox_list.append(self.uid.checkBox_DkcaI)
        self.mn_checkbox_list.append(self.uid.checkBox_DcaI)
        self.mn_checkbox_list.append(self.uid.checkBox_Dcam)
        self.mn_checkbox_list.append(self.uid.checkBox_Dcah)
        self.mn_checkbox_list.append(self.uid.checkBox_DhI)
        self.mn_checkbox_list.append(self.uid.checkBox_Dhm)
        self.mn_checkbox_list.append(self.uid.checkBox_DesynI)
        self.mn_checkbox_list.append(self.uid.checkBox_Dgesyn)
        self.mn_checkbox_list.append(self.uid.checkBox_DisynI)
        self.mn_checkbox_list.append(self.uid.checkBox_Dgisyn)
        self.mf_checkbox_list = []
        self.mf_checkbox_list.append(self.uid.checkBox_F)
        self.mf_checkbox_list.append(self.uid.checkBox_Spike)
        self.mf_checkbox_list.append(self.uid.checkBox_R)
        self.mf_checkbox_list.append(self.uid.checkBox_Xm)
        self.mf_checkbox_list.append(self.uid.checkBox_Vm)
        self.mf_checkbox_list.append(self.uid.checkBox_Am)
        self.mf_checkbox_list.append(self.uid.checkBox_CS)
        self.mf_checkbox_list.append(self.uid.checkBox_CaSR)
        self.mf_checkbox_list.append(self.uid.checkBox_CaSRCS)
        self.mf_checkbox_list.append(self.uid.checkBox_B)
        self.mf_checkbox_list.append(self.uid.checkBox_CaSP)
        self.mf_checkbox_list.append(self.uid.checkBox_T)
        self.mf_checkbox_list.append(self.uid.checkBox_CaSPB)
        self.mf_checkbox_list.append(self.uid.checkBox_CaSPT)
        self.mf_checkbox_list.append(self.uid.checkBox_A)
        self.mf_checkbox_list.append(self.uid.checkBox_XCE)
        self.mf_checkbox_list.append(self.uid.checkBox_As)
        
        # GUI operation-function 
        self.uid.buttonBox.connect(self.uid.buttonBox, SIGNAL("accepted()"), self.setScopeList)
        self.uid.buttonBox.connect(self.uid.applyButton, SIGNAL("clicked()"), self.setScopeList)

    def keyPressEvent(self, event):
        # operation of Key_Return or Key_Enter
        if(event.key() == 0x01000005 or event.key() == 0x01000004):
            pass
        else:
            # QDialog event
            super(OscilloscopeWindow, self).keyPressEvent(event) 
            
    def setObjectEnabled(self, from_self=True):
        # object enable setting
        if(self.MW.PSW.channel_enable['Snaf'] == True):
            self.uid.groupBox_Snaf.setEnabled(True)
        else:
            self.uid.groupBox_Snaf.setEnabled(False)
            self.uid.checkBox_SnaI.setCheckState(Qt.Unchecked)
            self.uid.checkBox_Snam.setCheckState(Qt.Unchecked)
            self.uid.checkBox_Snah.setCheckState(Qt.Unchecked)
            
        if(self.MW.PSW.channel_enable['Snap'] == True):
            self.uid.groupBox_Snap.setEnabled(True)
        else:
            self.uid.groupBox_Snap.setEnabled(False)
            self.uid.checkBox_SnapI.setCheckState(Qt.Unchecked)
            self.uid.checkBox_Snapm.setCheckState(Qt.Unchecked)
            
        if(self.MW.PSW.channel_enable['Skdr'] == True):
            self.uid.groupBox_Skdr.setEnabled(True)
        else:
            self.uid.groupBox_Skdr.setEnabled(False)
            self.uid.checkBox_SkdrI.setCheckState(Qt.Unchecked)
            self.uid.checkBox_Skdr.setCheckState(Qt.Unchecked)
            
        if(self.MW.PSW.channel_enable['Skca'] == True):
            self.uid.groupBox_Skca.setEnabled(True)
        else:
            self.uid.groupBox_Skca.setEnabled(False)
            self.uid.checkBox_SkcaI.setCheckState(Qt.Unchecked)
            
        if(self.MW.PSW.channel_enable['Scan'] == True):
            self.uid.groupBox_Scan.setEnabled(True)
        else:
            self.uid.groupBox_Scan.setEnabled(False)
            self.uid.checkBox_ScaI.setCheckState(Qt.Unchecked)
            self.uid.checkBox_Scam.setCheckState(Qt.Unchecked)
            self.uid.checkBox_Scah.setCheckState(Qt.Unchecked)
            
        if(self.MW.PSW.channel_enable['Sh'] == True):
            self.uid.groupBox_Sh.setEnabled(True)
        else:
            self.uid.groupBox_Sh.setEnabled(False)
            self.uid.checkBox_ShI.setCheckState(Qt.Unchecked)
            self.uid.checkBox_Shm.setCheckState(Qt.Unchecked)
            
        if(self.MW.PSW.channel_enable['Sesyn'] == True):
            self.uid.groupBox_Sesyn.setEnabled(True)
        else:
            self.uid.groupBox_Sesyn.setEnabled(False)
            self.uid.checkBox_SesynI.setCheckState(Qt.Unchecked)
            self.uid.checkBox_Sgesyn.setCheckState(Qt.Unchecked)
            
        if(self.MW.PSW.channel_enable['Sisyn'] == True):
            self.uid.groupBox_Sisyn.setEnabled(True)
        else:
            self.uid.groupBox_Sisyn.setEnabled(False)
            self.uid.checkBox_SisynI.setCheckState(Qt.Unchecked)
            self.uid.checkBox_Sgisyn.setCheckState(Qt.Unchecked)
            
        if(self.MW.PSW.channel_enable['Dcal'] == True):
            self.uid.groupBox_Dcal.setEnabled(True)
        else:
            self.uid.groupBox_Dcal.setEnabled(False)
            self.uid.checkBox_DcalI.setCheckState(Qt.Unchecked)
            self.uid.checkBox_Dcal.setCheckState(Qt.Unchecked)
            
        if(self.MW.PSW.channel_enable['Dnaf'] == True):
            self.uid.groupBox_Dnaf.setEnabled(True)
        else:
            self.uid.groupBox_Dnaf.setEnabled(False)
            self.uid.checkBox_DnaI.setCheckState(Qt.Unchecked)
            self.uid.checkBox_Dnam.setCheckState(Qt.Unchecked)
            self.uid.checkBox_Dnah.setCheckState(Qt.Unchecked)
            
        if(self.MW.PSW.channel_enable['Dnap'] == True):
            self.uid.groupBox_Dnap.setEnabled(True)
        else:
            self.uid.groupBox_Dnap.setEnabled(False)
            self.uid.checkBox_DnapI.setCheckState(Qt.Unchecked)
            self.uid.checkBox_Dnapm.setCheckState(Qt.Unchecked)
            
        if(self.MW.PSW.channel_enable['Dkdr'] == True):
            self.uid.groupBox_Dkdr.setEnabled(True)
        else:
            self.uid.groupBox_Dkdr.setEnabled(False)
            self.uid.checkBox_DkdrI.setCheckState(Qt.Unchecked)
            self.uid.checkBox_Dkdr.setCheckState(Qt.Unchecked)
            
        if(self.MW.PSW.channel_enable['Dkca'] == True):
            self.uid.groupBox_Dkca.setEnabled(True)
        else:
            self.uid.groupBox_Dkca.setEnabled(False)
            self.uid.checkBox_DkcaI.setCheckState(Qt.Unchecked)
            
        if(self.MW.PSW.channel_enable['Dcan'] == True):
            self.uid.groupBox_Dcan.setEnabled(True)
        else:
            self.uid.groupBox_Dcan.setEnabled(False)
            self.uid.checkBox_DcaI.setCheckState(Qt.Unchecked)
            self.uid.checkBox_Dcam.setCheckState(Qt.Unchecked)
            self.uid.checkBox_Dcah.setCheckState(Qt.Unchecked)
            
        if(self.MW.PSW.channel_enable['Dh'] == True):
            self.uid.groupBox_Dh.setEnabled(True)
        else:
            self.uid.groupBox_Dh.setEnabled(False)
            self.uid.checkBox_DhI.setCheckState(Qt.Unchecked)
            self.uid.checkBox_Dhm.setCheckState(Qt.Unchecked)
            
        if(self.MW.PSW.channel_enable['Desyn'] == True):
            self.uid.groupBox_Desyn.setEnabled(True)
        else:
            self.uid.groupBox_Desyn.setEnabled(False)
            self.uid.checkBox_DesynI.setCheckState(Qt.Unchecked)
            self.uid.checkBox_Dgesyn.setCheckState(Qt.Unchecked)
            
        if(self.MW.PSW.channel_enable['Disyn'] == True):
            self.uid.groupBox_Disyn.setEnabled(True)
        else:
            self.uid.groupBox_Disyn.setEnabled(False)
            self.uid.checkBox_DisynI.setCheckState(Qt.Unchecked)
            self.uid.checkBox_Dgisyn.setCheckState(Qt.Unchecked)
        
            self.setScopeList(from_self)
            
    def displayScope(self):
        if(self.MW.ModelType == 'Motoneuron' or self.MW.ModelType == 'Motor Unit'):        
            self.setObjectEnabled()        
            # checkbox check or uncheck           
            for i in range(self.mn_list_num):
                if(self.mn_checkbox_state[i] == False):
                    self.mn_checkbox_list[i].setCheckState(Qt.Unchecked)
                else: 
                    self.mn_checkbox_list[i].setCheckState(Qt.Checked)
        
        if(self.MW.ModelType == 'Muscle Fibers' or self.MW.ModelType == 'Motor Unit'):
            for i in range(self.mf_list_num):
                if(self.mf_checkbox_state[i] == False):
                    self.mf_checkbox_list[i].setCheckState(Qt.Unchecked)
                else:
                    self.mf_checkbox_list[i].setCheckState(Qt.Checked)
        
        # tab enable
        if(self.MW.ModelType == 'Motoneuron'):            
            self.uid.tabWidget.setTabEnabled(0, True)
            self.uid.tabWidget.setTabEnabled(1, False)
        elif(self.MW.ModelType == 'Muscle Fibers'):            
            self.uid.tabWidget.setTabEnabled(0, False)
            self.uid.tabWidget.setTabEnabled(1, True)       
        elif(self.MW.ModelType == 'Motor Unit'):            
            self.uid.tabWidget.setTabEnabled(0, True)
            self.uid.tabWidget.setTabEnabled(1, True)
            self.uid.tabWidget.setCurrentIndex(0)
        
        self.displayResult.setChecked(True)

    def setScopeList(self, from_self=True):
        # set the scope
        if(self.MW.ModelType == 'Motoneuron' or self.MW.ModelType == 'Motor Unit'):            
            for i in range(self.mn_list_num):                
                if(self.mn_checkbox_list[i].isChecked()):
                    self.mn_checkbox_state[i] = True
                    self.mn_ScopeList.add(self.mn_itemList[i])
                
                elif(not(self.mn_checkbox_list[i].isChecked()) and self.mn_checkbox_state[i] == True):
                    self.mn_checkbox_state[i] = False
                    self.mn_ScopeList.remove(self.mn_itemList[i])

        if(self.MW.ModelType == 'Muscle Fibers' or self.MW.ModelType == 'Motor Unit'):            
            for i in range(self.mf_list_num):                
                if(self.mf_checkbox_list[i].isChecked()):
                    self.mf_checkbox_state[i] = True
                    self.mf_ScopeList.add(self.mf_itemList[i])
                
                elif(not(self.mf_checkbox_list[i].isChecked()) and self.mf_checkbox_state[i] == True):
                    self.mf_checkbox_state[i] = False
                    self.mf_ScopeList.remove(self.mf_itemList[i])
       
        if(from_self==True):
            # Text message
            if(self.MW.ModelType == 'Motoneuron'):
                self.MW.setTextEdit("[Motoneuron] Output variables for plotting set.")
                item = "Item : "
                mn_num = len(self.mn_ScopeList)
                j = 1
                
                for i in self.mn_ScopeList:
                    if(j == mn_num):
                        item = item + str(i)
                    else:
                        item = item + str(i) + ", "
                    j += 1
                
                self.MW.uim.textEdit.append(item)
            
            elif(self.MW.ModelType == 'Muscle Fibers'):
                self.MW.setTextEdit("[Muscle Fibers] Output variables for plotting set.")
                item = "Item : "
                mf_num = len(self.mf_ScopeList)
                j = 1
                
                for i in self.mf_ScopeList:
                    if(j == mf_num):
                        item = item + str(i)
                    else:
                        item = item + str(i) + ", "
                    j += 1
                
                self.MW.uim.textEdit.append(item)
            
            elif(self.MW.ModelType == 'Motor Unit'):
                self.MW.setTextEdit("[Motor Unit] Output variables for plotting set.")
                item = "Item : [MN] "
                mn_num = len(self.mn_ScopeList)
                j = 1
                
                for i in self.mn_ScopeList:
                    if(j == mn_num):
                        item = item + str(i)
                    else:
                        item = item + str(i) + ", "
                    j += 1
                
                self.MW.uim.textEdit.append(item)
    
                item = "Item : [MF] "
                mf_num = len(self.mf_ScopeList)
                j = 1
                
                for i in self.mf_ScopeList:
                    if(j == mf_num):
                        item = item + str(i)
                    else:
                        item = item + str(i) + ", "
                    j += 1
                
                self.MW.uim.textEdit.append(item)
            
            # radio button set
            if(self.uid.radioButton_Individual.isChecked()):
                self.displayResult = self.uid.radioButton_Individual
                self.MW.uim.textEdit.append("Display : Separate") 
            elif(self.uid.radioButton_Combined.isChecked()):
                self.displayResult = self.uid.radioButton_Combined
                self.MW.uim.textEdit.append("Display : Overlapped")


if __name__  ==  '__main__':
    
    app=QApplication(sys.argv)

    # Main window instance
    MW = MainWindow()
    MW.show()
    app.exec_()
    
    del MW
    