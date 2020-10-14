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

import numpy as np
import matplotlib.pyplot as plt
import time
from math import log, tanh, exp, cosh
from scipy import integrate, signal
from PyQt4.QtCore import *

# Motoneuron class
class MotoNeuron:
    def __init__(self,uniqueNumber): 
        # initialization
        self.cellType='Motoneuron'
        self.uniqueNumber=uniqueNumber
        self.parameters=None
        self.ivalues=None
        self.cellState='Normal'
        self.figure = None
        self.Is=[]
        self.Se_syncon=[]
        self.Si_syncon=[]
        self.De_syncon=[]
        self.Di_syncon=[]
        self.SpikeTimes=[] # spike detection time array
        self.FiringRate=[] # spike rate array
        self.simulTime=0.

    # Motoneuron ODEs model
    def model(self, t, y):
        # motoneuron variables
        vs, sca, snam, snah, snapm, skdr, scam, scah, shm, vd,  dca, dcal, dnam, dnah, dnapm, dkdr, dcam, dcah, dhm = y
        # motoneuron constants
        rn,tm,VAsdDC,VAdsDC,VAsdAC,parea,svl,dvl,cv, sf,skca,salpha,sCAo,sEca, sgna,svna,sanamc,sanamv,sanama,sanamb,sbnamc,sbnamv,sbnama,sbnamb,snahth,snahslp,snahv,snaha,snahb,snahc, sgnap,svna,sanapmc,sanapmv,sanapma,sanapmb,sbnapmc,sbnapmv,sbnapma,sbnapmb, sgkdr,svk,skdrth,skdrslp,skdrv,skdra,skdrb,skdrc, sgkca,svk,skd, sgca,scamth,scamslp,scamtau,scahth,scahslp,scahtau, sgh,svh,shth,shslp,shtau, svesyn, svisyn, df,dkca,dalpha,dCAo,dEca, dgcal,dcalth,dcalslp,dcaltau, dgna,dvna,danamc,danamv,danama,danamb,dbnamc,dbnamv,dbnama,dbnamb,dnahth,dnahslp,dnahv,dnaha,dnahb,dnahc, dgnap,dvna,danapmc,danapmv,danapma,danapmb,dbnapmc,dbnapmv,dbnapma,dbnapmb, dgkdr,dvk,dkdrth,dkdrslp,dkdrv,dkdra,dkdrb,dkdrc,dgkca,dvk,dkd, dgca,dcamth,dcamslp,dcamtau,dcahth,dcahslp,dcahtau, dgh,dvh,dhth,dhslp,dhtau, dvesyn, dvisyn=self.parameters
        
        # Isoma
        if(self.IsSignalType=='Step'):
            i0, ip1, pon1, poff1, ip2, pon2, poff2, ip3, pon3, poff3, ip4, pon4, poff4, ip5, pon5, poff5, s = self.heav_param
            I_s = i0 + s*((self.heav(poff1-t)*self.heav(t-pon1)*ip1)
                         + (self.heav(poff2-t)*self.heav(t-pon2)*ip2)
                         + (self.heav(poff3-t)*self.heav(t-pon3)*ip3)
                         + (self.heav(poff4-t)*self.heav(t-pon4)*ip4)
                         + (self.heav(poff5-t)*self.heav(t-pon5)*ip5))
        elif(self.IsSignalType=='Ramp'):
            I_s=self.pv-((self.pv-self.iv)/self.p)*abs(t-self.p)
        elif(self.IsSignalType=='Import'):
            I_s= np.interp(t, self.times,  self.Is, 0, 0)
        I_s = I_s/3.1576 # unit conversion (nA -> mA/cm^2)

        # synaptic conductance signals
        if(self.se_ch == True):
            sgesyn = np.interp(t, self.sesyn_t, self.Se_syncon, 0, 0) # interpolation
        else:
            sgesyn = 0.
        if(self.si_ch == True):
            sgisyn = np.interp(t, self.sisyn_t, self.Si_syncon, 0, 0)
        else:
            sgisyn = 0.
        if(self.de_ch == True):
            dgesyn = np.interp(t, self.desyn_t, self.De_syncon, 0, 0)
        else:
            dgesyn = 0.
        if(self.di_ch == True):
            dgisyn = np.interp(t, self.disyn_t, self.Di_syncon, 0, 0)
        else:
            dgisyn = 0.
        
        ## [SOMA]
        # prevent a numerical error
        if(sca < 1.e-100):
            sca = 1.e-100
        # Ca2+ reversal potential
        if(self.const_sEca == True):
            svca = sEca
        else:
            svca = self.const*log(sCAo/sca)-70
        # N-Type Ca2+ channels current
        sica=sgca*scam**2*scah*(vs-svca)     
        # Fast Na+ channels current
        sina=sgna*snam**3*snah*(vs-svna)
        # Persistent Na+ channels current
        sinap=sgnap*snapm**3*(vs-svna)
        # HCN channels current
        sih=sgh*shm*(vs-svh)
        # Synaptic channels current
        siesyn=sgesyn*(vs-svesyn)
        siisyn=sgisyn*(vs-svisyn)
        sisyn=siesyn+siisyn
        # inward voltage-gated current 
        ins=-sina-sica-sinap-sih-sisyn
        # Delayed rectifier K+ channels current
        sikdr=sgkdr*skdr**4*(vs-svk)
        # Ca2+ dependent K+ channels current
        sikca=sgkca*(sca/(sca+skd))*(vs-svk)
        # outward voltage-gated current
        outs=-sikdr-sikca
        
        ## [Dendrite]
        # prevent a numerical error
        if(dca < 1.e-100):
            dca = 1.e-100        
        # Ca2+ reversal potential
        if(self.const_dEca == True):
            dvca = dEca
        else:
            dvca=self.const*log(dCAo/dca)-70
        # L-Type Ca2+ (CAv 1.3) channels current
        dical=dgcal*dcal*(vd-dvca) 
        # N-Type Ca2+ channels current
        dica=dgca*dcam**2*dcah*(vd-dvca)
        # Fast Na+ channels current
        dina=dgna*dnam**3*dnah*(vd-dvna)
        # Persistent Na+ channels current
        dinap=dgnap*dnapm**3*(vd-dvna)
        # HCN channels current 
        dih=dgh*dhm*(vd-dvh)
        # Synaptic channels current
        diesyn=dgesyn*(vd-dvesyn)
        diisyn=dgisyn*(vd-dvisyn)
        disyn=diesyn+diisyn
        # inward voltage-gated current
        ind=-dical-dina-dica-dinap-dih-disyn
        # Delayed rectifier K+ channels current
        dikdr=dgkdr*dkdr**4*(vd-dvk)
        # Ca2+ dependent K+ channels current
        dikca=dgkca*(dca/(dca+dkd))*(vd-dvk)
        # outward voltage-gated current 
        outd=-dikdr-dikca

        # Output from ODEs
        n = len(y)      
        dydt=range(n)
        ## [SOMA]
        # d(vs)/dt
        dydt[0]=(I_s+ins+outs-self.gms*(vs-svl)+self.gc*(vd-vs)/parea)/self.cms
        # d(sca)/dt
        dydt[1]=-sf*salpha*sica-sf*skca*sca
        # d(snam)/dt 
        alpha_snafm=sanamc*(vs-sanamv)/(exp(-(vs-sanamv)/sanama)+sanamb)
        beta_snafm=sbnamc*(vs-sbnamv)/(exp((vs-sbnamv)/sbnama)+sbnamb)
        dydt[2]=alpha_snafm*(1-snam)-beta_snafm*snam
        # d(snah)/dt
        snahinf=1.0/(1.0+exp((vs-snahth)/snahslp))
        snahtau=snahc/(exp((vs-snahv)/snaha)+exp(-(vs-snahv)/snahb))
        dydt[3]=(snahinf-snah)/snahtau
        # d(snapm))/dt
        alpha_snapm=sanapmc*(vs-sanapmv)/(exp(-(vs-sanapmv)/sanapma)+sanapmb)
        beta_snapm=sbnapmc*(vs-sbnapmv)/(exp((vs-sbnapmv)/sbnapma)+sbnapmb)
        dydt[4]=alpha_snapm*(1-snapm)-beta_snapm*snapm
        # d(skdr)/dt
        skdrinf=1.0/(1.0+exp(-(vs-skdrth)/skdrslp))
        skdrtau=skdrc/(exp((vs-skdrv)/skdra)+exp(-(vs-skdrv)/skdrb))
        dydt[5]=(skdrinf-skdr)/skdrtau
        # d(scam)/dt
        scaminf=1.0/(1.0+exp(-(vs-scamth)/scamslp))
        dydt[6]=(scaminf-scam)/scamtau
        # d(scah)/dt
        scahinf=1.0/(1.0+exp((vs-scahth)/scahslp))
        dydt[7]=(scahinf-scah)/scahtau
        # d(shm)/dt
        shminf=1/(1+exp((vs+shth)/shslp))
        dydt[8]=(shminf-shm)/shtau
        ## [DENDRITE]
        # d(vd)/dt 
        dydt[9]=(ind+outd-self.gmd*(vd-dvl)+self.gc*(vs-vd)/(1.0-parea))/self.cmd
        # sum of calcium current
        sdica = dical+dica
        # d(dca)/dt
        dydt[10]=-df*dalpha*sdica-df*dkca*dca
        # d(dcal)/dt
        dcalinf=1.0/(1.0+exp(-(vd-dcalth)/dcalslp))
        dydt[11]=(dcalinf-dcal)/dcaltau
        # d(dnam)/dt
        alpha_dnafm=danamc*(vd-danamv)/(exp(-(vd-danamv)/danama)+danamb)
        beta_dnafm=dbnamc*(vd-dbnamv)/(exp((vd-dbnamv)/dbnama)+dbnamb)
        dydt[12]=alpha_dnafm*(1-dnam)-beta_dnafm*dnam 
        # d(dnah)/dt
        dnahinf=1.0/(1.0+exp((vd-dnahth)/dnahslp))
        dnahtau=dnahc/(exp((vd-dnahv)/dnaha)+exp(-(vd-dnahv)/dnahb))
        dydt[13]=(dnahinf-dnah)/dnahtau
        # d(dnapm)/dt
        alpha_dnapm=danapmc*(vd-danapmv)/(exp(-(vd-danapmv)/danapma)+danapmb)
        beta_dnapm=dbnapmc*(vd-dbnapmv)/(exp((vd-dbnapmv)/dbnapma)+dbnapmb)
        dydt[14]=alpha_dnapm*(1-dnapm)-beta_dnapm*dnapm
        # d(dkdr)/dt
        dkdrinf=1.0/(1.0+exp(-(vd-dkdrth)/dkdrslp))
        dkdrtau=dkdrc/(exp((vd-dkdrv)/dkdra)+exp(-(vd-dkdrv)/dkdrb))
        dydt[15]=(dkdrinf-dkdr)/dkdrtau
        # d(dcam)/dt
        dcaminf=1.0/(1.0+exp(-(vd-dcamth)/dcamslp))
        dydt[16]=(dcaminf-dcam)/dcamtau
        # d(dcah)/dt
        dcahinf=1.0/(1.0+exp((vd-dcahth)/dcahslp))
        dydt[17]=(dcahinf-dcah)/dcahtau
        # d(dhm)/dt
        dhminf=1/(1+exp((vd+dhth)/dhslp))
        dydt[18]=(dhminf-dhm)/dhtau
       
        return dydt

    # Integration function
    def solModel(self, scope1, display_result):
        # motoneuron constants
        rn,tm,VAsdDC,VAdsDC,VAsdAC,parea,svl,dvl,cv, sf,skca,salpha,sCAo,sEca, sgna,svna,sanamc,sanamv,sanama,sanamb,sbnamc,sbnamv,sbnama,sbnamb,snahth,snahslp,snahv,snaha,snahb,snahc, sgnap,svna,sanapmc,sanapmv,sanapma,sanapmb,sbnapmc,sbnapmv,sbnapma,sbnapmb, sgkdr,svk,skdrth,skdrslp,skdrv,skdra,skdrb,skdrc, sgkca,svk,skd, sgca,scamth,scamslp,scamtau,scahth,scahslp,scahtau, sgh,svh,shth,shslp,shtau, svesyn, svisyn, df,dkca,dalpha,dCAo,dEca, dgcal,dcalth,dcalslp,dcaltau, dgna,dvna,danamc,danamv,danama,danamb,dbnamc,dbnamv,dbnama,dbnamb,dnahth,dnahslp,dnahv,dnaha,dnahb,dnahc, dgnap,dvna,danapmc,danapmv,danapma,danapmb,dbnapmc,dbnapmv,dbnapma,dbnapmb, dgkdr,dvk,dkdrth,dkdrslp,dkdrv,dkdra,dkdrb,dkdrc,dgkca,dvk,dkd, dgca,dcamth,dcamslp,dcamtau,dcahth,dcahslp,dcahtau, dgh,dvh,dhth,dhslp,dhtau, dvesyn, dvisyn=self.parameters        
        # Isoma
        if(self.IsSignalType=='Step'):
            i0, ip1, pon1, poff1, ip2, pon2, poff2, ip3, pon3, poff3, ip4, pon4, poff4, ip5, pon5, poff5, s = self.heav_param
        
        # initialize array
        self.SpikeTimes=[] # spike detection time array
        self.FiringRate=[] # spike rate array
        
        # set integrator
        r1= integrate.ode(self.model).set_integrator('vode')
        # initialize integrator
        r1.set_initial_value(self.ivalues, self.t_start)
        
        # Result arrays
        num_steps = int(np.floor((self.t_stop - self.t_start)/self.t_dt) + 1)
        T = np.zeros(num_steps)
        IS = np.zeros(num_steps)
        VS = np.zeros(num_steps)
        SCA = np.zeros(num_steps)
        SVCA = np.zeros(num_steps)
        SNAI = np.zeros(num_steps)
        SNAM = np.zeros(num_steps)
        SNAH = np.zeros(num_steps)
        SNAPI = np.zeros(num_steps)
        SNAPM = np.zeros(num_steps) 
        SKDRI = np.zeros(num_steps)
        SKDR = np.zeros(num_steps)
        SKCAI = np.zeros(num_steps)
        SCAI = np.zeros(num_steps)
        SCAM = np.zeros(num_steps) 
        SCAH = np.zeros(num_steps)
        SHI = np.zeros(num_steps)
        SHM = np.zeros(num_steps)
        VD = np.zeros(num_steps) 
        DCA = np.zeros(num_steps)
        DVCA = np.zeros(num_steps)
        DCALI = np.zeros(num_steps)
        DCAL = np.zeros(num_steps) 
        DNAI = np.zeros(num_steps)
        DNAM = np.zeros(num_steps)
        DNAH = np.zeros(num_steps)
        DNAPI = np.zeros(num_steps)
        DNAPM = np.zeros(num_steps)
        DKDRI = np.zeros(num_steps)
        DKDR = np.zeros(num_steps)
        DKCAI = np.zeros(num_steps)
        DCAI = np.zeros(num_steps)
        DCAM = np.zeros(num_steps)
        DCAH = np.zeros(num_steps)
        DHI = np.zeros(num_steps)
        DHM = np.zeros(num_steps)
        SESYNI = np.zeros(num_steps)
        SGESYN = np.zeros(num_steps)
        SISYNI = np.zeros(num_steps)
        SGISYN = np.zeros(num_steps)
        DESYNI = np.zeros(num_steps)
        DGESYN = np.zeros(num_steps)
        DISYNI = np.zeros(num_steps)
        DGISYN = np.zeros(num_steps)
            
        # set initial value 
        IS[0]= self.Is_0
        VS[0]=self.ivalues[0]
        SCA[0]=self.ivalues[1]
        if(self.const_sEca == True):
            SVCA[0] = sEca
        else:
            SVCA[0]=(self.const*log(sCAo/SCA[0]))-70
        SNAM[0]=self.ivalues[2]
        SNAH[0]=self.ivalues[3]
        SNAI[0] = sgna*SNAM[0]**3*SNAH[0]*(VS[0]-svna)
        SNAPM[0]=self.ivalues[4]
        SNAPI[0] = sgnap*SNAPM[0]**3*(VS[0]-svna)
        SKDR[0]=self.ivalues[5]
        SKDRI[0] = sgkdr*SKDR[0]**4*(VS[0]-svk)
        SKCAI[0] = sgkca*(SCA[0]/(SCA[0]+skd))*(VS[0]-svk)
        SCAM[0]=self.ivalues[6]
        SCAH[0]=self.ivalues[7]
        SCAI[0] = sgca*SCAM[0]**2*SCAH[0]*(VS[0]-SVCA[0])
        SHM[0]=self.ivalues[8]
        SHI[0] = sgh*SHM[0]*(VS[0]-svh)
        VD[0]=self.ivalues[9]
        DCA[0]=self.ivalues[10]
        if(self.const_dEca == True):
            DVCA[0] = dEca
        else:
            DVCA[0]=(self.const*log(dCAo/DCA[0]))-70
        DCAL[0]=self.ivalues[11]
        DCALI[0] = dgcal*DCAL[0]*(VD[0]-DVCA[0])
        DNAM[0]=self.ivalues[12]
        DNAH[0]=self.ivalues[13]
        DNAI[0] = dgna*DNAM[0]**3*DNAH[0]*(VD[0]-dvna)
        DNAPM[0]=self.ivalues[14]
        DNAPI[0] = dgnap*DNAPM[0]**3*(VD[0]-dvna)
        DKDR[0]=self.ivalues[15]
        DKDRI[0] = dgkdr*DKDR[0]**4*(VD[0]-dvk)
        DKCAI[0] = dgkca*(DCA[0]/(DCA[0]+dkd))*(VD[0]-dvk)
        DCAM[0]=self.ivalues[16]
        DCAH[0]=self.ivalues[17]
        DCAI[0] = dgca*DCAM[0]**2*DCAH[0]*(VD[0]-DVCA[0])
        DHM[0]=self.ivalues[18]
        DHI[0] = dgh*DHM[0]*(VD[0]-dvh)
        SGESYN[0] = self.Se_syncon[0]
        SESYNI[0] = SGESYN[0]*(VS[0]-svesyn)
        SGISYN[0] = self.Si_syncon[0]
        SISYNI[0] = SGISYN[0]*(VS[0]-svisyn)
        DGESYN[0] = self.De_syncon[0]
        DESYNI[0] = DGESYN[0]*(VD[0]-dvesyn)
        DGISYN[0] = self.Di_syncon[0]
        DISYNI[0] = DGISYN[0]*(VD[0]-dvisyn)
         
        ResultArrays={'Time': T,
                      'Is' : IS,
                      'V_soma': VS,
                      '[Ca]_soma': SCA,
                      'E_Ca_soma': SVCA,
                      'I_Naf_soma': SNAI,
                      'm_Naf_soma': SNAM,
                      'h_Naf_soma': SNAH,
                      'I_Nap_soma': SNAPI,
                      'm_Nap_soma': SNAPM,
                      'I_Kdr_soma': SKDRI,
                      'n_Kdr_soma' : SKDR,
                      'I_Kca_soma': SKCAI,
                      'I_Can_soma': SCAI,
                      'm_Can_soma': SCAM,
                      'h_Can_soma': SCAH,
                      'I_H_soma': SHI,
                      'm_H_soma' : SHM,    
                      'V_dend': VD,
                      '[Ca]_dend': DCA,
                      'E_Ca_dend': DVCA,
                      'I_Cal_dend': DCALI,
                      'l_Cal_dend': DCAL,
                      'I_Naf_dend': DNAI,
                      'm_Naf_dend': DNAM,
                      'h_Naf_dend': DNAH,
                      'I_Nap_dend': DNAPI,
                      'm_Nap_dend': DNAPM,
                      'I_Kdr_dend': DKDRI,
                      'n_Kdr_dend' : DKDR,
                      'I_Kca_dend': DKCAI,
                      'I_Can_dend': DCAI,
                      'm_Can_dend': DCAM,
                      'h_Can_dend': DCAH,
                      'I_H_dend': DHI,
                      'm_H_dend' : DHM,
                      'I_esyn_soma': SESYNI,
                      'G_esyn_soma': SGESYN,
                      'I_isyn_soma': SISYNI,
                      'G_isyn_soma': SGISYN,
                      'I_esyn_dend': DESYNI,
                      'G_esyn_dend': DGESYN,
                      'I_isyn_dend': DISYNI,
                      'G_isyn_dend': DGISYN }  
        
        k = 0
        # save integration start time
        start_time = time.time()
        
        # output graph creation
        if(len(scope1) != 0):
            self.createPlot(scope1, display_result)
            self.updatePlot(scope1, display_result, T, ResultArrays, k, num_steps)
        else:
            # figure instance creation for prevent GIL 
            self.figure=plt.figure(frameon=False)
            
        k = 1
        while r1.successful() and k < num_steps:
            # Integration with motoneuron
            r1.integrate(round(r1.t, 9) + self.t_dt)
            ## integration results
            T[k]=r1.t
            VS[k]=r1.y[0]
            SCA[k]=r1.y[1]
            SNAM[k]=r1.y[2]
            SNAH[k]=r1.y[3]
            SNAPM[k]=r1.y[4]
            SKDR[k]=r1.y[5]
            SCAM[k]=r1.y[6]
            SCAH[k]=r1.y[7]
            SHM[k]=r1.y[8]
            VD[k]=r1.y[9]
            DCA[k]=r1.y[10]
            DCAL[k]=r1.y[11]
            DNAM[k]=r1.y[12]
            DNAH[k]=r1.y[13]
            DNAPM[k]=r1.y[14]
            DKDR[k]=r1.y[15]
            DCAM[k]=r1.y[16]
            DCAH[k]=r1.y[17]
            DHM[k]=r1.y[18]
            
            ## re-calculate results about not state variable
            # Isoma            
            if(self.IsSignalType=='Step'):
                IS[k] = i0 + s*((self.heav(poff1-T[k])*self.heav(T[k]-pon1)*ip1)
                         + (self.heav(poff2-T[k])*self.heav(T[k]-pon2)*ip2)
                         + (self.heav(poff3-T[k])*self.heav(T[k]-pon3)*ip3)
                         + (self.heav(poff4-T[k])*self.heav(T[k]-pon4)*ip4)
                         + (self.heav(poff5-T[k])*self.heav(T[k]-pon5)*ip5))
            elif(self.IsSignalType=='Ramp'):
                IS[k] = self.pv-((self.pv-self.iv)/self.p)*abs(T[k]-self.p)
            elif(self.IsSignalType=='Import'):
                IS[k] = np.interp(T[k], self.times,  self.Is, 0, 0)
                
            # Spike detect & calculate Spike rate
            S_Detect = self.detect_Spike(T[k], VS[k-2], VS[k-1], VS[k])
            if(S_Detect == True):
                self.cal_FiringRate(self.SpikeTimes)
            
            ## [SOMA]
            if (SCA[k] < 1.e-100):
                SCA[k] = 1.e-100
            # Calcium dynamics
            if(self.const_sEca == True):
                SVCA[k] = sEca
            else:
                SVCA[k]=(self.const*log(sCAo/SCA[k]))-70
            # Can
            SCAI[k] = sgca*SCAM[k]**2*SCAH[k]*(VS[k]-SVCA[k])
            # Naf
            SNAI[k] = sgna*SNAM[k]**3*SNAH[k]*(VS[k]-svna)
            # Nap
            SNAPI[k] = sgnap*SNAPM[k]**3*(VS[k]-svna)
            # Kdr
            SKDRI[k] = sgkdr*SKDR[k]**4*(VS[k]-svk)
            # Kca
            SKCAI[k] = sgkca*(SCA[k]/(SCA[k]+skd))*(VS[k]-svk)
            # H
            SHI[k] = sgh*SHM[k]*(VS[k]-svh)
            # eSyn
            if(self.se_ch == True):
                SGESYN[k] = np.interp(T[k], self.sesyn_t, self.Se_syncon, 0, 0) # interpolation
            else:
                SGESYN[k] = 0.
            SESYNI[k] = SGESYN[k]*(VS[k]-svesyn)
            # iSyn
            if(self.si_ch == True):
                SGISYN[k] = np.interp(T[k], self.sisyn_t, self.Si_syncon, 0, 0)
            else:
                SGISYN[k] = 0.
            SISYNI[k] = SGISYN[k]*(VS[k]-svisyn)
            
            ## DENDRITE
            if (DCA[k] < 1.e-100):
                DCA[k] = 1.e-100
            # Calcium dynamics
            if(self.const_dEca == True):
                DVCA[k] = dEca
            else:
                DVCA[k]=(self.const*log(dCAo/DCA[k]))-70
            # Cal
            DCALI[k] = dgcal*DCAL[k]*(VD[k]-DVCA[k])
            # Naf
            DNAI[k] = dgna*DNAM[k]**3*DNAH[k]*(VD[k]-dvna)
            # Nap
            DNAPI[k] = dgnap*DNAPM[k]**3*(VD[k]-dvna)
            # Kdr
            DKDRI[k] = dgkdr*DKDR[k]**4*(VD[k]-dvk)
            # Kca
            DKCAI[k] = dgkca*(DCA[k]/(DCA[k]+dkd))*(VD[k]-dvk)
            # Can
            DCAI[k] = dgca*DCAM[k]**2*DCAH[k]*(VD[k]-DVCA[k])
            # H
            DHI[k] = dgh*DHM[k]*(VD[k]-dvh)
            # eSyn
            if(self.de_ch == True):
                DGESYN[k] = np.interp(T[k], self.desyn_t, self.De_syncon, 0, 0) # interpolation
            else:
                DGESYN[k] = 0.
            DESYNI[k] = DGESYN[k]*(VD[k]-dvesyn)
            # iSyn            
            if(self.di_ch == True):
                DGISYN[k] = np.interp(T[k], self.disyn_t, self.Di_syncon, 0, 0)
            else:
                DGISYN[k] = 0.
            DISYNI[k] = DGISYN[k]*(VD[k]-dvisyn)

            # output plotting
            if(len(scope1) != 0):
                self.updatePlot(scope1, display_result, T, ResultArrays, k, num_steps)
                
            else:
                # GUI event update
                self.figure.canvas.flush_events()
            
            k += 1
            
            # If user pushes Stop button, cancel the integration
            if (self.cellState=='Stop'):
                # save the results so far
                T = T[:k]
                IS = IS[:k]
                VS = VS[:k]
                SCA = SCA[:k]
                SVCA = SVCA[:k]
                SNAI = SNAI[:k]
                SNAM = SNAM[:k]
                SNAH = SNAH[:k]
                SNAPI = SNAPI[:k]
                SNAPM = SNAPM[:k]
                SKDRI = SKDRI[:k]
                SKDR = SKDR[:k]
                SKCAI = SKCAI[:k]
                SCAI = SCAI[:k]
                SCAM = SCAM[:k]
                SCAH = SCAH[:k]
                SHI = SHI[:k]
                SHM = SHM[:k]
                VD = VD[:k]
                DCA = DCA[:k]
                DVCA = DVCA[:k]
                DCALI = DCALI[:k]
                DCAL = DCAL[:k]
                DNAI = DNAI[:k]
                DNAM = DNAM[:k]
                DNAH = DNAH[:k]
                DNAPI = DNAPI[:k]
                DNAPM = DNAPM[:k]
                DKDRI = DKDRI[:k]
                DKDR = DKDR[:k]
                DKCAI = DKCAI[:k]
                DCAI = DCAI[:k]
                DCAM = DCAM[:k]
                DCAH = DCAH[:k]
                DHI = DHI[:k]
                DHM = DHM[:k]
                SESYNI = SESYNI[:k]
                SGESYN = SGESYN[:k]
                SISYNI = SISYNI[:k]
                SGISYN = SGISYN[:k]
                DESYNI = DESYNI[:k]
                DGESYN = DGESYN[:k]
                DISYNI = DISYNI[:k]
                DGISYN = DGISYN[:k]

                ResultArrays={'Time': T,
                      'Is' : IS,
                      'V_soma': VS,
                      '[Ca]_soma': SCA,
                      'E_Ca_soma': SVCA,
                      'I_Naf_soma': SNAI,
                      'm_Naf_soma': SNAM,
                      'h_Naf_soma': SNAH,
                      'I_Nap_soma': SNAPI,
                      'm_Nap_soma': SNAPM,
                      'I_Kdr_soma': SKDRI,
                      'n_Kdr_soma' : SKDR,
                      'I_Kca_soma': SKCAI,
                      'I_Can_soma': SCAI,
                      'm_Can_soma': SCAM,
                      'h_Can_soma': SCAH,
                      'I_H_soma': SHI,
                      'm_H_soma' : SHM,    
                      'V_dend': VD,
                      '[Ca]_dend': DCA,
                      'E_Ca_dend': DVCA,
                      'I_Cal_dend': DCALI,
                      'l_Cal_dend': DCAL,
                      'I_Naf_dend': DNAI,
                      'm_Naf_dend': DNAM,
                      'h_Naf_dend': DNAH,
                      'I_Nap_dend': DNAPI,
                      'm_Nap_dend': DNAPM,
                      'I_Kdr_dend': DKDRI,
                      'n_Kdr_dend' : DKDR,
                      'I_Kca_dend': DKCAI,
                      'I_Can_dend': DCAI,
                      'm_Can_dend': DCAM,
                      'h_Can_dend': DCAH,
                      'I_H_dend': DHI,
                      'm_H_dend' : DHM,
                      'I_esyn_soma': SESYNI,
                      'G_esyn_soma': SGESYN,
                      'I_isyn_soma': SISYNI,
                      'G_isyn_soma': SGISYN,
                      'I_esyn_dend': DESYNI,
                      'G_esyn_dend': DGESYN,
                      'I_isyn_dend': DISYNI,
                      'G_isyn_dend': DGISYN }
                      
                if(len(scope1) == 0):
                    self.figure = None
                # calculate Elapsed time  
                #self.simulTime=time.time() - start_time    
                break
        
        # calculate Elapsed time   
        self.simulTime=time.time() - start_time
        
        if(len(scope1) == 0):
            self.figure = None
        else:
            # last graph update
            if(display_result == 'Individual'):
                for i in range(len(scope1)):
                    self.lines[i].set_animated(False)
                    self.background[i] = None
                    self.ax[i].relim()
                    self.ax[i].autoscale_view()
            elif(display_result == 'Combined'):
                for i in range(len(scope1)):
                    self.lines[i].set_animated(False)
                self.background = None
                self.ax.relim()
                self.ax.autoscale_view()
            self.figure.canvas.draw()
            # window close button Enable
            self.figure.canvas.parent().setWindowFlags(
            Qt.Window)
            self.figure.show()

        return ResultArrays
    
    def createPlot(self, scope1, display_result):
        self.scopeLength=len(scope1)
        self.sampling_rate=self.t_pt/self.t_dt
        self.xdata=range(self.scopeLength)
        self.ydata=range(self.scopeLength)
        self.lines=range(self.scopeLength)
        self.background=range(self.scopeLength)
        plt.rc('figure', figsize=(7, 7))
        
        # close previous figure
        if(self.figure != None):
            plt.close(self.figure)
        # new figure
        self.figure = plt.figure()
        self.figure.canvas.set_window_title('MOTONEURON')
            
        # Individual
        if(display_result == 'Individual'):
            self.ax=range(self.scopeLength)
            
            for i in range(self.scopeLength):
                var = scope1[i]
                # add axes
                self.ax[i] = self.figure.add_subplot(self.scopeLength, 1, i+1)
                self.ax[i].set_xlim(self.t_start, self.t_stop)
                self.ax[i].set_autoscaley_on(True)
                self.ax[i].grid()
                # Y label
                ylabel = var
                if(var=='V_soma' or var=='V_dend' or var=='E_Ca_soma' or var=='E_Ca_dend'):
                    ylabel += ' (mV)'
                elif(var=='Firing_rate'):
                    ylabel += ' (Hz)'
                elif(var=='Is'):
                    ylabel += ' (nA)'                
                elif(var=='G_esyn_soma' or var=='G_isyn_soma' or var=='G_esyn_dend' or var=='G_isyn_dend'):
                    ylabel += ' (mS/cm^2)'
                elif(var=='[Ca]_soma' or var=='[Ca]_dend'):
                    ylabel += ' (mM)'
                elif(var=='I_Naf_soma' or var=='I_Nap_soma' or var=='I_Kdr_soma' or var=='I_Kca_soma' or var=='I_Can_soma' or var=='I_H_soma' or var=='I_esyn_soma' or var=='I_isyn_soma' or var=='I_Cal_dend' or var=='I_Naf_dend' or var=='I_Nap_dend' or var=='I_Kdr_dend' or var=='I_Kca_dend' or var=='I_Can_dend' or var=='I_H_dend' or var=='I_esyn_dend' or var=='I_isyn_dend'):
                    ylabel += ' (mA/cm^2)'
                self.ax[i].set_ylabel(ylabel)
               
               # new Line with dot
                if(var=='Firing_rate'):
                    self.lines[i], = self.ax[i].plot([],[], '.', color='b', animated=True)
                # new Line with line
                else:
                    self.lines[i], = self.ax[i].plot([],[], color='b', animated=True)
                self.xdata[i]=[]
                self.ydata[i]=[]
            
            # X label
            if(i == (self.scopeLength-1)):
                self.ax[i].set_xlabel('Time (ms)')
            self.figure.canvas.draw()
            
            for i in range(self.scopeLength):
                # cache the background
                self.background[i] = self.figure.canvas.copy_from_bbox(self.ax[i].bbox)

        # Conbined
        elif(display_result == 'Combined'):
            # add axes
            self.ax = self.figure.add_subplot(1, 1, 1)
            self.ax.set_xlim(self.t_start, self.t_stop)
            self.ax.set_autoscaley_on(True)
            self.ax.grid()
            # X label
            self.ax.set_xlabel('Time (ms)')
            
            for i in range(self.scopeLength):
                var = scope1[i]
               # new Line with dot
                if(var=='Firing_rate'):
                    self.lines[i], = self.ax.plot([],[], '.', label=var, animated=True)
                # new Line with line
                else:
                    self.lines[i],=self.ax.plot([],[], label=var, animated=True)
                self.xdata[i]=[]
                self.ydata[i]=[]
                
                ## Y label
                # one variable
                if(self.scopeLength == 1):
                    if(var=='V_soma' or var=='V_dend' or var=='E_Ca_soma' or var=='E_Ca_dend'):
                        ylabel = ' (mV)'
                    elif(var=='Firing_rate'):
                        ylabel = ' (Hz)'
                    elif(var=='Is'):
                        ylabel = ' (nA)'                
                    elif(var=='G_esyn_soma' or var=='G_isyn_soma' or var=='G_esyn_dend' or var=='G_isyn_dend'):
                        ylabel = ' (mS/cm^2)'
                    elif(var=='[Ca]_soma' or var=='[Ca]_dend'):
                        ylabel = ' (mM)'
                    elif(var=='I_Naf_soma' or var=='I_Nap_soma' or var=='I_Kdr_soma' or var=='I_Kca_soma' or var=='I_Can_soma' or var=='I_H_soma' or var=='I_esyn_soma' or var=='I_isyn_soma' or var=='I_Cal_dend' or var=='I_Naf_dend' or var=='I_Nap_dend' or var=='I_Kdr_dend' or var=='I_Kca_dend' or var=='I_Can_dend' or var=='I_H_dend' or var=='I_esyn_dend' or var=='I_isyn_dend'):
                        ylabel = ' (mA/cm^2)'
                    self.ax.set_ylabel(ylabel)
            # two variable        
            if(self.scopeLength == 2):
                if((scope1[0]=='V_soma' and scope1[1] == 'V_dend') or (scope1[1]=='V_soma' and scope1[0] == 'V_dend')):                     
                    ylabel = 'V_soma & V_dend (mV)'
                    self.ax.set_ylabel(ylabel)
            
            self.ax.legend(loc='best')
            self.figure.canvas.draw()
            # cache the background
            self.background = self.figure.canvas.copy_from_bbox(self.ax.bbox)

        # window close button disable
        self.figure.canvas.parent().setWindowFlags(
        Qt.WindowTitleHint | 
        Qt.Dialog | 
        Qt.WindowMinimizeButtonHint | 
        Qt.CustomizeWindowHint)
        self.figure.show()

    def updatePlot(self, scope1, display_result, T, ResultArrays, k, num_steps):
        # Individual
        if(display_result == 'Individual'):
            for i in range(self.scopeLength):
                var = scope1[i]
                # x, y data update
                if(var=='Firing_rate'):
                    self.xdata[i] = self.SpikeTimes[1:] # except first value
                    self.ydata[i] = self.FiringRate
                else:
                    self.xdata[i].append(T[k])
                    self.ydata[i].append(ResultArrays[var][k])
                
                # plotting every t_pt step and last k.
                if(k%self.sampling_rate==0 or k==num_steps-1): 
                    self.lines[i].set_xdata(self.xdata[i]) 
                    self.lines[i].set_ydata(self.ydata[i])
                    ymin, ymax = self.ax[i].get_ylim()
                    
                    # scale update
                    if(var=='Firing_rate'):
                        if(len(self.ydata[i]) > 1):
                            if((self.ydata[i][-1] > ymax) or (self.ydata[i][-1] < ymin)):
                                self.ax[i].relim()
                                self.ax[i].autoscale_view()
                                self.figure.canvas.draw()
                    else:
                        if((ResultArrays[var][k] > ymax) or (ResultArrays[var][k] < ymin)):
                            self.ax[i].relim()
                            self.ax[i].autoscale_view()
                            self.figure.canvas.draw()
                    
                    # restore background
                    self.figure.canvas.restore_region(self.background[i])
                    self.ax[i].draw_artist(self.lines[i])
                    # fill in the axes rectangle (blit)
                    self.figure.canvas.blit(self.ax[i].bbox)
                    # figure update
                    self.ax[i].relim()
                    self.ax[i].autoscale_view()
                    self.figure.canvas.flush_events()
                    
        # Conbined
        elif(display_result == 'Combined'):
            for i in range(self.scopeLength):
                var = scope1[i]
                # x, y data update
                if(var=='Firing_rate'):
                    self.xdata[i] = self.SpikeTimes[1:]
                    self.ydata[i] = self.FiringRate
                else:
                    self.xdata[i].append(T[k]) 
                    self.ydata[i].append(ResultArrays[var][k])       
                
                # last plotting for resting part.
                if(k%self.sampling_rate==0 or k==num_steps-1):
                    self.lines[i].set_xdata(self.xdata[i])
                    self.lines[i].set_ydata(self.ydata[i])
                    ymin, ymax = self.ax.get_ylim()
                    
                    # scale update
                    if(var=='Firing_rate'):
                        if(len(self.ydata[i]) > 1):
                            if((self.ydata[i][-1] > ymax) or (self.ydata[i][-1] < ymin)): 
                                self.ax.relim()
                                self.ax.autoscale_view()
                                self.figure.canvas.draw()
                    else:
                        if((ResultArrays[var][k] > ymax) or (ResultArrays[var][k] < ymin)):
                            self.ax.relim()
                            self.ax.autoscale_view()
                            self.figure.canvas.draw()
                            
                    if(i == 0):
                        # restore background
                        self.figure.canvas.restore_region(self.background)
                    self.ax.draw_artist(self.ax.lines[i])

                    if(i == self.scopeLength-1):
                        # fill in the axes rectangle (blit)
                        self.figure.canvas.blit(self.ax.bbox)
                        # figure update
                        self.ax.relim()
                        self.ax.autoscale_view()
                        self.figure.canvas.flush_events()
                        
    def detect_Spike(self, t, Vs2, Vs1, Vs):
        # If rate of change of Vs bigger than Vth, it regards the spike occured.
        Vth=200
        dt=self.t_dt
        V1=(Vs1-Vs2)/dt
        V=(Vs-Vs1)/dt
        
        if((V1<Vth) and (V>Vth)):
            # save Spike Time
            self.SpikeTimes.append(t)
            S_Detect=True
        else:
            S_Detect=False

        return S_Detect
        
    def cal_FiringRate(self, SpikeTimes):
        i = len(SpikeTimes)-1
        if(i >= 1): 
            s = 1000
            # calculate Spike rate from Spike Time
            T1 = SpikeTimes[i-1]
            T = SpikeTimes[i]
            FiringRate = 1/(T-T1)*s
            self.FiringRate.append(FiringRate)

    def heav(self, x):
        # heavian function
        return np.where(x < 0, 0, 1)
        
    def setModelParam(self, parameters, const_sEca, const_dEca):
        # create 1D array
        temp_arr = []
        for item in parameters:
            temp_arr = temp_arr + item
        self.parameters = temp_arr
        parameters = temp_arr

        # Constants
        self.const_sEca = const_sEca
        self.const_dEca = const_dEca
        R=8.31441
        Temp=309.15
        Zca=2
        Fe=96485.309
        self.const = 1000*R*Temp/Zca/Fe
        w=1570.796
        rn = parameters[0]
        tm = parameters[1]
        VAsdDC = parameters[2]
        VAdsDC = parameters[3]
        VAsdAC = parameters[4]
        parea = parameters[5]
        rn *= 0.31576
        
        # Conductance Inverse equations with DDVA properties
        gms = (1.-VAdsDC)/(rn*(1.-VAsdDC*VAdsDC))
        gmd = (parea*VAdsDC*(1.-VAsdDC))/((1.-parea)*rn*VAsdDC*(1.-VAsdDC*VAdsDC))
        gc = (parea*VAdsDC)/(rn*(1.-VAsdDC*VAdsDC))
        cmd = (1./(w*(1.-parea)))*np.sqrt(((gc**2)/(VAsdAC**2))-((gc+gmd*(1.-parea))**2))
        cms = (tm*(parea*(1.-parea)*tm*gms*gmd+parea*gms*(tm*gc-cmd)+(parea**2)*gms*cmd+(1.-parea)*(tm*gc*gmd-gc*cmd)))/(parea*((1.-parea)*(tm*gmd-cmd)+(tm*gc)))
        self.gms = gms*(1.e-1)
        self.gmd = gmd*(1.e-1)
        self.gc = gc*(1.e-1)
        self.cmd = cmd*(1.e+2)
        self.cms = cms*(1.e+2)
      
    def setInputSignal(self, IsSignalType, IsiValue, IspValue, Is_0, heav_param, IsPeriod=0, times=0, Is=0):
        # set the Isoma parameters    
        self.IsSignalType=IsSignalType
        self.iv=IsiValue
        self.pv=IspValue
        self.p=IsPeriod
        self.Is_0 = Is_0
        self.heav_param = heav_param
        self.times = times
        self.Is=Is
        
    def setSynConSignal(self, se_ch, si_ch, de_ch, di_ch, Se_times, Si_times, De_times, Di_times, Se_syncon, Si_syncon, De_syncon, Di_syncon):
        # set the Isyn parameters
        self.Se_syncon=Se_syncon
        self.Si_syncon=Si_syncon
        self.De_syncon=De_syncon
        self.Di_syncon=Di_syncon
        if(se_ch == False):
            Se_syncon *= 0.
        if(si_ch == False):
            Si_syncon *= 0.
        if(de_ch == False):
            De_syncon *= 0.
        if(di_ch == False):
            Di_syncon *= 0.
        self.se_ch = se_ch
        self.si_ch = si_ch
        self.de_ch = de_ch
        self.di_ch = di_ch
        self.sesyn_t = Se_times
        self.sisyn_t = Si_times
        self.desyn_t = De_times
        self.disyn_t = Di_times
        
    def setInitialValues(self, ivalues):
        # create 1D array        
        temp_arr = []
        for item in ivalues:
            temp_arr.extend(item)
        self.ivalues = temp_arr
    
    def setIntegrationEnv(self, t_start, t_stop, t_dt, t_pt):
        self.t_start=t_start
        self.t_stop=t_stop
        self.t_dt=t_dt
        self.t_pt=t_pt
     
# Muscle fibers class 
class MuscleFibers:
    def __init__(self, uniqueNumber): 
        # initialization
        self.cellType='Muscle Fibers'
        self.uniqueNumber=uniqueNumber
        self.parameters=None
        self.ivalues=None
        self.sd = 0 # spike delay (ms)
        self.spike=[]
        self.spike_idx=[]
        self.xm = []
        self.vm = []
        self.am = []
        self.cellState='Normal'
        self.figure = None
        self.SpikeTimes=[] # spike detection time array
        self.simulTime=0.

    # Muscle fibers ODEs model    
    def model(self, t, y):
        # Muscle fibers variables
        CS,CaSR,CaSRCS,B,CaSP,T,CaSPB,CaSPT,A,XCE = y       
        # Muscle fibers constants
        K1, K2, K3, K4, K5i, K6i, K, Pmax, Umax, t1, t2,  u1, u2, u3, u4, C1, C2, C3, C4, C5, KSE,  P0, g1, g2, a0, b0, c0, d0  = self.parameters  
        ms=0.001
        b0, d0= b0*ms, d0*ms 
        
        # calculate Xm, Vm, Am
        Xm, Vm, Am=self.get_Xm_Vm_Am(t)
        
        K6=(K6i/(1+5*A))
                
        if(Xm<=-8):
            uXm=(u1*Xm)+u2
        if(Xm>-8):
            uXm=(u3*Xm)+u4
        
        K5=uXm*K5i
           
        # Output from ODEs
        n = len(y)      
        dydt=range(n)
        
        ## Module 1    
        # d(CS)/dt
        dydt[0]=-K1*CS*CaSR+K2*CaSRCS
        R=self.get_R(t, CaSR)
        U=Umax*(((CaSP**2)*(K**2))/(1+CaSP*K+(CaSP**2)*(K**2)))**2
        dydt[1]=-K1*CS*CaSR+K2*CaSRCS-R+U
        dydt[2]=K1*CS*CaSR-K2*CaSRCS
        dydt[3]=-K3*CaSP*B+K4*CaSPB
        dydt[4]=-K5*CaSP*T+K6*CaSPT-K3*CaSP*B+K4*CaSPB+R-U
        dydt[5]=-K5*CaSP*T+K6*CaSPT
        dydt[6]=K3*CaSP*B-K4*CaSPB
        dydt[7]=K5*CaSP*T-K6*CaSPT
        ## Module 2
        # d(A)/dt
        Ainf=0.5*(1+tanh(((CaSPT/(CaSPT+T))-C1)/C2))
        TA=C3/(cosh(((CaSPT/(CaSPT+T))-C4)/(2*C5)))
        dydt[8]=(Ainf-A)/TA
        # d(F)/dt
        gXm=exp(-((Xm-g1)/g2)**2)
        As=self.get_As(Xm,A,Vm,Am,t)
        F=self.get_F(t, XCE, Xm)
        Fc=P0*gXm*As
        if(F<=Fc):
            dydt[9]=(-b0*(Fc-F))/(F+a0*Fc/P0)
        else:
            gainlength=(-d0*(Fc-F))/(2*Fc-F+c0*Fc/P0)
            if(gainlength<=0):
                dydt[9]=(1.e-3)*(1.e5)
            else:
                dydt[9]=gainlength

        return dydt

    def setIntegrationEnv(self, t_start, t_stop, t_dt, t_pt):
        self.t_start=t_start
        self.t_stop=t_stop
        self.t_dt=t_dt
        self.t_pt=t_pt    

    def get_Xm_Vm_Am(self, t):        
        # interpolate Xm, Vm, Am
        if(self.XmSignalType=='Isometric'):
            xm=self.xm[0]
            vm=0.
            am=0.
        elif(self.XmSignalType=='Isokinetic'):
            xm=np.interp(t, self.times,  self.xm, -8, self.xm[-1])
            vm=np.interp(t, self.times,  self.vm, -1, -8)
            am=np.interp(t, self.times,  self.am, -1, -8)
        elif(self.XmSignalType=='Dynamic'):
            xm=np.interp(t, self.times,  self.xm)
            vm=np.interp(t, self.times,  self.vm)
            am=np.interp(t, self.times,  self.am)
        elif(self.XmSignalType=='Exp'):
            xm=np.interp(t, self.times,  self.xm, -8, self.xm[-1])
            vm=np.interp(t, self.times,  self.vm, -1, -8)
            am=np.interp(t, self.times,  self.am, -1, -8)
     
        return xm, vm, am        
         
    def get_F(self, t, XCE, Xm):
        # constants
        xm_init=self.xm[0]   
        xce_init=self.ivalues[9]
        d_xce=XCE-xce_init
        d_xm=Xm-xm_init
        d_se=d_xm-d_xce
        P0=self.parameters[21]
        KSE=self.parameters[20]

        # fix F initial value for preventing initial integration error.
        if(t<0.002): 
            F=1.e-5
        elif(d_se<=0):   
            F=0.
        else:
            F=P0*KSE*(d_se)

        return F
        
    def get_R(self, t ,CaSR):
        # constants
        Pmax = self.R_Pmax
        t1 = self.R_t1
        t2 = self.R_t2
        R=0.
        
        # calculate R with spikes time array
        for i in self.SpikeTimes:
            if(t>=i):
                R+=CaSR*Pmax*(1-exp(-(t-i)/t1))*exp(-(t-i)/t2)
       
        return R

    def get_As(self, Xm,A,Vm,Am,t):
        # constants   
        L0=-8
        beta=0.47
        ai=2.
        a1, a2, a3= [4.77, 400., 160.]
        r=0.001
        u1, u2, u3, u4 = self.parameters[11:15]
        
        if(Xm<=-8):
            uxm=(u1*Xm)+u2
        if(Xm>-8):
            uxm=(u3*Xm)+u4
        
        # calculate A_tilde
        if(Am == 0.):
            As=A**ai
        else:
            alpha=a1*(1+tanh((t-a2)/a3))+ai
            At=A**(alpha)
            
            if((Xm<L0) and (Vm>0)):
                Ad=(1+beta*uxm)*(1+r*Vm)
                As=At/Ad
            else:
                As=At
                
        if(As<1.e-3):
            As=1.e-3
           
        return As
    
    def detect_Spike(self, t, Vs2, Vs1, Vs):
        # If rate of change of Vs bigger than Vth, it regards the spike occured.
        Vth=200
        dt=self.t_dt
        V1=(Vs1-Vs2)/dt
        V=(Vs-Vs1)/dt
        
        if((V1<Vth) and (V>Vth)):
            # save Spike Time
            self.SpikeTimes.append(t+self.sd) # apply the spike delay
            S_Detect=True
        else:
            S_Detect=False

        return S_Detect
        
    def setModelParam(self, parameters):
        # create 1D array
        temp_arr = []
        for item in parameters:
            temp_arr = temp_arr + item
        self.parameters = temp_arr
        
        # set constants related to R
        self.R_Pmax=self.parameters[7]
        self.R_t1, self.R_t2=self.parameters[9], self.parameters[10]
        
    def setInitialValues(self, ivalues):
        # create 1D array
        temp_arr = []
        for item in ivalues:
            temp_arr = temp_arr + item
        self.ivalues = temp_arr
    
    def setSpikeSignal(self, spike, spike_idx, SpikeTimes):
        self.spike = spike
        self.SpikeTimes = SpikeTimes
        self.spike_idx = spike_idx

    def setXmSignal(self, XmSignalType, times, Xm):
        self.XmSignalType = XmSignalType
        self.times = times
        self.xm = Xm
        lengh = len(Xm)
        dt = self.t_dt 
        ms=0.001
        
        self.vm = np.zeros(lengh)
        self.am = np.zeros(lengh)
        self.vm[0] = 0.
        self.am[0] = 0.

        # calculate Vm, Am with Xm (Euler Method, t-1, t)
        if(self.XmSignalType=='Isometric'):
            return
        
        else:
            if(self.XmSignalType=='Exp' or self.XmSignalType=='Dynamic'):
                dt = times[1] - times[0]
                    
            for t in xrange(1, lengh):
                xm = Xm[t]
                if(t == 1):
                    xm_1 = Xm[0]
                    self.vm[t] = (xm-xm_1)/(dt*ms)                
                    self.am[t] = self.vm[t]/dt/ms
                    self.am[t] = 0
                else:
                    xm_1 = Xm[t-1]
                    xm_2 = Xm[t-2]
                    self.vm[t] = (xm-xm_1)/(dt*ms)
                    vm_1 = (xm_1-xm_2)/(dt*ms)
                    self.am[t] = (self.vm[t]-vm_1)/(dt*ms)
                
                # Filtering Am by 0.
                if(abs(self.am[t]) < 1.e-3):
                    self.am[t] = 0.
                    
    def createPlot(self, scope2, display_result):    
        self.scopeLength=len(scope2)
        self.sampling_rate=self.t_pt/self.t_dt
        self.lines=range(self.scopeLength)
        self.xdata=range(self.scopeLength)
        self.ydata=range(self.scopeLength)
        self.background=range(self.scopeLength)
        plt.rc('figure', figsize=(7, 7))
        
        # close previous figure
        if(self.figure != None):
            plt.close(self.figure)
        # new figure
        self.figure = plt.figure()
        self.figure.canvas.set_window_title('MUSCLE FIBERS')
            
        # Individual
        if(display_result == 'Individual'):
            self.ax=range(self.scopeLength)
            for i in range(self.scopeLength):
                var = scope2[i]
                # add axes
                self.ax[i] = self.figure.add_subplot(self.scopeLength, 1, i+1)
                self.ax[i].set_xlim(self.t_start, self.t_stop)
                self.ax[i].set_autoscaley_on(True)
                self.ax[i].grid()
                # Y label
                ylabel = var
                if(var=='F'):
                    ylabel += ' (N)'
                elif(var=='Xm'):
                    ylabel += ' (mm)'
                elif(var=='Vm'):
                    ylabel += ' (mm/s)'
                elif(var=='Am'):
                    ylabel += ' (mm/s^2)'
                elif(var=='XCE'):
                    ylabel += ' (mm)'
                elif(var=='Cs' or var=='CaSR' or var=='CaSRCS' or var=='B' or var=='CaSP' or var=='T' or var=='CaSPB' or var=='CaSPT'):
                    ylabel += ' (M)'
                self.ax[i].set_ylabel(ylabel)
                
                # new Line
                self.lines[i], = self.ax[i].plot([],[], animated=True)
                self.xdata[i]=[]
                self.ydata[i]=[]
                
            # X label
            if(i == (self.scopeLength-1)):
                self.ax[i].set_xlabel('Time (ms)')
            self.figure.canvas.draw()
            for i in range(self.scopeLength):
                # cache the background
                self.background[i] = self.figure.canvas.copy_from_bbox(self.ax[i].bbox) 

        # Conbined
        elif(display_result == 'Combined'):                
            # add axes
            self.ax = self.figure.add_subplot(1, 1, 1)
            self.ax.set_xlim(self.t_start, self.t_stop)
            self.ax.set_autoscaley_on(True)
            self.ax.grid()
            # X label
            self.ax.set_xlabel('Time (ms)')
            
            for i in range(self.scopeLength):
                var = scope2[i]
                # new Line
                self.lines[i],=self.ax.plot([],[], label=var, animated=True)
                self.xdata[i]=[]
                self.ydata[i]=[]
                
            self.ax.legend(loc='best') 
            self.figure.canvas.draw()
            # cache the background
            self.background = self.figure.canvas.copy_from_bbox(self.ax.bbox)

        # window close button disable
        self.figure.canvas.parent().setWindowFlags(
        Qt.WindowTitleHint | 
        Qt.Dialog | 
        Qt.WindowMinimizeButtonHint | 
        Qt.CustomizeWindowHint)
        self.figure.show()
    
    def updatePlot(self, scope2, display_result, T, ResultArrays2, k, num_steps):
        # Individual
        if(display_result == 'Individual'):      
            for i in range(self.scopeLength):
                var = scope2[i]
                # x, y data update
                self.xdata[i].append(T[k]) 
                self.ydata[i].append(ResultArrays2[var][k])       
                
                # plotting every t_pt step and last k.
                if(k%self.sampling_rate==0 or k==num_steps-1):
                    self.lines[i].set_xdata(self.xdata[i])
                    self.lines[i].set_ydata(self.ydata[i])
                    ymin, ymax = self.ax[i].get_ylim()
                    
                    # scale update
                    if((ResultArrays2[var][k] > ymax) or (ResultArrays2[var][k] < ymin)):
                        self.ax[i].relim()
                        self.ax[i].autoscale_view()
                        self.figure.canvas.draw()
                        
                    # restore background
                    self.figure.canvas.restore_region(self.background[i])
                    self.ax[i].draw_artist(self.lines[i])
                    # fill in the axes rectangle (blit)
                    self.figure.canvas.blit(self.ax[i].bbox)
                    # figure update
                    self.ax[i].relim()
                    self.ax[i].autoscale_view()
                    self.figure.canvas.flush_events()
            
        # Conbined
        elif(display_result == 'Combined'):
            for i in range(self.scopeLength):
                var = scope2[i]
                # x, y data update
                self.xdata[i].append(T[k]) 
                self.ydata[i].append(ResultArrays2[var][k])     
                
                # last plotting for resting part.
                if(k%self.sampling_rate==0 or k==num_steps-1):     
                    self.lines[i].set_xdata(self.xdata[i])
                    self.lines[i].set_ydata(self.ydata[i])
                    ymin, ymax = self.ax.get_ylim()
                    
                    # scale update
                    if((ResultArrays2[var][k] > ymax) or (ResultArrays2[var][k] < ymin)):
                        self.ax.relim()
                        self.ax.autoscale_view()
                        self.figure.canvas.draw()
                    
                    if(i == 0):
                        # restore background
                        self.figure.canvas.restore_region(self.background)
                    self.ax.draw_artist(self.ax.lines[i])

                    if(i == self.scopeLength-1):
                        # fill in the axes rectangle
                        self.figure.canvas.blit(self.ax.bbox)
                        # figure update
                        self.ax.relim()
                        self.ax.autoscale_view()
                        self.figure.canvas.flush_events()
                    
    # Integration function  
    def solModel(self, scope2, display_result):
        # set integrator
        r2= integrate.ode(self.model).set_integrator('lsoda', max_step=10*self.t_dt)
        # initialize integrator
        r2.set_initial_value(self.ivalues, self.t_start) 
        
        # Result arrays
        num_steps = int(np.floor((self.t_stop - self.t_start)/self.t_dt) + 1)  
        T2= np.zeros(num_steps)
        F = np.zeros(num_steps)
        CS= np.zeros(num_steps)
        CaSR=np.zeros(num_steps)
        CaSRCS=np.zeros(num_steps)
        B=np.zeros(num_steps)
        CaSP=np.zeros(num_steps)
        T=np.zeros(num_steps)
        CaSPB=np.zeros(num_steps)
        CaSPT=np.zeros(num_steps)
        A=np.zeros(num_steps)
        XCE=np.zeros(num_steps)
        XM=np.zeros(num_steps)
        VM=np.zeros(num_steps)
        AM=np.zeros(num_steps)
        R=np.zeros(num_steps)
        AS=np.zeros(num_steps)
        SP=np.zeros(num_steps)
        
        # set initial value
        CS[0]=self.ivalues[0]
        CaSR[0]=self.ivalues[1]
        CaSRCS[0]=self.ivalues[2]
        B[0]=self.ivalues[3]
        CaSP[0]=self.ivalues[4]
        T[0]=self.ivalues[5]
        CaSPB[0]=self.ivalues[6]
        CaSPT[0]=self.ivalues[7]
        A[0]=self.ivalues[8]
        XCE[0]=self.ivalues[9]
        SP[0] = self.spike[0]
        XM[0]=self.xm[0]
        VM[0]=0.
        AM[0]=0.
        R[0]=self.get_R(T2[0], CaSR[0])
        F[0]=self.get_F(T2[0], XCE[0], XM[0]) 
        AS[0]=self.get_As(XM[0],A[0],VM[0], AM[0],T[0])

        for i in self.SpikeTimes:
            k = int(round(i / (self.t_dt)))
            SP[k] = 1

        ResultArrays2={ 'Time': T2,
                        'Xm': XM,
                       'Vm'  : VM,
                       'Am' : AM,
                       'R'  : R,
                       'Spike': SP,
                        'F': F,
                        'A': AS,
                       'Cs': CS,
                      'CaSR': CaSR,
                      'CaSRCS': CaSRCS,
                      'B': B,
                      'CaSP': CaSP,
                      'T': T,
                      'CaSPB': CaSPB,
                      'CaSPT': CaSPT,
                      'A_tilde': A,
                      'XCE': XCE}

        k = 0
        # save integration start time
        start_time = time.time()

        # output graph creation
        if(len(scope2) != 0):
            self.createPlot(scope2, display_result)
            self.updatePlot(scope2, display_result, T2, ResultArrays2, k, num_steps)
        else:
            # figure instance creation for prevent GIL
            self.figure=plt.figure(frameon=False)
            
        k = 1
        while r2.successful() and k < num_steps:
            # Integration with muscle fibers
            r2.integrate(round(r2.t, 9) + self.t_dt)
            # integration results 
            T2[k]=r2.t
            CS[k]=r2.y[0]
            CaSR[k]=r2.y[1]
            CaSRCS[k]=r2.y[2]
            B[k]=r2.y[3]
            CaSP[k]=r2.y[4]
            T[k]=r2.y[5]
            CaSPB[k]=r2.y[6]
            CaSPT[k]=r2.y[7]
            A[k]=r2.y[8]
            XCE[k]=r2.y[9]
            
            ## re-calculate results about not state variable 
            # R
            R[k]=self.get_R(T2[k], CaSR[k])

            # Xm, Vm, Am
            XM[k], VM[k], AM[k] =self.get_Xm_Vm_Am(T2[k])
            # A_tilde
            AS[k]=self.get_As(XM[k],A[k],VM[k],AM[k], T[k])
            # F
            F[k]=self.get_F(T2[k], XCE[k], XM[k])
            
            # output plotting
            if(len(scope2) != 0):
                self.updatePlot(scope2, display_result, T2, ResultArrays2, k, num_steps)
            else:
                # GUI event update
                self.figure.canvas.flush_events()

            k += 1
            
            # If user pushes Stop button, cancel the integration
            if(self.cellState=='Stop'):
                # save the results so far
                T2 = T2[:k]
                F = F[:k]
                CS = CS[:k]
                CaSR = CaSR[:k]
                CaSRCS = CaSRCS[:k]
                B = B[:k]
                CaSP = CaSP[:k]
                T = T[:k]
                CaSPB = CaSPB[:k]
                CaSPT = CaSPT[:k]
                A = A[:k]
                XCE = XCE[:k]
                XM = XM[:k]
                VM = VM[:k]
                AM = AM[:k]
                R = R[:k]
                AS = AS[:k]
                SP = SP[:k]
                
                ResultArrays2={ 'Time': T2,
                        'Xm': XM,
                       'Vm'  : VM,
                       'Am' : AM,
                       'R'  : R,
                       'Spike': SP,
                        'F': F,
                        'A': AS,
                       'Cs': CS,
                      'CaSR': CaSR,
                      'CaSRCS': CaSRCS,
                      'B': B,
                      'CaSP': CaSP,
                      'T': T,
                      'CaSPB': CaSPB,
                      'CaSPT': CaSPT,
                      'A_tilde': A,
                      'XCE': XCE
                      }
                     
                if(len(scope2) == 0):
                    self.figure = None
                
                break
        
        # calculate Elapsed time 
        self.simulTime = time.time() - start_time
        
        if(len(scope2) == 0):
            self.figure = None
        else:
            # last graph update
            if(display_result == 'Individual'):
                for i in range(len(scope2)):
                    self.lines[i].set_animated(False)
                    self.background[i] = None
                    self.ax[i].relim()
                    self.ax[i].autoscale_view()
            elif(display_result == 'Combined'):
                for i in range(len(scope2)):
                    self.lines[i].set_animated(False)
                self.background = None
                self.ax.relim()
                self.ax.autoscale_view()
            self.figure.canvas.draw()
            # window close button Enable
            self.figure.canvas.parent().setWindowFlags(
            Qt.Window)
            self.figure.show()

        return ResultArrays2

# Motor Unit class      
class Motorunit:
    def __init__(self, uniqueNumber, MN, MF):
        # initialization
        self.cellType='Motor Unit'
        self.uniqueNumber=uniqueNumber
        self.MN=MN
        self.MF=MF
        self.cellState='Normal'
        self.figure = None
        self.simulTime=0.
        
    def setIntegrationEnv(self, t_start, t_stop, t_dt, t_pt):
        self.t_start=t_start
        self.t_stop=t_stop
        self.t_dt=t_dt
        self.t_pt=t_pt
        
        self.MN.setIntegrationEnv(t_start, t_stop, t_dt, t_pt)
        self.MF.setIntegrationEnv(t_start, t_stop, t_dt, t_pt)
    
    # Integration function    
    def solModel(self, scope1, scope2, display_result):
        # motoneuron constants
        rn,tm,VAsdDC,VAdsDC,VAsdAC,parea,svl,dvl,cv, sf,skca,salpha,sCAo,sEca, sgna,svna,sanamc,sanamv,sanama,sanamb,sbnamc,sbnamv,sbnama,sbnamb,snahth,snahslp,snahv,snaha,snahb,snahc, sgnap,svna,sanapmc,sanapmv,sanapma,sanapmb,sbnapmc,sbnapmv,sbnapma,sbnapmb, sgkdr,svk,skdrth,skdrslp,skdrv,skdra,skdrb,skdrc, sgkca,svk,skd, sgca,scamth,scamslp,scamtau,scahth,scahslp,scahtau, sgh,svh,shth,shslp,shtau, svesyn, svisyn, df,dkca,dalpha,dCAo,dEca, dgcal,dcalth,dcalslp,dcaltau, dgna,dvna,danamc,danamv,danama,danamb,dbnamc,dbnamv,dbnama,dbnamb,dnahth,dnahslp,dnahv,dnaha,dnahb,dnahc, dgnap,dvna,danapmc,danapmv,danapma,danapmb,dbnapmc,dbnapmv,dbnapma,dbnapmb, dgkdr,dvk,dkdrth,dkdrslp,dkdrv,dkdra,dkdrb,dkdrc,dgkca,dvk,dkd, dgca,dcamth,dcamslp,dcamtau,dcahth,dcahslp,dcahtau, dgh,dvh,dhth,dhslp,dhtau, dvesyn, dvisyn=self.MN.parameters
        # initialize array
        self.MF.SpikeTimes=[] # spike detection time array
        self.MN.FiringRate=[] # spike rate array
        const = self.MN.const
        # Isoma heavian parameter
        if(self.MN.IsSignalType=='Step'):
            i0, ip1, pon1, poff1, ip2, pon2, poff2, ip3, pon3, poff3, ip4, pon4, poff4, ip5, pon5, poff5, s = self.MN.heav_param

        ivalues=self.MN.ivalues
        # set MN integrator (vode)
        r1= integrate.ode(self.MN.model).set_integrator('vode')
        # initialize MN integrator
        r1.set_initial_value(ivalues, self.t_start)
        
        # Result arrays
        num_steps = int(np.floor((self.t_stop - self.t_start)/self.t_dt) + 1)
        FiringRate = np.zeros(num_steps)
        T3 = np.zeros(num_steps)
        IS = np.zeros(num_steps)
        VS = np.zeros(num_steps)
        SCA = np.zeros(num_steps)
        SVCA = np.zeros(num_steps)
        SNAI = np.zeros(num_steps)
        SNAM = np.zeros(num_steps)
        SNAH = np.zeros(num_steps)
        SNAPI = np.zeros(num_steps)
        SNAPM = np.zeros(num_steps) 
        SKDRI = np.zeros(num_steps)
        SKDR = np.zeros(num_steps)
        SKCAI = np.zeros(num_steps)
        SCAI = np.zeros(num_steps)
        SCAM = np.zeros(num_steps) 
        SCAH = np.zeros(num_steps)
        SHI = np.zeros(num_steps)
        SHM = np.zeros(num_steps)
        VD = np.zeros(num_steps) 
        DCA = np.zeros(num_steps)
        DVCA = np.zeros(num_steps)
        DCALI = np.zeros(num_steps)
        DCAL = np.zeros(num_steps) 
        DNAI = np.zeros(num_steps)
        DNAM = np.zeros(num_steps)
        DNAH = np.zeros(num_steps)
        DNAPI = np.zeros(num_steps)
        DNAPM = np.zeros(num_steps)
        DKDRI = np.zeros(num_steps)
        DKDR = np.zeros(num_steps)
        DKCAI = np.zeros(num_steps)
        DCAI = np.zeros(num_steps)
        DCAM = np.zeros(num_steps)
        DCAH = np.zeros(num_steps)
        DHI = np.zeros(num_steps)
        DHM = np.zeros(num_steps)
        SESYNI = np.zeros(num_steps)
        SGESYN = np.zeros(num_steps)
        SISYNI = np.zeros(num_steps)
        SGISYN = np.zeros(num_steps)
        DESYNI = np.zeros(num_steps)
        DGESYN = np.zeros(num_steps)
        DISYNI = np.zeros(num_steps)
        DGISYN = np.zeros(num_steps)
            
        # set initial value
        IS[0]=self.MN.Is_0
        VS[0]=ivalues[0]
        SCA[0]=ivalues[1]
        if(self.MN.const_sEca == True):
            SVCA[0] = sEca
        else:
            SVCA[0]=(const*log(sCAo/SCA[0]))-70
        SNAM[0]=ivalues[2]
        SNAH[0]=ivalues[3]
        SNAI[0] = sgna*SNAM[0]**3*SNAH[0]*(VS[0]-svna)
        SNAPM[0]=ivalues[4]
        SNAPI[0] = sgnap*SNAPM[0]**3*(VS[0]-svna)
        SKDR[0]=ivalues[5]
        SKDRI[0] = sgkdr*SKDR[0]**4*(VS[0]-svk)
        SKCAI[0] = sgkca*(SCA[0]/(SCA[0]+skd))*(VS[0]-svk)
        SCAM[0]=ivalues[6]
        SCAH[0]=ivalues[7]
        SCAI[0] = sgca*SCAM[0]**2*SCAH[0]*(VS[0]-SVCA[0])
        SHM[0]=ivalues[8]
        SHI[0] = sgh*SHM[0]*(VS[0]-svh)
        VD[0]=ivalues[9]
        DCA[0]=ivalues[10]
        if(self.MN.const_dEca == True):
            DVCA[0] = dEca
        else:
            DVCA[0]=(const*log(dCAo/DCA[0]))-70
        DCAL[0]=ivalues[11]
        DCALI[0] = dgcal*DCAL[0]*(VD[0]-DVCA[0])
        DNAM[0]=ivalues[12]
        DNAH[0]=ivalues[13]
        DNAI[0] = dgna*DNAM[0]**3*DNAH[0]*(VD[0]-dvna)
        DNAPM[0]=ivalues[14]
        DNAPI[0] = dgnap*DNAPM[0]**3*(VD[0]-dvna)
        DKDR[0]=ivalues[15]
        DKDRI[0] = dgkdr*DKDR[0]**4*(VD[0]-dvk)
        DKCAI[0] = dgkca*(DCA[0]/(DCA[0]+dkd))*(VD[0]-dvk)
        DCAM[0]=ivalues[16]
        DCAH[0]=ivalues[17]
        DCAI[0] = dgca*DCAM[0]**2*DCAH[0]*(VD[0]-DVCA[0])
        DHM[0]=ivalues[18]
        DHI[0] = dgh*DHM[0]*(VD[0]-dvh)
        SGESYN[0] = self.MN.Se_syncon[0]
        SESYNI[0] = SGESYN[0]*(VS[0]-svesyn)
        SGISYN[0] = self.MN.Si_syncon[0]
        SISYNI[0] = SGISYN[0]*(VS[0]-svisyn)
        DGESYN[0] = self.MN.De_syncon[0]
        DESYNI[0] = DGESYN[0]*(VD[0]-dvesyn)
        DGISYN[0] = self.MN.Di_syncon[0]
        DISYNI[0] = DGISYN[0]*(VD[0]-dvisyn)
        
        ResultArrays={'Time': T3,
                      'Is' : IS,
                      'FiringRate' : FiringRate,
                      'V_soma': VS,
                      '[Ca]_soma': SCA,
                      'E_Ca_soma': SVCA,
                      'I_Naf_soma': SNAI,
                      'm_Naf_soma': SNAM,
                      'h_Naf_soma': SNAH,
                      'I_Nap_soma': SNAPI,
                      'm_Nap_soma': SNAPM,
                      'I_Kdr_soma': SKDRI,
                      'n_Kdr_soma' : SKDR,
                      'I_Kca_soma': SKCAI,
                      'I_Can_soma': SCAI,
                      'm_Can_soma': SCAM,
                      'h_Can_soma': SCAH,
                      'I_H_soma': SHI,
                      'm_H_soma' : SHM,    
                      'V_dend': VD,
                      '[Ca]_dend': DCA,
                      'E_Ca_dend': DVCA,
                      'I_Cal_dend': DCALI,
                      'l_Cal_dend': DCAL,
                      'I_Naf_dend': DNAI,
                      'm_Naf_dend': DNAM,
                      'h_Naf_dend': DNAH,
                      'I_Nap_dend': DNAPI,
                      'm_Nap_dend': DNAPM,
                      'I_Kdr_dend': DKDRI,
                      'n_Kdr_dend' : DKDR,
                      'I_Kca_dend': DKCAI,
                      'I_Can_dend': DCAI,
                      'm_Can_dend': DCAM,
                      'h_Can_dend': DCAH,
                      'I_H_dend': DHI,
                      'm_H_dend' : DHM,
                      'I_esyn_soma': SESYNI,
                      'G_esyn_soma': SGESYN,
                      'I_isyn_soma': SISYNI,
                      'G_isyn_soma': SGISYN,
                      'I_esyn_dend': DESYNI,
                      'G_esyn_dend': DGESYN,
                      'I_isyn_dend': DISYNI,
                      'G_isyn_dend': DGISYN }  
        

        ivalues=self.MF.ivalues
        # set MF integrator (lsoda) max_step prevent unpredictable jumping of lsoda.
        r2= integrate.ode(self.MF.model).set_integrator('lsoda', max_step=10*self.t_dt) 
        # initialize MF integrator
        r2.set_initial_value(ivalues, self.t_start)
        
        # Result arrays
        R=np.zeros(num_steps)
        F= np.zeros(num_steps)
        CS= np.zeros(num_steps)
        CaSR=np.zeros(num_steps)
        CaSRCS=np.zeros(num_steps)
        B=np.zeros(num_steps)
        CaSP=np.zeros(num_steps)
        T=np.zeros(num_steps)
        CaSPB=np.zeros(num_steps)
        CaSPT=np.zeros(num_steps)
        A=np.zeros(num_steps)
        XCE=np.zeros(num_steps)
        XM=np.zeros(num_steps)
        VM=np.zeros(num_steps)
        AM=np.zeros(num_steps)
        AS=np.zeros(num_steps)
        SP=np.zeros(num_steps)
        
        # set initial value
        CS[0]=ivalues[0]
        CaSR[0]=ivalues[1]
        CaSRCS[0]=ivalues[2]
        B[0]=ivalues[3]
        CaSP[0]=ivalues[4]
        T[0]=ivalues[5]
        CaSPB[0]=ivalues[6]
        CaSPT[0]=ivalues[7]
        A[0]=ivalues[8]
        XCE[0]=ivalues[9]
        XM[0]=self.MF.xm[0]
        VM[0]=0.
        AM[0]=0.
        F[0]=self.MF.get_F(0, XCE[0], XM[0]) 
        R[0]=self.MF.get_R(0, CaSR[0])
        AS[0]=self.MF.get_As(XM[0],A[0],VM[0], AM[0], T[0])
        
        ResultArrays2={ 'A': AS,
                        'Spike': SP,
                        'Xm': XM,
                       'Vm'  : VM,
                       'Am' : AM,
                       'R'  : R,
                        'F': F,
                       'Cs': CS,
                      'CaSR': CaSR,
                      'CaSRCS': CaSRCS,
                      'B': B,
                      'CaSP': CaSP,
                      'T': T,
                      'CaSPB': CaSPB,
                      'CaSPT': CaSPT,
                      'A_tilde': A,
                      'XCE': XCE}
        ResultArrays.update(ResultArrays2)
    
        k = 0
        # save integration start time
        start_time = time.time()
        
        # output graph creation
        if((len(scope1) != 0) or (len(scope2) != 0)):
            self.createPlot(scope1, scope2, display_result)
            self.updatePlot(T3, display_result, ResultArrays, k, num_steps)

        else:
            # figure instance creation for prevent GIL 
            self.figure=plt.figure(frameon=False)
            
        k = 1
        while r2.successful() and k < num_steps:            
            # Integration with motoneuron
            r1.integrate(round(r1.t, 9) + self.t_dt)
            ## integration results
            T3[k]=r1.t
            VS[k]=r1.y[0]
            SCA[k]=r1.y[1]
            SNAM[k]=r1.y[2]
            SNAH[k]=r1.y[3]
            SNAPM[k]=r1.y[4]
            SKDR[k]=r1.y[5]
            SCAM[k]=r1.y[6]
            SCAH[k]=r1.y[7]
            SHM[k]=r1.y[8]         
            VD[k]=r1.y[9]
            DCA[k]=r1.y[10]
            DCAL[k]=r1.y[11]
            DNAM[k]=r1.y[12]
            DNAH[k]=r1.y[13]
            DNAPM[k]=r1.y[14]
            DKDR[k]=r1.y[15]
            DCAM[k]=r1.y[16]
            DCAH[k]=r1.y[17]
            DHM[k]=r1.y[18]
                        
            ## re-calculate results about not state variable   
            # Isoma
            if(self.MN.IsSignalType=='Step'):
                IS[k] = i0 + s*((self.MN.heav(poff1-T3[k])*self.MN.heav(T3[k]-pon1)*ip1)
                         + (self.MN.heav(poff2-T3[k])*self.MN.heav(T3[k]-pon2)*ip2)
                         + (self.MN.heav(poff3-T3[k])*self.MN.heav(T3[k]-pon3)*ip3)
                         + (self.MN.heav(poff4-T3[k])*self.MN.heav(T3[k]-pon4)*ip4)
                         + (self.MN.heav(poff5-T3[k])*self.MN.heav(T3[k]-pon5)*ip5))
            elif(self.MN.IsSignalType=='Ramp'):
                IS[k] = self.MN.pv-((self.MN.pv-self.MN.iv)/self.MN.p)*abs(T3[k]-self.MN.p)
            elif(self.MN.IsSignalType=='Import'):
                IS[k] = np.interp(T3[k], self.MN.times,  self.MN.Is, 0, 0)
            
            ## [SOMA]
            if (SCA[k] < 1.e-100):
                SCA[k] = 1.e-100
            # Calcium dynamics
            if(self.MN.const_sEca == True):
                SVCA[k] = sEca
            else:
                SVCA[k]=(const*log(sCAo/SCA[k]))-70
            # Naf
            SNAI[k] = sgna*SNAM[k]**3*SNAH[k]*(VS[k]-svna)
            # Nap
            SNAPI[k] = sgnap*SNAPM[k]**3*(VS[k]-svna)
            # Kdr
            SKDRI[k] = sgkdr*SKDR[k]**4*(VS[k]-svk)
            # Kca
            SKCAI[k] = sgkca*(SCA[k]/(SCA[k]+skd))*(VS[k]-svk)
            # Can
            SCAI[k] = sgca*SCAM[k]**2*SCAH[k]*(VS[k]-SVCA[k])
            # H
            SHI[k] = sgh*SHM[k]*(VS[k]-svh)
            # eSyn
            if(self.MN.se_ch == True):
                SGESYN[k] = np.interp(T3[k], self.MN.sesyn_t, self.MN.Se_syncon, 0, 0)
            else:
                SGESYN[k] = 0.
            SESYNI[k] = SGESYN[k]*(VS[k]-svesyn)
            # iSyn
            if(self.MN.si_ch == True):
                SGISYN[k] = np.interp(T3[k], self.MN.sisyn_t, self.MN.Si_syncon, 0, 0)
            else:
                SGISYN[k] = 0.
            SISYNI[k] = SGISYN[k]*(VS[k]-svisyn)
            
            ## [DENDRITE]
            if (DCA[k] < 1.e-100):
                DCA[k] = 1.e-100
            # Calcium dynamics 
            if(self.MN.const_dEca == True):
                DVCA[k] = dEca
            else:
                DVCA[k]=(const*log(dCAo/DCA[k]))-70
            # Cal
            DCALI[k] = dgcal*DCAL[k]*(VD[k]-DVCA[k])
            # Naf
            DNAI[k] = dgna*DNAM[k]**3*DNAH[k]*(VD[k]-dvna)
            # Nap
            DNAPI[k] = dgnap*DNAPM[k]**3*(VD[k]-dvna)
            # Kdr
            DKDRI[k] = dgkdr*DKDR[k]**4*(VD[k]-dvk)
            # Kca
            DKCAI[k] = dgkca*(DCA[k]/(DCA[k]+dkd))*(VD[k]-dvk)
            # Can
            DCAI[k] = dgca*DCAM[k]**2*DCAH[k]*(VD[k]-DVCA[k])
            # H
            DHI[k] = dgh*DHM[k]*(VD[k]-dvh)
            # eSyn
            if(self.MN.de_ch == True):
                DGESYN[k] = np.interp(T3[k], self.MN.desyn_t, self.MN.De_syncon, 0, 0)
            else:
                DGESYN[k] = 0.
            DESYNI[k] = DGESYN[k]*(VD[k]-dvesyn)
            # iSyn            
            if(self.MN.di_ch == True):
                DGISYN[k] = np.interp(T3[k], self.MN.disyn_t, self.MN.Di_syncon, 0, 0)
            else:
                DGISYN[k] = 0.
            DISYNI[k] = DGISYN[k]*(VD[k]-dvisyn)
            
            # Spike detect & calculate Spike rate
            if(k >= 1):
                S_Detect = self.MF.detect_Spike(T3[k], VS[k-2], VS[k-1], VS[k])
            if(S_Detect==True):
                self.MN.cal_FiringRate(self.MF.SpikeTimes)
                SP[k]=1
            
            # Integration with muscle fibers
            r2.integrate(round(r2.t, 9) + self.t_dt)
            T3[k]=r2.t
            CS[k]=r2.y[0]
            CaSR[k]=r2.y[1]
            CaSRCS[k]=r2.y[2]
            B[k]=r2.y[3]
            CaSP[k]=r2.y[4]
            T[k]=r2.y[5]
            CaSPB[k]=r2.y[6]
            CaSPT[k]=r2.y[7]
            A[k]=r2.y[8]
            XCE[k]=r2.y[9]
            
            ## re-calculate results about not state variable
            # R
            R[k]=self.MF.get_R(T3[k], CaSR[k])
            # Xm, Vm, Am
            XM[k], VM[k], AM[k] =self.MF.get_Xm_Vm_Am(T3[k]) 
            # A_tilde
            AS[k]=self.MF.get_As(XM[k],A[k],VM[k], AM[k], T[k])
            # F
            F[k]=self.MF.get_F(T3[k], XCE[k], XM[k])
            
            # output plotting
            if((len(scope1) != 0) or (len(scope2) != 0)):
                self.updatePlot(T3, display_result, ResultArrays, k, num_steps)
            else:
                # GUI event update
                self.figure.canvas.flush_events()
            
            k += 1
            
            # If user pushes Stop button, cancel the integration
            if (self.cellState=='Stop'):
                # save the results so far
                T3 = T3[:k]
                IS = IS[:k]
                VS = VS[:k]
                SCA = SCA[:k]
                SVCA = SVCA[:k]
                SNAI = SNAI[:k]
                SNAM = SNAM[:k]
                SNAH = SNAH[:k]
                SNAPI = SNAPI[:k]
                SNAPM = SNAPM[:k]
                SKDRI = SKDRI[:k]
                SKDR = SKDR[:k]
                SKCAI = SKCAI[:k]
                SCAI = SCAI[:k]
                SCAM = SCAM[:k]
                SCAH = SCAH[:k]
                SHI = SHI[:k]
                SHM = SHM[:k]
                VD = VD[:k]
                DCA = DCA[:k]
                DVCA = DVCA[:k]
                DCALI = DCALI[:k]
                DCAL = DCAL[:k]
                DNAI = DNAI[:k]
                DNAM = DNAM[:k]
                DNAH = DNAH[:k]
                DNAPI = DNAPI[:k]
                DNAPM = DNAPM[:k]
                DKDRI = DKDRI[:k]
                DKDR = DKDR[:k]
                DKCAI = DKCAI[:k]
                DCAI = DCAI[:k]
                DCAM = DCAM[:k]
                DCAH = DCAH[:k]
                DHI = DHI[:k]
                DHM = DHM[:k]
                SESYNI = SESYNI[:k]
                SGESYN = SGESYN[:k]
                SISYNI = SISYNI[:k]
                SGISYN = SGISYN[:k]
                DESYNI = DESYNI[:k]
                DGESYN = DGESYN[:k]
                DISYNI = DISYNI[:k]
                DGISYN = DGISYN[:k]
                F = F[:k]
                CS = CS[:k]
                CaSR = CaSR[:k]
                CaSRCS = CaSRCS[:k]
                B = B[:k]
                CaSP = CaSP[:k]
                T = T[:k]
                CaSPB = CaSPB[:k]
                CaSPT = CaSPT[:k]
                A = A[:k]
                XCE = XCE[:k]
                XM = XM[:k]
                VM = VM[:k]
                AM = AM[:k]
                R = R[:k]
                AS = AS[:k]
                SP = SP[:k]
                
                ResultArrays={'Time': T3,
                      'Is' : IS,
                      'V_soma': VS,
                      '[Ca]_soma': SCA,
                      'E_Ca_soma': SVCA,
                      'I_Naf_soma': SNAI,
                      'm_Naf_soma': SNAM,
                      'h_Naf_soma': SNAH,
                      'I_Nap_soma': SNAPI,
                      'm_Nap_soma': SNAPM,
                      'I_Kdr_soma': SKDRI,
                      'n_Kdr_soma' : SKDR,
                      'I_Kca_soma': SKCAI,
                      'I_Can_soma': SCAI,
                      'm_Can_soma': SCAM,
                      'h_Can_soma': SCAH,
                      'I_H_soma': SHI,
                      'm_H_soma' : SHM,    
                      'V_dend': VD,
                      '[Ca]_dend': DCA,
                      'E_Ca_dend': DVCA,
                      'I_Cal_dend': DCALI,
                      'l_Cal_dend': DCAL,
                      'I_Naf_dend': DNAI,
                      'm_Naf_dend': DNAM,
                      'h_Naf_dend': DNAH,
                      'I_Nap_dend': DNAPI,
                      'm_Nap_dend': DNAPM,
                      'I_Kdr_dend': DKDRI,
                      'n_Kdr_dend' : DKDR,
                      'I_Kca_dend': DKCAI,
                      'I_Can_dend': DCAI,
                      'm_Can_dend': DCAM,
                      'h_Can_dend': DCAH,
                      'I_H_dend': DHI,
                      'm_H_dend' : DHM,
                      'I_esyn_soma': SESYNI,
                      'G_esyn_soma': SGESYN,
                      'I_isyn_soma': SISYNI,
                      'G_isyn_soma': SGISYN,
                      'I_esyn_dend': DESYNI,
                      'G_esyn_dend': DGESYN,
                      'I_isyn_dend': DISYNI,
                      'G_isyn_dend': DGISYN,
                      'Xm': XM,
                      'Vm'  : VM,
                      'Am' : AM,
                      'R'  : R,
                      'Spike': SP,
                      'F': F,
                      'A': AS,
                      'Cs': CS,
                      'CaSR': CaSR,
                      'CaSRCS': CaSRCS,
                      'B': B,
                      'CaSP': CaSP,
                      'T': T,
                      'CaSPB': CaSPB,
                      'CaSPT': CaSPT,
                      'A_tilde': A,
                      'XCE': XCE
                      }
                      
                if(len(scope1) == 0 and len(scope2) == 0):
                    self.figure = None
                break
        
        k = 0
        for i in self.MF.SpikeTimes[1:]:
            FiringRate[int(round(i/self.t_dt, 9))] = self.MN.FiringRate[k]
            k+=1
        
        # calculate Elapsed time      
        self.simulTime=time.time() - start_time

        if(len(scope1) == 0 and len(scope2) == 0):
            self.figure = None
        
        else:
            # last graph update
            if(display_result == 'Individual'):
                for i in range(self.scopeLength):
                    self.lines[i].set_animated(False)
                    self.background[i] = None
                    self.ax[i].relim()
                    self.ax[i].autoscale_view()
            elif(display_result == 'Combined'):
                for i in range(self.scopeLength):
                    self.lines[i].set_animated(False)
                self.background = None
                self.ax.relim()
                self.ax.autoscale_view()
            self.figure.canvas.draw()
            # window close button Enable
            self.figure.canvas.parent().setWindowFlags(
            Qt.Window)
            self.figure.show()
        
        return ResultArrays
    
    def createPlot(self, scope1, scope2, display_result):
        scopeLength1=len(scope1)
        scopeLength2=len(scope2)        
        self.scopeLength = scopeLength1 + scopeLength2
        self.lines=range(self.scopeLength)
        self.xdata=range(self.scopeLength)
        self.ydata=range(self.scopeLength)
        self.scope = scope1 + scope2
        self.background=range(self.scopeLength)
        self.sampling_rate=self.t_pt/self.t_dt
        plt.rc('figure', figsize=(7, 7))
        
        # close previous figure
        if(self.figure != None):
            plt.close(self.figure)
            
        # new figure
        self.figure = plt.figure()
        self.figure.canvas.set_window_title('MOTOR UNIT')
                
        # Individual
        if(display_result == 'Individual'):
            self.ax=range(self.scopeLength)
            
            for i in range(self.scopeLength):
                var = self.scope[i]
                # add axes
                self.ax[i] = self.figure.add_subplot(self.scopeLength, 1, i+1)
                self.ax[i].set_xlim(self.t_start, self.t_stop)
                self.ax[i].set_autoscaley_on(True)
                self.ax[i].grid()
                # Y label
                ylabel = var
                if(var=='V_soma' or var=='V_dend' or var=='E_Ca_soma' or var=='E_Ca_dend'):
                    ylabel += ' (mV)'
                elif(var=='Firing_rate'):
                    ylabel += ' (Hz)'
                elif(var=='Is'):
                    ylabel += ' (nA)'
                elif(var=='G_esyn_soma' or var=='G_isyn_soma' or var=='G_esyn_dend' or var=='G_isyn_dend'):
                    ylabel +=' (mS/cm^2)'
                elif(var=='[Ca]_soma' or var=='[Ca]_dend'):
                    ylabel += ' (mM)'
                elif(var=='I_Naf_soma' or var=='I_Nap_soma' or var=='I_Kdr_soma' or var=='I_Kca_soma' or var=='I_Can_soma' or var=='I_H_soma' or var=='I_esyn_soma' or var=='I_isyn_soma' or var=='I_Cal_dend' or var=='I_Naf_dend' or var=='I_Nap_dend' or var=='I_Kdr_dend' or var=='I_Kca_dend' or var=='I_Can_dend' or var=='I_H_dend' or var=='I_esyn_dend' or var=='I_isyn_dend'):
                    ylabel += ' (mA/cm^2)'
                elif(var=='F'):
                    ylabel += ' (N)'
                elif(var=='Xm'):
                    ylabel += ' (mm)'
                elif(var=='Vm'):
                    ylabel += ' (mm/s)'
                elif(var=='Am'):
                    ylabel += ' (mm/s^2)'
                elif(var=='XCE'):
                    ylabel += ' (mm)'
                elif(var=='Cs' or var=='CaSR' or var=='CaSRCS' or var=='B' or var=='CaSP' or var=='T' or var=='CaSPB' or var=='CaSPT'):
                    ylabel += ' (M)'
                self.ax[i].set_ylabel(ylabel)
                
                # new Line with dot
                if(var=='Firing_rate'):
                    self.lines[i], = self.ax[i].plot([],[], '.', animated=True)
                # new Line with line
                else:
                    self.lines[i], = self.ax[i].plot([],[], animated=True)
                self.xdata[i]=[]
                self.ydata[i]=[]
            
            # X label
            if(i == (self.scopeLength-1)):
                self.ax[i].set_xlabel('Time (ms)')
            self.figure.canvas.draw()

            for i in range(self.scopeLength):
                # cache the background
                self.background[i] = self.figure.canvas.copy_from_bbox(self.ax[i].bbox)
               
        # Conbined
        elif(display_result == 'Combined'):                
            # add axes
            self.ax = self.figure.add_subplot(1, 1, 1)
            self.ax.set_xlim(self.t_start, self.t_stop)
            self.ax.set_autoscaley_on(True)
            self.ax.grid()
            # X label
            self.ax.set_xlabel('Time (ms)')
            
            for i in range(self.scopeLength):
                var = self.scope[i]
                # new Line with dot
                if(var=='Firing_rate'):
                    self.lines[i], = self.ax.plot([],[], '.', label=var, animated=True)
                # new Line with line
                else:
                    self.lines[i],=self.ax.plot([],[], label=var, animated=True)
                self.xdata[i]=[]
                self.ydata[i]=[]
                
                ## Y label
                # one variable
                if(self.scopeLength == 1):           
                    if(var=='V_soma' or var=='V_dend' or var=='E_Ca_soma' or var=='E_Ca_dend'):
                        ylabel = ' (mV)'
                    elif(var=='Firing_rate'):
                        ylabel = ' (Hz)'
                    elif(var=='Is'):
                        ylabel = ' (nA)'
                    elif(var=='G_esyn_soma' or var=='G_isyn_soma' or var=='G_esyn_dend' or var=='G_isyn_dend'):
                        ylabel =' (mS/cm^2)'
                    elif(var=='[Ca]_soma' or var=='[Ca]_dend'):
                        ylabel = ' (mM)'
                    elif(var=='I_Naf_soma' or var=='I_Nap_soma' or var=='I_Kdr_soma' or var=='I_Kca_soma' or var=='I_Can_soma' or var=='I_H_soma' or var=='I_esyn_soma' or var=='I_isyn_soma' or var=='I_Cal_dend' or var=='I_Naf_dend' or var=='I_Nap_dend' or var=='I_Kdr_dend' or var=='I_Kca_dend' or var=='I_Can_dend' or var=='I_H_dend' or var=='I_esyn_dend' or var=='I_isyn_dend'):
                        ylabel = ' (mA/cm^2)'
                    elif(var=='F'):
                        ylabel = ' (N)'
                    elif(var=='Xm'):
                        ylabel = ' (mm)'
                    elif(var=='Vm'):
                        ylabel = ' (mm/s)'
                    elif(var=='Am'):
                        ylabel = ' (mm/s^2)'
                    elif(var=='XCE'):
                        ylabel = ' (mm)'
                    elif(var=='Cs' or var=='CaSR' or var=='CaSRCS' or var=='B' or var=='CaSP' or var=='T' or var=='CaSPB' or var=='CaSPT'):
                        ylabel = ' (M)'
                        self.ax.set_ylabel(ylabel)
            # two variable  
            if(self.scopeLength == 2):
                if((self.scope[0]=='V_soma' and self.scope[1] == 'V_dend') or (self.scope[1]=='V_soma' and self.scope[0] == 'V_dend')):                     
                    ylabel = 'V_soma & V_dend (mV)'
                    self.ax.set_ylabel(ylabel)
                    
            self.ax.legend(loc='best')
            self.figure.canvas.draw()
           # cache the background
            self.background = self.figure.canvas.copy_from_bbox(self.ax.bbox)
        
        # window close button disable
        self.figure.canvas.parent().setWindowFlags(
        Qt.WindowTitleHint | 
        Qt.Dialog | 
        Qt.WindowMinimizeButtonHint| 
        Qt.CustomizeWindowHint)
        self.figure.show()  
    
    def updatePlot(self, T, display_result, ResultArrays, k, num_steps):
        # Individual
        if(display_result == 'Individual'):    
            for i in range(self.scopeLength):
                var = self.scope[i]
                # x, y data update
                if(var=='Firing_rate'):
                    self.xdata[i] = self.MF.SpikeTimes[1:] # except first value
                    self.ydata[i] = self.MN.FiringRate
                else:
                    self.xdata[i].append(T[k])
                    self.ydata[i].append(ResultArrays[var][k]) 
                
                # plotting every t_pt step and last k.
                if(k%self.sampling_rate==0 or k==num_steps-1): 
                    self.lines[i].set_xdata(self.xdata[i])
                    self.lines[i].set_ydata(self.ydata[i])
                    ymin, ymax = self.ax[i].get_ylim()
                    
                    # scale update
                    if(var=='Firing_rate'):
                        if(len(self.ydata[i]) > 1):
                            if((self.ydata[i][-1] > ymax) or (self.ydata[i][-1] < ymin)):
                                self.ax[i].relim()
                                self.ax[i].autoscale_view()
                                self.figure.canvas.draw()
                    else:
                        if((ResultArrays[var][k] > ymax) or (ResultArrays[var][k] < ymin)):
                            self.ax[i].relim()
                            self.ax[i].autoscale_view()
                            self.figure.canvas.draw()
                    
                    # restore background
                    self.figure.canvas.restore_region(self.background[i])
                    self.ax[i].draw_artist(self.lines[i])
                    # fill in the axes rectangle (blit)
                    self.figure.canvas.blit(self.ax[i].bbox)
                    ### figure update
                    self.ax[i].relim()
                    self.ax[i].autoscale_view()
                    self.figure.canvas.flush_events()
            
        # Conbined
        elif(display_result == 'Combined'):                   
            for i in range(self.scopeLength):
                var = self.scope[i]
                # x, y data update
                if(var=='Firing_rate'):
                    self.xdata[i] = self.MF.SpikeTimes[1:]
                    self.ydata[i] = self.MN.FiringRate
                else:
                    self.xdata[i].append(T[k]) 
                    self.ydata[i].append(ResultArrays[var][k])     

                # last plotting for resting part.
                if(k%self.sampling_rate==0 or k==num_steps-1):
                    self.lines[i].set_xdata(self.xdata[i])
                    self.lines[i].set_ydata(self.ydata[i])
                    ymin, ymax = self.ax.get_ylim()
                
                    # scale update
                    if(var=='Firing_rate'):
                        if(len(self.ydata[i]) > 1):
                            if((self.ydata[i][-1] > ymax) or (self.ydata[i][-1] < ymin)):
                                self.ax.relim()
                                self.ax.autoscale_view()
                                self.figure.canvas.draw()
                    else:
                        if((ResultArrays[var][k] > ymax) or (ResultArrays[var][k] < ymin)):
                            self.ax.relim()
                            self.ax.autoscale_view()
                            self.figure.canvas.draw()
                    
                    if(i == 0):
                        # restore background
                        self.figure.canvas.restore_region(self.background)
                    self.ax.draw_artist(self.ax.lines[i])

                    if(i == self.scopeLength-1):
                        # fill in the axes rectangle (blit)
                        self.figure.canvas.blit(self.ax.bbox)
                        # figure update
                        self.ax.relim()
                        self.ax.autoscale_view()
                        self.figure.canvas.flush_events()
        
    def setInputSignal(self, IsSignalType, IsiValue, IspValue, Is_0, heav_param, IsPeriod=0, times=0, Is=0):
        self.MN.setInputSignal(IsSignalType, IsiValue, IspValue, Is_0, heav_param, IsPeriod, times, Is)
        
    def setSynConSignal(self, se_ch, si_ch, de_ch, di_ch, Se_times, Si_times, De_times, Di_times, Se_syncon, Si_syncon, De_syncon, Di_syncon):
        self.MN.setSynConSignal(se_ch, si_ch, de_ch, di_ch, Se_times, Si_times, De_times, Di_times, Se_syncon, Si_syncon, De_syncon, Di_syncon)
        
    def setXmSignal(self, XmSignalType, times, Xm):
        self.MF.setXmSignal(XmSignalType, times, Xm)
    
    def setModelParam(self, param1, param2, const_sEca, const_dEca): 
        self.MN.setModelParam(param1, const_sEca, const_dEca)
        self.MF.setModelParam(param2)
        cv = param1[0][8] # axon conduction velocity
        ms = 1000
        self.MF.sd = 1/cv *ms # spike delay (ms) = axon length 1m / cv (m/s) * ms
        
    def setInitialValues(self, ivalues1,ivalues2):
        self.MN.setInitialValues(ivalues1)
        self.MF.setInitialValues(ivalues2)
            
# Isoma signal generator class       
class InputSignalGenerator:
    def __init__(self, signalType, t_final, t_dt, heav_param=[], iValue=0, pValue=0, period=0, exp_time=0, exp_Is=0):
        # initialization
        self.signalType=signalType
        self.t_final=t_final
        self.t_dt=t_dt
        self.iValue=iValue
        self.pValue=pValue
        self.period=period
        self.heav_param=heav_param
        self.Is = []
        self.times = []
        self.exp_time=exp_time
        self.exp_Is=exp_Is
        
    def setValue(self, signalType, t_final, t_dt, heav_param=[], iValue=0, pValue=0, period=0, exp_time=0, exp_Is=0):
        # initialization
        self.signalType=signalType
        self.t_final=t_final
        self.t_dt=t_dt
        self.iValue=iValue
        self.pValue=pValue
        self.period=period
        self.heav_param=heav_param
        self.exp_time=exp_time
        self.exp_Is=exp_Is
    
    def heav(self, x):
        # heavian function
        return np.where(x < 0, 0, 1)
        
    def genSignal(self):
        num_steps = int(np.floor((self.t_final/self.t_dt)+1))
        t = np.linspace(0, self.t_final, num_steps)
        self.times = t
        Is = np.zeros(num_steps)
        
        # generate Isoma signal
        if(self.signalType=='Step'):
            i0, ip1, pon1, poff1, ip2, pon2, poff2, ip3, pon3, poff3, ip4, pon4, poff4, ip5, pon5, poff5, s = self.heav_param
            Is = i0 + s*((self.heav(poff1-t)*self.heav(t-pon1)*ip1)
                         + (self.heav(poff2-t)*self.heav(t-pon2)*ip2)
                         + (self.heav(poff3-t)*self.heav(t-pon3)*ip3)
                         + (self.heav(poff4-t)*self.heav(t-pon4)*ip4)
                         + (self.heav(poff5-t)*self.heav(t-pon5)*ip5))
        elif(self.signalType=='Ramp'):
            iv=self.iValue
            pv=self.pValue
            p=self.period
            # raise error if the value is '0' or negative number
            if(np.sign(p) != 1):
                raise Exception
            Is = pv-((pv-iv)/p)*np.abs(t-p)
        elif(self.signalType=='Import'):
            # raise Error if data is empty.
            if(self.exp_Is == []):
                raise Exception
            Is = self.exp_Is
            self.times= self.exp_time
            for i in self.times: 
                # raise Error if time is negative number.
                if(np.sign(i) == -1):
                        raise Exception
        
        self.Is = Is # generated Isoma signal
        self.Is_0 = Is[0] # Isoma initial value

# Simulation control class
class Oscilloscope():
    def __init__(self):
        # MN output list
        self.mn_scopeList=[]
        # MF output list
        self.mf_scopeList=[]
        # Model
        self.cell = None
        # array for simulation results
        self.mn_ResultArrays = None
        self.mf_ResultArrays = None
        self.mu_ResultArrays = None
    
    def setModel(self, cell):
        self.cell = cell
    
    def setScope(self, mn_scopeList=[], mf_scopeList=[], display_result='Individual'):
        # set the output list
        self.mn_scopeList=mn_scopeList
        self.mf_scopeList=mf_scopeList
        self.display_result = display_result

    def run(self):
        # execute the integration function
        if(self.cell.cellType=='Motoneuron'):
            self.mn_ResultArrays = self.cell.solModel(self.mn_scopeList, self.display_result)
        elif(self.cell.cellType=='Muscle Fibers'):
            self.mf_ResultArrays = self.cell.solModel(self.mf_scopeList, self.display_result)
        elif(self.cell.cellType=='Motor Unit'):            
            self.mu_ResultArrays = self.cell.solModel(self.mn_scopeList, self.mf_scopeList, self.display_result)

# Xm signal generator class
class XmSignalGenerator:
    def __init__(self, signalType, t_final=0, t_dt=0, value1=0, value2=0, time1=0, time2=0, exp_time=0, exp_Xm=0):
        # initialization
        self.signalType=signalType
        self.t_final=t_final
        self.t_dt=t_dt
        self.value1=value1
        self.value2=value2
        self.time1=time1 
        self.time2=time2
        self.exp_time=exp_time
        self.exp_Xm=exp_Xm
        self.xm=[]
        self.times=[]
            
    def setValue(self, signalType, t_final=0, t_dt=0, value1=0, value2=0, time1=0, time2=0, exp_time=0, exp_Xm=0):
        # initialization
        self.signalType=signalType
        self.t_final=t_final
        self.t_dt=t_dt
        self.value1=value1
        self.value2=value2
        self.time1=time1 
        self.time2=time2
        self.exp_time=exp_time
        self.exp_Xm=exp_Xm
            
    def genSignal(self):
        num_steps = int(np.floor(self.t_final/self.t_dt)+1)
        times = np.linspace(0, self.t_final, num_steps)
        self.times = times
        self.xm = np.zeros(num_steps)
        i = 0
            
        # generate Xm signal
        if(self.signalType=='Isometric'):
            for t in times:
                self.xm[i] = self.value1
                i += 1
        elif(self.signalType=='Isokinetic'):
            # raise Error about time
            if(self.time1 >= self.time2 or np.sign(self.time1) == -1 or np.sign(self.time2) == -1):
                raise Exception
            
            vm = (self.value2-self.value1)/(self.time2-self.time1)
            
            for t in times:
                if(t <= self.time1):
                    self.xm[i] = self.value1
                elif(t > self.time1 and t < self.time2):
                    self.xm[i] = vm*(t-self.time1)+self.value1
                else:
                    self.xm[i] = self.value2
                i += 1
        elif(self.signalType=='Exp'):
            # raise Error if data is empty.
            if(self.exp_Xm == []):
                raise Exception
            self.xm = self.exp_Xm
            self.times= self.exp_time
            for i in self.times: 
                # raise Error if time is negative number.
                if(np.sign(i) == -1):
                        raise Exception
        elif(self.signalType=='Dynamic'):
            ms=0.001
            sample_rate=10000 # sampling rate 0.0001 sec
            t_dt=1/float(sample_rate)
            
            filter_order=3000*(self.t_final)*ms
            # Length of the filter
            numtaps=filter_order+1
            # The first N-1 samples are "corrupted" by the initial conditions
            warmup = filter_order/2
            # The phase delay of the filtered signal
            delay = (warmup) / (sample_rate*ms)
            
            t_start=0
            t_stop=self.t_final+(delay) # msec
            num=int(((t_stop - t_start)*ms)/t_dt)+1
            times=np.linspace(t_start,t_stop,num,endpoint=True)
            # white uniform 10khz sequence using a random number generator
            n=25
            low=8*n # mm
            high=-8*n # mm
            y=np.random.uniform(low, high, num)
            
            # The Nyquist rate of the signal.
            nyq_rate = sample_rate / 2.
            # The cutoff frequency of the filter
            cutoff_hz = 5
            # Use firwin with a blackman window to create a lowpass FIR filter.
            fir_coeff =signal.firwin(numtaps, cutoff_hz/nyq_rate, window='blackman')
            # Use lfilter to filter x with the FIR filter.
            filtered_signal = signal.lfilter(fir_coeff, 1.0, y)
            # adjust filtered signal
            xm=filtered_signal[int(warmup):]-8
            times_xm=times[int(warmup):]-delay
            
            self.times = times_xm
            self.xm = xm # generated Xm signal

            return
            
# Iaxon signal generator class
class SpikeSignalGenerator:
    def __init__(self, signalType, t_final, t_dt, time1=0, time2=0, hz=0, scale=0, exp_SpikeTimes=0):
        # initialization        
        self.signalType=signalType
        self.time1=time1 
        self.time2=time2
        self.hz=hz
        self.scale=scale
        self.t_final=t_final
        self.t_dt=t_dt
        self.exp_SpikeTimes = exp_SpikeTimes    
        self.spike = []
        self.SpikeTimes = []
        self.spike_idx = []
        
    def setValue(self, signalType, t_final, t_dt, time1=0, time2=0, hz=0, scale=0, exp_SpikeTimes=0):
        # initialization        
        self.signalType=signalType
        self.time1=time1
        self.time2=time2
        self.hz=hz
        self.scale=scale
        self.t_final=t_final
        self.t_dt=t_dt
        self.exp_SpikeTimes = exp_SpikeTimes
        
    def genSignal(self):
        # generate Iaxon signal
        if(self.signalType == 'User'):
            # raise Error about time
            if(self.time1 >= self.time2 or np.sign(self.time1) == -1 or np.sign(self.time2) == -1):
                raise Exception
                
            num_steps = int(np.floor((self.t_final)/self.t_dt)+1)  
            t = np.linspace(0, self.t_final, num_steps)
            
            self.spike_idx = t
            self.spike = np.zeros(num_steps)
            if(self.hz == 0):
                self.SpikeTimes = []
                return
            ms = 1000
            t1 = self.time1
            t2 = self.time2
            period = 1/self.hz*ms
            times = np.arange(t1, t2, period)
            
            if(t2 == times[-1]+period):
                times = np.append(times, t2)

            # create random spike from times
            num = len(times)
            rn = np.zeros(num)
            k = 0
            if(self.scale != 0):
                for t in times:
                    if(k == 0):
                        rn[k] = times[k]
                    else:
                        rn[k] = np.random.normal(t, self.scale, 1)
                    k += 1
            else:
                rn = times
        
            rn = np.around(rn, decimals = 1)
            rn.sort()
            self.SpikeTimes = rn # spike times array

            # fit to simulation time
            for i in rn:
                k = int(round(i / (self.t_dt)))
                if (k > num_steps):
                    break
                else:
                    self.spike[k] = 1   
        
        elif(self.signalType == 'Exp'):
            # raise Error if data is empty.
            if(self.exp_SpikeTimes == []):
                raise Exception
            num_steps = int(np.floor((self.t_final)/self.t_dt)+1)  
            t = np.linspace(0, self.t_final, num_steps)
            self.spike_idx = t
            self.spike = np.zeros(num_steps)
            
            # fit to simulation time
            for i in self.exp_SpikeTimes: 
                if(np.sign(i) == -1):
                    raise Exception
                k = int(round(i/(self.t_dt)))
                if (k > num_steps):
                    break
                else:
                    self.spike[k] = 1
                    
            self.SpikeTimes = self.exp_SpikeTimes
            
# synaptic conductance generator class            
class SynConSignalGenerator:
    def __init__(self, synType, signalType, t_final, t_dt, heav_param=[], iValue=0, pValue=0, period=0, tau=0., std_max=0, noise=False, exp_time=0, exp_G=0):  
        # initialization
        self.synType=synType
        self.t_final=t_final
        self.t_dt=t_dt
        self.e_times = []
        self.i_times = []
        self.G_e = []
        self.G_i = []     
        if(self.synType=='Excitatory'):
            self.e_signalType=signalType
            self.e_iValue=iValue
            self.e_pValue=pValue
            self.e_heav_param=heav_param
            self.e_period=period
            self.e_tau=tau
            self.e_std_max=std_max
            self.e_noise=noise
            self.e_exp_time = exp_time
            self.e_exp_G = exp_G
        elif(self.synType=='Inhibitory'):
            self.i_signalType=signalType
            self.i_iValue=iValue
            self.i_pValue=pValue
            self.i_heav_param=heav_param
            self.i_period=period
            self.i_tau=tau
            self.i_std_max=std_max
            self.i_noise=noise
            self.i_exp_time = exp_time
            self.i_exp_G = exp_G

    def setValue(self, synType, signalType, t_final, t_dt, heav_param=[], iValue=0, pValue=0, period=0, tau=0., std_max=0, noise=False, exp_time=0, exp_G=0):
        # initialization
        self.synType=synType
        self.t_final=t_final
        self.t_dt=t_dt
        if(self.synType=='Excitatory'):
            self.e_signalType=signalType
            self.e_iValue=iValue
            self.e_pValue=pValue
            self.e_heav_param=heav_param
            self.e_period=period
            self.e_tau=tau
            self.e_std_max=std_max
            self.e_noise=noise
            self.e_exp_time = exp_time
            self.e_exp_G = exp_G
        elif(self.synType=='Inhibitory'):
            self.i_signalType=signalType
            self.i_iValue=iValue
            self.i_pValue=pValue
            self.i_heav_param=heav_param
            self.i_period=period
            self.i_tau=tau
            self.i_std_max=std_max
            self.i_noise=noise
            self.i_exp_time = exp_time
            self.i_exp_G = exp_G

    def heav(self, x):
        # heavian function
        return np.where(x < 0, 0, 1)
        
    def genSignal(self):
        dt = self.t_dt
        num_steps = int(np.floor(self.t_final/self.t_dt+1))
        t = np.linspace(0, self.t_final, num_steps)
        
        # Excitatory synapse channel
        if(self.synType=='Excitatory'):
            G_e = np.zeros(num_steps)
            self.e_times = t
            
            if(self.e_signalType=='Import'):
                G_e = self.e_exp_G
                self.e_times= self.e_exp_time
                for i in self.e_times: 
                    # raise Error if time is negative number.
                    if(np.sign(i) == -1):
                            raise Exception
                # negative values are not allowed.
                for i in range(len(self.e_times)):
                    if(G_e[i] < 0.):
                        G_e[i] = 0.
                self.G_e = G_e
                return
            
            ## generate Isyn signal
            # No noisy signal (linear signal)  
            if(self.e_noise == False):
                if(self.e_signalType=='Step'):    
                    i0, ip1, pon1, poff1, ip2, pon2, poff2, ip3, pon3, poff3, ip4, pon4, poff4, ip5, pon5, poff5, s = self.e_heav_param
                    G_e = i0 + s*((self.heav(poff1-t)*self.heav(t-pon1)*ip1)
                            + (self.heav(poff2-t)*self.heav(t-pon2)*ip2)
                            + (self.heav(poff3-t)*self.heav(t-pon3)*ip3)
                            + (self.heav(poff4-t)*self.heav(t-pon4)*ip4)
                            + (self.heav(poff5-t)*self.heav(t-pon5)*ip5))
                elif(self.e_signalType=='Ramp'):
                    iv=self.e_iValue
                    pv=self.e_pValue
                    p=self.e_period
                    G_e = pv-((pv-iv)/p)*np.abs(t-p)

            # noisy signal (Ornstein–Uhlenbeck process - exact update rule)
            else:
                G_e1 = np.zeros(num_steps)
                amp_e = np.zeros(num_steps)
                tau_e = self.e_tau
                std_e = self.e_std_max
                
                if(tau_e !=0):
                    if(self.e_signalType=='Step'):
                        i0, ip1, pon1, poff1, ip2, pon2, poff2, ip3, pon3, poff3, ip4, pon4, poff4, ip5, pon5, poff5, s = self.e_heav_param
                        G_e0 = max([i0, ip1, ip2, ip3, ip4, ip5])
                        # calculate normalized value (0~1)
                        norm_e = i0 + s*((self.heav(poff1-t)*self.heav(t-pon1)*ip1)
                                + (self.heav(poff2-t)*self.heav(t-pon2)*ip2)
                                + (self.heav(poff3-t)*self.heav(t-pon3)*ip3)
                                + (self.heav(poff4-t)*self.heav(t-pon4)*ip4)
                                + (self.heav(poff5-t)*self.heav(t-pon5)*ip5))
                        norm_e /= G_e0
                    
                    elif(self.e_signalType=='Ramp'):
                        norm_e = np.zeros(num_steps)
                        pv=self.e_pValue
                        iv = self.e_iValue
                        p = self.e_period
                        ## calculate normalized value (0~1)
                        # Ramp
                        if(pv >= iv):
                            G_e0 = pv
                            if(pv == 0):
                                norm_e = np.ones(num_steps)
                            else:
                                x = iv/pv
                                norm_e =1-((1-x)/p)*np.abs(t-p) # x->1->x
                        # inverse Ramp
                        else:
                            G_e0 = iv
                            x = pv/iv
                            norm_e =x-((x-1)/p)*np.abs(t-p) # 1->x->1
                    
                    # negative values are not allowed.
                    for i in xrange(0, num_steps):
                        if(norm_e[i] < 0):
                            norm_e[i] = 0.
                                
                    exp_e = exp(-dt/tau_e)
                    amp_e = np.sqrt(norm_e)*std_e * np.sqrt( (1-np.exp(-2*dt/tau_e)) )
                    
                    # exact update rule
                    for i in xrange(1, num_steps):
                        G_e1[i] = exp_e * G_e1[i - 1] + amp_e[i - 1] * np.random.normal(loc=0.0, scale=1.0)
                    G_e = norm_e * G_e0 + G_e1
                # White noise
                else:
                    for i in xrange(0, num_steps):
                        G_e[i] = std_e * np.random.normal(loc=0.0, scale=1.0)
            
            # negative values are not allowed.
            for i in xrange(0, num_steps):
                if(G_e[i] < 0):
                    G_e[i] = 0.
            
            self.G_e = G_e
                
        # Inhibitory synapse channel
        if(self.synType=='Inhibitory'):
            G_i = np.zeros(num_steps)
            self.i_times = t
            
            if(self.i_signalType=='Import'):
                G_i = self.i_exp_G
                self.i_times= self.i_exp_time
                for i in self.i_times: 
                    # raise Error if time is negative number.
                    if(np.sign(i) == -1):
                            raise Exception
                # negative values are not allowed.
                for i in range(len(self.i_times)):
                    if(G_i[i] < 0):
                        G_i[i] = 0.
                self.G_i = G_i
                return
            
            ## generate Isyn signal
            # No noisy signal (linear signal)        
            if(self.i_noise == False):
                if(self.i_signalType=='Step'):
                    i0, ip1, pon1, poff1, ip2, pon2, poff2, ip3, pon3, poff3, ip4, pon4, poff4, ip5, pon5, poff5, s = self.i_heav_param
                    G_i = i0 + s*((self.heav(poff1-t)*self.heav(t-pon1)*ip1)
                            + (self.heav(poff2-t)*self.heav(t-pon2)*ip2)
                            + (self.heav(poff3-t)*self.heav(t-pon3)*ip3)
                            + (self.heav(poff4-t)*self.heav(t-pon4)*ip4)
                            + (self.heav(poff5-t)*self.heav(t-pon5)*ip5))
                elif(self.i_signalType=='Ramp'):
                    iv=self.i_iValue
                    pv=self.i_pValue
                    p=self.i_period
                    G_i = pv-((pv-iv)/p)*np.abs(t-p)

            # noisy signal (Ornstein–Uhlenbeck process - exact update rule)
            else:
                G_i1 = np.zeros(num_steps)
                amp_i = np.zeros(num_steps)
                tau_i = self.i_tau
                std_i = self.i_std_max
                
                if(tau_i !=0):
                    if(self.i_signalType=='Step'):
                        i0, ip1, pon1, poff1, ip2, pon2, poff2, ip3, pon3, poff3, ip4, pon4, poff4, ip5, pon5, poff5, s = self.i_heav_param
                        G_i0 = max([i0, ip1, ip2, ip3, ip4, ip5])
                        # calculate normalized value (0~1)
                        norm_i = i0 + s*((self.heav(poff1-t)*self.heav(t-pon1)*ip1)
                                + (self.heav(poff2-t)*self.heav(t-pon2)*ip2)
                                + (self.heav(poff3-t)*self.heav(t-pon3)*ip3)
                                + (self.heav(poff4-t)*self.heav(t-pon4)*ip4)
                                + (self.heav(poff5-t)*self.heav(t-pon5)*ip5))
                        norm_i /= G_i0
                    
                    elif(self.i_signalType=='Ramp'):
                        norm_i = np.zeros(num_steps)
                        pv=self.i_pValue
                        iv = self.i_iValue
                        p = self.i_period
                        ## calculate normalized value (0~1)
                        # Ramp
                        if(pv >= iv):
                            G_i0 = pv
                            if(pv == 0):
                                norm_i = np.ones(num_steps)
                            else:
                                x = iv/pv
                                norm_i =1-((1-x)/p)*np.abs(t-p) # x->1->x
                        # inverse Ramp
                        else:
                            G_i0 = iv
                            x = pv/iv
                            norm_i =x-((x-1)/p)*np.abs(t-p) # 1->x->1
                            
                    # negative values are not allowed.
                    for i in xrange(0, num_steps):
                        if(norm_i[i] < 0):
                            norm_i[i] = 0.
                    
                    exp_i = exp(-dt/tau_i)
                    amp_i = np.sqrt(norm_i)*std_i * np.sqrt( (1-np.exp(-2*dt/tau_i)) )
                    
                    # exact update rule
                    for i in xrange(1, num_steps):
                        G_i1[i] = exp_i * G_i1[i - 1] + amp_i[i - 1] * np.random.normal(loc=0.0, scale=1.0)
                    G_i = norm_i * G_i0 + G_i1
                # White noise
                else:
                    for i in xrange(0, num_steps):
                        G_i[i] = std_i * np.random.normal(loc=0.0, scale=1.0)
            
            # negative values are not allowed.
            for i in xrange(0, num_steps):
                if(G_i[i] < 0):
                    G_i[i] = 0.
            
            self.G_i = G_i
            
