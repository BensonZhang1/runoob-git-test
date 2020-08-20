# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import copy
import argparse
import json
import sys

class MeterConnector(object):
    """接线判断"""
    arrow_fmt = {'a': [dict(arrowstyle='-|>', linewidth=2, connectionstyle="arc3", color="yellow"), 'yellow'],
                                'b': [dict(arrowstyle='-|>', linewidth=2, connectionstyle="arc3", color="green"), 'green'],
                                'c': [dict(arrowstyle='-|>', linewidth=2, connectionstyle="arc3", color="red"), 'red']
                                }
    deg2rad = np.pi/180.0
    rad2deg = 180.0/np.pi
    sections = {"section11":[-30,0], "section12":[-60,-30],"section21":[-90,-60],"section22":[-120,-90],"section31":[-150,-120],"section32":[-180,-150],
                            "section41":[150,180], "section42":[120,150],"section51":[90,120],"section52":[60,90],"section61":[30,60],"section62":[0,30]}
    def __init__(self, data):
        """meter 2:
        data = {'ua':ua, 'ub':ub, 'uc':uc, 'ia':ia, 'ib':ib, 'ic':ic,
                'ctnum':ctnum, 'angleAB':angleAB, 'angleAC':angleAC, 
                'phA':phA, 'phB':phB, 'phC':phC}
      meter 2:
        data = {'uab':uab, 'ucb':ucb, 'ia':ia, 'ic':ic, 'ctnum':ctnum,
                'angleUabUcb':angleUabUcb, 'phUabIa':phUabIa, 'phUcbIc':phUcbIc}
      msg={"ua":"正常", "ub":"正常", "uc":"正常", "pf_gtet_05":{"ia":"正常", "ib":"正常", "ic":"正常"}, "pf_lt_05":{"ia":"正常", "ib":"正常", "ic":"正常"}}
    """
        self.data = data
        self.volt_fault = [0, 0, 0, 0, 0]
        self.current_small = [0, 0, 0]
        self.derection = True
        self.r_uiabc = None
        self.phi_uiabc = None
        self.msg = {"derection":"逆时针", "ua":"正常", "ub":"正常", "uc":"正常", "uab":"正常", "ucb":"正常", 
                                  "ia_is_small": "正常", "ib_is_small": "正常", "ic_is_small": "正常",
                                "pf_gtet_05":{"ia":"正常", "ib":"正常", "ic":"正常"}, 
                                "pf_lt_05":{"ia":"正常", "ib":"正常", "ic":"正常"}
                                }
        
    def check_meter_connect(self):
        msg1 = ("\n方向:{derection}\n"+
                        "ua幅值:{ua}, ub幅值:{ub}, uc幅值:{uc}\n" + 
                        "ia幅值:{ia_is_small}, ib幅值:{ib_is_small}, ic幅值:{ic_is_small}\n" ).format(**self.msg)
        msg2 = ("\n功率因数大于等于0.5时:\nia:{ia}\nib:{ib}\nic:{ic}\n" ).format(**self.msg["pf_gtet_05"])                       
        msg3 = ("\n功率因数小于0.5时:\nia:{ia}\nib:{ib}\nic:{ic}\n" ).format(**self.msg["pf_lt_05"])
        msg = msg1 + msg2 + msg3
        return msg
        
    def _is_positive_phase_sequence(self):
        if self.data["ctnum"] == 3:
            if self.data["angleAB"]>=0 and self.data["angleAC"]<0:
                self.derection = True
            elif self.data["angleAB"]<0 and self.data["angleAC"]>=0:
                self.derection = False
            elif self.data["angleAB"]>=0 and self.data["angleAC"]>=0 and self.data["angleAC"]-self.data["angleAB"]<0:
                self.derection = False
            elif self.data["angleAB"]>=0 and self.data["angleAC"]>=0 and self.data["angleAC"]-self.data["angleAB"]>=0:
                 self.derection = True
            elif self.data["angleAB"]<0 and self.data["angleAC"]<0 and self.data["angleAC"]-self.data["angleAB"]<0:
                self.derection = False
            elif self.data["angleAB"]<0 and self.data["angleAC"]<0 and self.data["angleAC"]-self.data["angleAB"]>=0:
                self.derection = True
        else:
            if self.data["angleUabUcb"]<0 :
                self.derection = True
            else:
                self.derection = False
        if self.derection:
            self.msg["derection"] = "逆时针"
        else:
            self.msg["derection"] = "顺时针"
            
    def _is_volt_fault(self):
        if self.data["ctnum"]==3:
            self.msg["uab"] = ""
            self.msg["ucb"] = ""
            if self.data["ua"] < 40 or \
                 (self.derection and (np.abs(self.data["angleAB"]-120)>20)and(np.abs(self.data["angleAC"]+120)>20))or\
                 ((not self.derection) and (np.abs(self.data["angleAB"]+120)>20)and(np.abs(self.data["angleAC"]-120)>20)):
                 self.volt_fault[0] = 1
                 self.msg["ua"] = "虚接或电压过小"
            else:
                self.volt_fault[0] = 0
                self.msg["ua"] = "正常"
                
            if self.data["ub"] < 40 or \
                 (self.derection and (np.abs(self.data["angleAB"]-120)>20))or\
                 ((not self.derection) and (np.abs(self.data["angleAB"]+120)>20)):
                 self.volt_fault[1] = 1
                 self.msg["ub"] = "虚接或电压过小"
            else:
                self.volt_fault[1] = 0
                self.msg["ub"] = "正常"
                
            if self.data["uc"] < 40 or \
                 (self.derection and (np.abs(self.data["angleAC"]+120)>20))or\
                 ((not self.derection) and (np.abs(self.data["angleAC"]-120)>20)):
                 self.volt_fault[2] = 1
                 self.msg["uc"] = "虚接或电压过小"
            else:
                self.volt_fault[2] = 0
                self.msg["uc"] = "正常"
                
        else:
            self.msg["ua"] = ""
            self.msg["ub"] = ""
            self.msg["uc"] = ""
            if self.data["uab"] <40 and \
            ((self.derection and np.abs(self.data["angleUabUcb"] + 60)>20) or 
            ((not self.derection) and np.abs(self.data["angleUabUcb"] - 60)>20)):
                self.volt_fault[3] = 1
                self.msg["uab"] = "虚接或电压过小"
            else:
                self.volt_fault[3] = 0
                self.msg["uab"] = "正常"
                
            if self.data["ucb"] <40 and \
            ((self.derection and np.abs(self.data["angleUabUcb"] + 60)>20) or
            ((not self.derection) and np.abs(self.data["angleUabUcb"] - 60)>20)):
                self.volt_fault[4] = 1
                self.msg["ucb"] = "虚接或电压过小"
            else:
                self.volt_fault[4] = 0
                self.msg["ucb"] = "正常"
                
    def _is_current_small(self):
        
        if self.data["ia"]<0.1:
            self.current_small[0] = 1
            self.msg["ia_is_small"] = "电流太小,角度判断不准"
        else:
            self.current_small[0] = 0
            self.msg["ia_is_small"] = "正常"
            
        if self.data["ctnum"]==3:      
            if self.data["ib"]<0.1:
                self.current_small[1] = 1
                self.msg["ib_is_small"] = "电流太小,角度判断不准"
            else:
                self.current_small[1] = 0
                self.msg["ib_is_small"] = "正常"
        else:
            self.msg["ib_is_small"] = ""
                
        if self.data["ic"]<0.1:
            self.current_small[2] = 1
            self.msg["ic_is_small"] = "电流太小,角度判断不准"
        else:
            self.current_small[2] = 0 
            self.msg["ic_is_small"] = "正常"
                
    @staticmethod
    def _pu_UI(uabc, iabc, ctnum):
        """ uabc:[ua, ub, uc] or [uab, ucb]
       iabc:[ia, ib, ic] or [ia, ic] """
        rst = {}
        max_U = np.max(uabc)
        max_I = np.max(iabc)
        if ctnum == 3:
            rst['ua'] = uabc[0]/(max_U+0.01)
            rst['ub'] = uabc[1]/(max_U+0.01)
            rst['uc'] = uabc[2]/(max_U+0.01)
            
            rst['ia'] = iabc[0]/(max_I+0.01)
            rst['ib'] = iabc[1]/(max_I+0.01)
            rst['ic'] = iabc[2]/(max_I+0.01)
        else:
            rst['uab'] = uabc[0]/(max_U+0.01)
            rst['ucb'] = uabc[1]/(max_U+0.01)
            
            rst['ia'] = iabc[0]/(max_I+0.01)
            rst['ic'] = iabc[1]/(max_I+0.01)     
        return rst
        
    def _phi_cal(self, phi_uabc, phi_iabc, ctnum):
        """ phi_uabc:[angleAB, angleAC] or [angleUabUcb]
       phi_iabc:[phA, phB, phC] or [phUabIa, phUcbIc] """    
        phi = {}
        if ctnum == 3:
            phi['ua'] = 0
            phi['ub'] = -phi_uabc[0]*self.deg2rad
            phi['uc'] = -phi_uabc[1]*self.deg2rad
            
            phi['ia'] = -phi_iabc[0]*self.deg2rad
            phi['ib'] = phi['ub'] - phi_iabc[1]*self.deg2rad
            phi['ic'] = phi['uc'] - phi_iabc[2]*self.deg2rad
        else:
            if phi_uabc[0] > 0:
                ph_uab_init = -30*self.deg2rad
            else:
                ph_uab_init = 30*self.deg2rad
            phi['uab'] = ph_uab_init
            phi['ucb'] = phi['uab'] - phi_uabc[0]*self.deg2rad
            
            phi['ia'] = phi['uab']  - phi_iabc[0]*self.deg2rad
            phi['ic'] = phi['ucb'] - phi_iabc[1]*self.deg2rad
        return phi
        
    @staticmethod    
    def _deg_180to180(phi):
        phi %= 360
        if phi>180:
            phi -= 360
        elif phi<= -180:
            phi += 360
        return phi
    
    def _find_section(self,phi):
        phi = self._deg_180to180(phi*self.rad2deg)
        
        if phi>=self.sections["section11"][0] and phi<self.sections["section11"][1]:
            section_name = ("section11", "section1")
        elif phi>=self.sections["section12"][0] and phi<self.sections["section12"][1]:
            section_name = ("section12", "section1")
        elif phi>=self.sections["section21"][0] and phi<self.sections["section21"][1]:
            section_name = ("section21", "section2")
        elif phi>=self.sections["section22"][0] and phi<self.sections["section22"][1]:
            section_name = ("section22", "section2")            
        elif phi>=self.sections["section31"][0] and phi<self.sections["section31"][1]:
            section_name = ("section31", "section3")
        elif phi>=self.sections["section32"][0] and phi<self.sections["section32"][1]:
            section_name = ("section32", "section3")
        elif phi>=self.sections["section41"][0] and phi<self.sections["section41"][1]:
            section_name = ("section41", "section4")
        elif phi>=self.sections["section42"][0] and phi<self.sections["section42"][1]:
            section_name = ("section42", "section4")
        elif phi>=self.sections["section51"][0] and phi<self.sections["section51"][1]:
            section_name = ("section51", "section5")
        elif phi>=self.sections["section52"][0] and phi<self.sections["section52"][1]:
            section_name = ("section52", "section5")            
        elif phi>=self.sections["section61"][0] and phi<self.sections["section61"][1]:
            section_name = ("section61", "section6")   
        elif phi>=self.sections["section62"][0] and phi<self.sections["section62"][1]:
           section_name = ("section62", "section6")  
        return section_name

    def _positive_sequence_ct(self):
        #功率因数大于等于0.5,
        pos_seq_ct_table1 = {
                                                    "ia":{"section1":"正常",                                                                                   "section2":"采集设备A相CT>---反向---<<---<工厂C相",
                                                                   "section3":"采集设备A相CT>---正向--->>---<工厂B相", "section4":"采集设备A相CT>---反向---<<---<工厂A相",
                                                                   "section5":"采集设备A相CT>---正向--->>---<工厂C相", "section6":"采集设备A相CT>---反向---<<---<工厂B相"}, 
                                                                
                                                    "ib":{"section1":"采集设备B相CT>---正向--->>---<工厂A相", "section2":"采集设备B相CT>---反向---<<---<工厂C相",
                                                                   "section3":"正常",                                                                                    "section4":"采集设备B相CT>---反向---<<---<工厂A相",
                                                                   "section5":"采集设备B相CT>---正向--->>---<工厂C相", "section6":"采集设备B相CT>---反向---<<---<工厂B相"},
                                                                
                                                    "ic":{"section1":"采集设备C相CT>---正向--->>---<工厂A相", "section2":"采集设备C相CT>---反向---<<---<工厂C相",
                                                                   "section3":"采集设备C相CT>---正向--->>---<工厂B相", "section4":"采集设备C相CT>---反向---<<---<工厂A相",
                                                                   "section5":"正常",                                                                                    "section6":"采集设备C相CT>---反向---<<---<工厂B相"}
                                                    }
                                                    
        ia_section = self._find_section(self.phi_uiabc["ia"])
        if self.data["ctnum"] == 3:
            ib_section = self._find_section(self.phi_uiabc["ib"])
        ic_section = self._find_section(self.phi_uiabc["ic"])

        self.msg["pf_gtet_05"]["ia"] = pos_seq_ct_table1["ia"][ia_section[1]]
        if self.data["ctnum"] == 3:
            self.msg["pf_gtet_05"]["ib"] = pos_seq_ct_table1["ib"][ib_section[1]]
        else:
            self. msg["pf_gtet_05"]["ib"] = ""
        self.msg["pf_gtet_05"]["ic"] = pos_seq_ct_table1["ic"][ic_section[1]]
        #功率因数小于0.5，
        pos_seq_ct_table2 = {
                                                    "ia":{"section11":"采集设备A相CT>---反向---<<---<工厂B相", "section21":"正常",
                                                                   "section31":"采集设备A相CT>---反向---<<---<工厂C相", "section41":"采集设备A相CT>---正向--->>---<工厂B相",
                                                                   "section51":"采集设备A相CT>---反向---<<---<工厂A相", "section61":"采集设备A相CT>---正向--->>---<工厂C相"}, 
                                                                
                                                    "ib":{"section11":"采集设备B相CT>---反向---<<---<工厂B相", "section21":"采集设备B相CT>---正向--->>---<工厂A相",
                                                                   "section31":"采集设备B相CT>---反向---<<---<工厂C相", "section41":"正常",
                                                                   "section51":"采集设备B相CT>---反向---<<---<工厂A相", "section61":"采集设备B相CT>---正向--->>---<工厂C相"},
                                                                
                                                    "ic":{"section11":"采集设备C相CT>---反向---<<---<工厂B相", "section21":"采集设备C相CT>---正向--->>---<工厂A相",
                                                                   "section31":"采集设备C相CT>---反向---<<---<工厂C相", "section41":"采集设备C相CT>---正向--->>---<工厂B相",
                                                                   "section51":"采集设备C相CT>---反向---<<---<工厂A相", "section61":"正常"}
                                                    }       

        if ia_section[0] in ["section12", "section22", "section32", "section42", "section52", "section62"]:
           self. msg["pf_lt_05"]["ia"] = "ia功率因数不可能小于0.5"
        else:
           self. msg["pf_lt_05"]["ia"] = pos_seq_ct_table2["ia"][ia_section[0]]
            
        if self.data["ctnum"] == 3:
            if ib_section[0] in ["section12", "section22", "section32", "section42", "section52", "section62"]:
               self. msg["pf_lt_05"]["ib"] = "ib功率因数不可能小于0.5"
            else:
               self. msg["pf_lt_05"]["ib"] = pos_seq_ct_table2["ib"][ib_section[0]]
        else:
            self. msg["pf_lt_05"]["ib"] = ""
                
        if ic_section[0] in ["section12", "section22", "section32", "section42", "section52", "section62"]:
           self. msg["pf_lt_05"]["ic"] = "ic功率因数不可能小于0.5"
        else:
            self.msg["pf_lt_05"]["ic"] = pos_seq_ct_table2["ic"][ic_section[0]]
            
    def _negative_sequence_ct(self):
        #功率因数大于等于0.5,
       neg_seq_ct_table1 = {
                                                    "ia":{"section1":"正常",                                                                                   "section2":"采集设备A相CT>---反向---<<---<工厂B相",
                                                                   "section3":"采集设备A相CT>---正向--->>---<工厂C相", "section4":"采集设备A相CT>---反向---<<---<工厂A相",
                                                                   "section5":"采集设备A相CT>---正向--->>---<工厂B相", "section6":"采集设备A相CT>---反向---<<---<工厂C相"}, 
                                                                
                                                    "ib":{"section1":"采集设备B相CT>---正向--->>---<工厂A相", "section2":"采集设备B相CT>---反向---<<---<工厂B相",
                                                                   "section3":"采集设备B相CT>---正向--->>---<工厂C相", "section4":"采集设备B相CT>---反向---<<---<工厂A相",
                                                                   "section5":"正常",                                                                                    "section6":"采集设备B相CT>---反向---<<---<工厂C相"},
                                                                
                                                    "ic":{"section1":"采集设备C相CT>---正向--->>---<工厂A相", "section2":"采集设备C相CT>---反向---<<---<工厂B相",
                                                                   "section3":"正常",                                                                                   "section4":"采集设备C相CT>---反向---<<---<工厂A相",
                                                                   "section5":"采集设备C相CT>---正向--->>---<工厂B相", "section6":"采集设备C相CT>---反向---<<---<工厂C相"}
                                                    }
                                                    
       ia_section = self._find_section(self.phi_uiabc["ia"])
       if self.data["ctnum"] == 3:
          ib_section = self._find_section(self.phi_uiabc["ib"])
       ic_section = self._find_section(self.phi_uiabc["ic"])

       self.msg["pf_gtet_05"]["ia"] = neg_seq_ct_table1["ia"][ia_section[1]]
       if self.data["ctnum"] == 3:
           self.msg["pf_gtet_05"]["ib"] = neg_seq_ct_table1["ib"][ib_section[1]]
       else:
            self. msg["pf_gtet_05"]["ib"] = ""
       self.msg["pf_gtet_05"]["ic"] = neg_seq_ct_table1["ic"][ic_section[1]]
        #功率因数小于0.5，
       neg_seq_ct_table2 = {
                                                    "ia":{"section11":"采集设备A相CT>---反向---<<---<工厂C相", "section21":"正常",
                                                                   "section31":"采集设备A相CT>---反向---<<---<工厂B相", "section41":"采集设备A相CT>---正向--->>---<工厂C相",
                                                                   "section51":"采集设备A相CT>---反向---<<---<工厂A相", "section61":"采集设备A相CT>---正向--->>---<工厂B相"}, 
                                                                
                                                    "ib":{"section11":"采集设备B相CT>---反向---<<---<工厂C相", "section21":"采集设备B相CT>---正向--->>---<工厂A相",
                                                                   "section31":"采集设备B相CT>---反向---<<---<工厂B相", "section41":"采集设备B相CT>---正向--->>---<工厂C相",
                                                                   "section51":"采集设备B相CT>---反向---<<---<工厂A相", "section61":"正常"},
                                                                
                                                    "ic":{"section11":"采集设备C相CT>---反向---<<---<工厂C相", "section21":"采集设备C相CT>---正向--->>---<工厂A相",
                                                                   "section31":"采集设备C相CT>---反向---<<---<工厂B相", "section41":"正常",
                                                                   "section51":"采集设备C相CT>---反向---<<---<工厂A相", "section61":"采集设备C相CT>---正向--->>---<工厂B相"}
                                                    }       

       if ia_section[0] in ["section12", "section22", "section32", "section42", "section52", "section62"]:
           self.msg["pf_lt_05"]["ia"] = "ia功率因数不可能小于0.5"
       else:
           self.msg["pf_lt_05"]["ia"] = neg_seq_ct_table2["ia"][ia_section[0]]
            
       if self.data["ctnum"] == 3:
           if ib_section[0] in ["section12", "section22", "section32", "section42", "section52", "section62"]:
              self. msg["pf_lt_05"]["ib"] = "ib功率因数不可能小于0.5"
           else:
               self.msg["pf_lt_05"]["ib"] = neg_seq_ct_table2["ib"][ib_section[0]]
       else:
            self. msg["pf_lt_05"]["ib"] = ""
            
       if ic_section[0] in ["section12", "section22", "section32", "section42", "section52", "section62"]:
           self.msg["pf_lt_05"]["ic"] = "ic功率因数不可能小于0.5"
       else:
           self.msg["pf_lt_05"]["ic"] = neg_seq_ct_table2["ic"][ic_section[0]]
     
    def output(self):
        self. _is_positive_phase_sequence()
        self._is_volt_fault()
        self._is_current_small()
        if self.data["ctnum"] == 3:
            self.r_uiabc = self._pu_UI([self.data['ua'], self.data['ub'], self.data['uc']], 
                                                                            [self.data['ia'], self.data['ib'], self.data['ic']], 
                                                                            self.data['ctnum'])       
            self.phi_uiabc = self._phi_cal([self.data['angleAB'], self.data['angleAC']],
                                                                    [self.data['phA'], self.data['phB'], self.data['phC']], 
                                                                    self.data['ctnum'])
        else:
            self.r_uiabc = self._pu_UI([self.data['uab'], self.data['ucb']], 
                                                                            [self.data['ia'], self.data['ic']], 
                                                                            self.data['ctnum'])
            self.phi_uiabc = self._phi_cal([self.data['angleUabUcb']],
                                                                [self.data['phUabIa'], self.data['phUcbIc']], 
                                                                self.data['ctnum'])
        
        if self.derection :
            self._positive_sequence_ct()
        else:
            self._negative_sequence_ct()
            
    @staticmethod
    def _plot_setting(style='dark_background', facecolor='k'):
        plt.style.use(style)
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='polar', facecolor=facecolor)
        plt.grid(color='w', linestyle='--', linewidth=2, alpha=0.4)
        ax.tick_params(axis='both',width=2,colors='w')
        ax.set_rgrids(np.arange(0.3, 1.51, 0.4))
        ax.set_thetagrids(np.arange(0.0, 360.0, 30), [theta for theta in range(0, 365, 30)])
        return ax
     
    @staticmethod
    def _polar_var(phi, r, name, fmt):
        plt.annotate('', xy=(phi, r), xytext=(0,0), arrowprops = fmt[0])
        plt.text(phi + 0.0, r + 0.05, name, color=fmt[1], fontsize=14)
        
    def plot(self,style='dark_background', facecolor='k'):
        arrow_rad = 0.3
        ax = self._plot_setting()
      
        if self.data['ctnum'] == 3:
            if not self.derection:
                ax.set_theta_direction(-1)
                arrow_rad = -arrow_rad
            for item in zip(['ua', 'ub', 'uc'], [self.arrow_fmt['a'], self.arrow_fmt['b'], self.arrow_fmt['c']]):
                self._polar_var(self.phi_uiabc[item[0]], self.r_uiabc[item[0]], item[0], item[1])
                
            for item in zip(['ia', 'ib', 'ic'], [self.arrow_fmt['a'], self.arrow_fmt['b'], self.arrow_fmt['c']]):
                self._polar_var(self.phi_uiabc[item[0]], self.r_uiabc[item[0]], item[0], item[1])

        else:                                                              
            if not self.derection:
                ax.set_theta_direction(-1)
                delt_phi = -30*self.deg2rad
                arrow_rad = -arrow_rad
            else:
                delt_phi = 30*self.deg2rad
 
            for item in zip(['uab', 'ucb'], [self.arrow_fmt['a'], self.arrow_fmt['c']]):
                self._polar_var(self.phi_uiabc[item[0]], self.r_uiabc[item[0]], item[0], item[1]) 

            for item in zip(['ia', 'ic'], [self.arrow_fmt['a'], self.arrow_fmt['c']]):
                self._polar_var(self.phi_uiabc[item[0]], self.r_uiabc[item[0]], item[0], item[1]) 
            
            self._polar_var(self.phi_uiabc['uab']-delt_phi, self.r_uiabc['uab']/1.732, 'ua', self.arrow_fmt['a'])
            self._polar_var(-self.phi_uiabc['uab']-self.phi_uiabc['ucb'], self.r_uiabc['uab']/1.732, 'ub', self.arrow_fmt['b'])
            self._polar_var(self.phi_uiabc['ucb']+delt_phi, self.r_uiabc['ucb']/1.732, 'uc', self.arrow_fmt['c'])

        plt.annotate('', xy=(30*self.deg2rad, 1.4), xytext=(0,1.4),
                                      arrowprops = dict(arrowstyle='-|>', linewidth=2, 
                                      connectionstyle="arc3,rad=%f"%arrow_rad, color="red"))
        plt.show()
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "json_file_path",
        metavar='json_file_path',
        type=str,
        help=("pull the path of json file")
    )
    data = {'ctnum':2, 'angleUabUcb':-60, 'ia':0.896, 'ic':0.896, 'phUabIa':50.799, 'phUcbIc':-9.4, 'uab':103.79, 'ucb':103.42}
    #data = {'ctnum':2, 'angleUabUcb':60, 'ia':3.009, 'ic':3.071, 'phUabIa':4.7, 'phUcbIc':62.5, 'uab':103.79, 'ucb':103.42}
    # data = {'ctnum':3, 'ua':60, 'ub':60,'uc':60,'ia':1.707, 'ib':1.35, 'ic':1.53, 'angleAB':120, 'angleAC':-120, 'phA':33.09, 'phB':28.9,'phC':37}
    #data = {'ctnum':3, 'ua':60, 'ub':60,'uc':60,'ia':2.27, 'ib':2.23, 'ic':2.35, 'angleAB':-120, 'angleAC':120, 'phA':11.9, 'phB':10.08,'phC':12.1}
    #file = r"E:\\工作资料\\查看采样是否正确\\15-38_1586931398.json"
    file = sys.argv[1]
    f =  open(file, 'r')
    json_data = json.loads(f.read())
    data1 = json_data["data"]
    data_var2 = ['ctnum', 'angleUabUcb', 'ia', 'ic', 'phUabIa', 'phUcbIc', 'uab', 'ucb']
    data_var3 = ['ctnum', 'ua', 'ub','uc','ia', 'ib', 'ic', 'angleAB', 'angleAC', 'phA', 'phB','phC']
    data = {}
    if data1["ctnum"] == 3:
        for key in data_var3:
            data[key] = data1[key]
    else:
        for key in data_var2:
            data[key] = data1[key]
    data = MeterConnector(data)
    data.output()
    print(data.check_meter_connect())
    data.plot()