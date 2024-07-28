import tkinter as tk
from tkinter.font import Font
import csv
from tkinter import Canvas
import socket, threading
import time
import numpy as np
import datetime
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import urllib.request
import os
from pathlib import WindowsPath, Path
from datetime import datetime
from decimal import Decimal
import ROOT

class xpm_analysis(tk.Frame):
    def elifetimetrend(self) :
        #plt.clf()
        ROOT.gErrorIgnoreLevel = 6001
        self.tree.Reset()
        try :
            self.tree.ReadFile(path+file_name+'.txt','Tc:Ta:TcRise:TaRise:cat:an:offst:datime:IR:UV:chi2:nAvg:nTrig:cat_ll:cat_ul:an_ll:an_ul',',')
            self.tree.Draw('Entry$:datime','','goff')
            entry0 = self.tree.GetV1()[0]
        except ReferenceError :
            try :
                self.tree = ROOT.TTree('xpmdata','')
                self.tree.ReadFile(path+file_name+'.txt','Tc:Ta:TcRise:TaRise:cat:an:offst:datime:IR:UV:chi2:nAvg:nTrig',',')
                self.tree.Draw('Entry$:datime','','goff')
                entry0 = self.tree.GetV1()[0]
            except ReferenceError :
                self.tree = ROOT.TTree('xpmdata','')
                self.tree.ReadFile(path+file_name+'.txt','Tc:Ta:TcRise:TaRise:cat:an:offst:datime:IR:UV',',')

        self.tree.Draw('(Tc-Ta)/log(an/cat):datime','','goff')
        lf_time = vtoa( self.tree.GetV1() , self.tree.GetSelectedRows() - 1 )
        t0 = self.tree.GetV2()[0]
        t1 = self.tree.GetV2()[self.tree.GetEntries()-1]
        tau = vtoa( self.tree.GetV1() , self.tree.GetEntries() )
        maxtau = np.max( tau )
        bin_by_fibersave = False
        avg_samples = int(self.obsPerBin.get())
        binning_mode = (self.binbyfibersave.get())
        if binning_mode == 1 :
            ##################### binning by Fiber Save group ################## 
            bin_by_fibersave = True
            timebucket = []
            tlow = []
            tcent = []
            for entry in range(0,self.tree.GetSelectedRows()-1) :
                ind = entry - 1
                if ind < 0 :
                    ind = 0
                dtl = self.tree.GetV2()[entry] - self.tree.GetV2()[ind]
                dtr = self.tree.GetV2()[entry+1] - self.tree.GetV2()[entry]
                if  dtl < 1800.0 :
                    timebucket.append(self.tree.GetV2()[entry])
                elif len(timebucket) < avg_samples :
                    timebucket.append(self.tree.GetV2()[entry])
                else :
                    try :
                        tlow.append( np.array(timebucket).min()-1.0 )
                        tlow.append( np.array(timebucket).max()+1.0 )
                        tcent.append( np.array(timebucket).mean() )
                        #print(len(timebucket),(self.tree.GetV2()[entry-1]-t0)/3600.0)
                        timebucket = []
                        timebucket.append(self.tree.GetV2()[entry])
                    except ValueError :
                        pass
            tlow.append( np.array(timebucket).min() )
            tlow.append( np.array(timebucket).max() )
            tcent.append( np.array(timebucket).mean() )
            tlow = (np.array(tlow) - t0)/3600.0
            tcent = (np.array(tcent) - t0)/3600.0
            print('number of bins',len(tlow))
            ####################################################################
            self.myhist = ROOT.TH2F('myhist','',len(tlow)-1,tlow,200,0.0,maxtau)            
        else :
            nbinsX = int( self.tree.GetEntries()/avg_samples )
            self.myhist=ROOT.TH2F('myhist','',nbinsX,0.0,(t1-t0)/3600.0,int(maxtau/100.0),0.0,maxtau)
         
        self.tree.Draw('(Tc-Ta)/log(an/cat):(datime-'+str(t0)+')/3600.0','','goff')
        norm_time_clean = vtoa( self.tree.GetV2() , self.tree.GetSelectedRows() - 1 )
        com = '(Tc-Ta)/log(an/cat):(datime-'+str(t0)+')/3600.0>>myhist'
        mindiff = float(self.minacdiff.get())
        tcut_base = 'UV>30 && cat>0 && an>0 && an<(cat-'+str(mindiff)+')'
        tcut = tcut_base
        myhist_el = self.myhist
        myhist_eh = self.myhist
        print(self.tree.Draw(com,tcut,'goff'))
        myprof = self.myhist.ProfileX()
        myprof_noweights = myprof
        if self.tree.GetNbranches() > 13 : #point-by-point statistical errors have been computed 
            tcut = '(' + tcut_base + ')*TMath::Power(1.0/((Ta-Tc)/log(cat_ll/an_ul) - (Ta-Tc)/log(cat_ul/an_ll)),2)' 
            print(self.tree.Draw(com,tcut,'goff'))
            myhist_el=ROOT.TH2F('myhist_el','',nbinsX,0.0,(t1-t0)/3600.0,int(100.0),0.0,30.0)
            myhist_eh=ROOT.TH2F('myhist_eh','',nbinsX,0.0,(t1-t0)/3600.0,int(100.0),0.0,30.0)
            print(self.tree.Draw(com,tcut_base,'goff'))
            myprof_noweights = self.myhist.ProfileX()

        x = []
        y = []
        e = []
        ex = []
        for bin in range(1,myprof.GetNbinsX()+1) :
            if( myprof_noweights.GetBinEntries(bin)<=1 ): continue
            x.append(myprof.GetBinCenter(bin))
            y.append(myprof.GetBinContent(bin))
            e.append(myprof.GetBinError(bin))
            ex.append(myprof.GetBinWidth(bin)/2.0)
            if self.tree.GetNbranches() > 13 : #point-by-point statistical errors have been computed 
                eh_cut = ROOT.TCut('(fabs( (datime-'+str(t0)+')/3600.0 - ' + str(x[-1]) + ') <= ' + str(myprof.GetBinWidth(bin)/2.0) + ')')
                eh_cut = eh_cut and ROOT.TCut('(Ta-Tc)/log(cat/an) < ' + str(y[-1]))
                eh_cut = eh_cut and ROOT.TCut( tcut_base )
                eh_cut = eh_cut * ROOT.TCut('TMath::Power(1.0/((Ta-Tc)/log(cat_ll/an_ul) - (Ta-Tc)/log(cat/an)),2)') 
                com_eh = str(y[-1]) + '- (Ta-Tc)/log(cat/an)' 
                com_eh = com_eh + ':(datime-'+str(t0)+')/3600.0>>myhist_eh'
                self.tree.Draw(com_eh,eh_cut,'goff')
                prof_eh = myhist_eh.ProfileX()
                variance = myhist_eh.ProjectionY().GetRMS()**2
                mean = prof_eh.GetBinContent(bin)
                ehy.append( np.sqrt( mean**2 + variance ) )
    
                el_cut = ROOT.TCut('(fabs( (datime-'+str(t0)+')/3600.0 - ' + str(x[-1]) + ') <= ' + str(myprof.GetBinWidth(bin)/2.0) + ')')
                el_cut = el_cut and ROOT.TCut('(Ta-Tc)/log(cat/an) > ' + str(y[-1]))
                el_cut = el_cut and ROOT.TCut( tcut_base )
                el_cut = el_cut * ROOT.TCut('TMath::Power(1.0/((Ta-Tc)/log(cat/an) - (Ta-Tc)/log(cat_ul/an_ll)),2)') 
                print(el_cut.GetTitle())
                com_el = '(Ta-Tc)/log(cat/an) - ' + str(y[-1]) 
                com_el = com_el + ':(datime-'+str(t0)+')/3600.0>>myhist_el'
                self.tree.Draw(com_el,el_cut,'goff')
                prof_el = myhist_el.ProfileX()
                variance = myhist_el.ProjectionY().GetRMS()**2
                mean = prof_el.GetBinContent(bin)
                ely.append( np.sqrt( mean**2 + variance ) )

        if self.isrational == 0 :
            fitfunc = 'E'
        else :
            fitfunc = 'R'
        
        fexp = ROOT.TF1('fexp','expo + [2]',x[0],x[-1])
        row_labels=['$\chi^2$/ndf','Constant','Slope','Baseline']
        if fitfunc=='R' :
            fexp = ROOT.TF1('fexp','[0]/([1]+[2]*x)',x[0],x[-1])
            fexp.SetParameters(myprof.GetBinContent(1)*2,2.0,0.005) 
            row_labels=['$\chi^2$/ndf','A','B','C']
        try :
          mi,ma = str(self.fitdomain).split()
          fexp.SetRange(float(mi),float(ma))
          myprof.Fit(fexp,'NRE')
        except ValueError:
          myprof.Fit(fexp,'N')
          fexp.SetRange(x[0],x[-1]) 
          mi = x[0]
          ma = x[-1]
        #myprof.Fit(fexp,'NM')
        fpts = 500
        yfit = []
        xfit = []
        for pt in range(0,fpts):
            x_pt =  x[0] + pt*(x[-1] - x[0])/fpts
            if x_pt >= float(mi) and x_pt <= float(ma) :
              xfit.append( x_pt )
              yfit.append( fexp.Eval( x_pt ) )
        col_labels=['Fit']
        table_vals=[['{0:.3g}'.format(fexp.GetChisquare())+'/'+'{0:.3g}'.format(fexp.GetNDF())],
                    ['{0:.3g}'.format(fexp.GetParameter(0))+'$\pm$'+'{0:.3g}'.format(fexp.GetParError(0))],
                    ['{0:.3g}'.format(fexp.GetParameter(1))+'$\pm$'+'{0:.3g}'.format(fexp.GetParError(1))],
                    ['{0:.3g}'.format(fexp.GetParameter(2))+'$\pm$'+'{0:.3g}'.format(fexp.GetParError(2))]]
        # the rectangle is where I want to place the table
        #plt.show()
        ########
        subsample_norm_time = norm_time_clean[::avg_samples]
        #lf_time = (raw_input_file_clean[:,1]-raw_input_file_clean[:,0])/np.log(raw_input_file_clean[:,4]/raw_input_file_clean[:,5])
        avg,err = smooth(lf_time,avg_samples)
        llsel = self.plotting_option.get()[0] 
        if llsel == str(1): #1 - Average only
            plt.close()
            the_table = plt.table(cellText=table_vals,
                                  rowLabels=row_labels,
                                  colLabels=col_labels,
                                  fontsize=8,
                                  colWidths=[0.3]*3,
                                  loc='upper right')
            #plt.title(tcut+' ; '+ str(avg_samples)+' samples/bin')
            plt.title(tcut+'\n'+ str(np.around(myprof.GetEntries()/myprof.GetNbinsX(),1))+' post-cut samples/bin')
            
            if self.tree.GetNbranches() > 13 : #point-by-point statistical errors have been computed 
                errormatrix=np.array([np.array(ely),np.array(ehy)])
                plt.errorbar(x,y,xerr=ex,yerr=errormatrix,fmt='r.')
            else :
                plt.errorbar(x,y,e,ex,fmt='ro')
            
            plt.plot(xfit,yfit,'b-')
            plt.xlabel('Hours')
            plt.annotate('Average', xy=(0.1, 0.95), xycoords='axes fraction',color='red',weight='bold')
            plt.ylabel('e$^{-}$ lifetime [$\mu$s]')
            #plt.show(block=False)
        elif llsel == str(2): #2 - Scatter points only
            self.tree.Draw(com,tcut,'goff')
            cut_time = vtoa( self.tree.GetV2() , self.tree.GetSelectedRows() - 1 )
            cut_tau = vtoa( self.tree.GetV1() , self.tree.GetSelectedRows() - 1 )
            plt.close()
            plt.plot(norm_time_clean,lf_time,'go',markersize=1.5,zorder=-32)           
            plt.plot(cut_time,cut_tau,'m^',markersize=1.5,zorder=-32)           
            plt.xlabel('Hours')
            plt.ylabel('e$^{-}$ lifetime [$\mu$s]')
            plt.annotate('Scatter plot only', xy=(0.1, 0.95), xycoords='axes fraction',color='red',weight='bold')
            #plt.show(block=False)
        elif llsel == str(3):#3 - Scatter+Average
            plt.close()
            the_table = plt.table(cellText=table_vals,
                                  rowLabels=row_labels,
                                  colLabels=col_labels,
                                  fontsize=8,
                                  colWidths=[0.3]*3,
                                  loc='upper right')
            #plt.title(tcut+' ; '+ str(avg_samples)+' samples/bin')
            plt.title(tcut+'\n'+ str(np.around(myprof.GetEntries()/myprof.GetNbinsX(),1))+' post-cut samples/bin')
            plt.errorbar(x,y,e,ex,fmt='ro')
            plt.plot(xfit,yfit,'b-')
            self.tree.Draw(com,tcut,'goff')
            cut_time = vtoa( self.tree.GetV2() , self.tree.GetSelectedRows() - 1 )
            cut_tau = vtoa( self.tree.GetV1() , self.tree.GetSelectedRows() - 1 )
            plt.plot(norm_time_clean,lf_time,'go',markersize=1.5,zorder=-32)           
            plt.plot(cut_time,cut_tau,'m^',markersize=1.5,zorder=-32)           
            #plt.errorbar(subsample_norm_time,avg,yerr=err/np.sqrt(avg_samples-1),fmt='ro')
            #plt.title('Raw Lifetime and Lifetime Averaged Over %i Samples'%avg_samples)
            if time_choice == 'day': plt.xticks(rotation=25)
            if time_choice == 'hour': plt.xlabel('Hours') 
            plt.ylabel('e$^{-}$ lifetime [$\mu$s]')
            plt.annotate('Average', xy=(0.1, 0.95), xycoords='axes fraction',color='red',weight='bold')
            plt.show(block=False)
        elif llsel == str(4):#4 - Median (lognormal)
            plt.close()
            self.tree.Draw('log((Tc-Ta)/log(an/cat)):datime','','goff')
            t0 = self.tree.GetV2()[0]
            t1 = self.tree.GetV2()[self.tree.GetEntries()-1]
            lnhist=ROOT.TH2F()
            tau = vtoa( self.tree.GetV1() , self.tree.GetEntries() )
            maxtau = np.max( tau )
            #maxtau = 50000.0
            lnhist=ROOT.TH2F('lnhist','',nbinsX,0.0,(t1-t0)/3600.0,int(np.exp(maxtau)/100.0),0.0,maxtau)
            com = 'log((Tc-Ta)/log(an/cat)):(datime-'+str(t0)+')/3600.0>>lnhist'
            tcut = 'UV>30 && cat>0 && an>0 && an<(cat-'+str(mindiff)+')'
            print(self.tree.Draw(com,tcut,'goff'))
            lnprof = lnhist.ProfileX()
            lnprof.SetMarkerStyle(20)
            x = []
            y = []
            ely = []
            ehy = []
            ex = []
            pxx = lnhist.ProjectionX('pxx',1,2)
            lngraph=ROOT.TGraphAsymmErrors()
            for bin in range(1, lnprof.GetNbinsX()+1 ) :
              if( lnprof.GetBinEntries(bin)<=1 ): continue
              x.append(lnprof.GetBinCenter(bin))
              y.append(np.exp(lnprof.GetBinContent(bin)))
              ely.append(np.exp(lnprof.GetBinContent(bin)) - np.exp(lnprof.GetBinContent(bin) - lnprof.GetBinError(bin)))
              ehy.append(np.exp(lnprof.GetBinContent(bin) + lnprof.GetBinError(bin)) - np.exp(lnprof.GetBinContent(bin)) )
              ex.append(lnprof.GetBinWidth(bin)/2.0)
              lngraph.SetPoint(lngraph.GetN(),x[-1],y[-1])
              lngraph.SetPointError(lngraph.GetN()-1,0.0,0.0,ely[-1],ehy[-1])
              #print(ely[bin-1],ehy[bin-1])
            #fexp = ROOT.TF1('fexp','expo + [2]',x[0],x[-1])
            lngraph.Fit(fexp,'NO')
            fpts = 500
            yfit = []
            xfit = []
            for pt in range(0,fpts):
              xfit.append( x[0] + pt*(x[-1] - x[0])/fpts )
              yfit.append( fexp.Eval( xfit[pt] ) )
            col_labels=['Fit']
            #row_labels=['$\chi^2$/ndf','Constant','Slope','Baseline']
            table_vals=[['{0:.3g}'.format(fexp.GetChisquare())+'/'+'{0:.3g}'.format(fexp.GetNDF())],
                        ['{0:.3g}'.format(fexp.GetParameter(0))+'$\pm$'+'{0:.3g}'.format(fexp.GetParError(0))],
                        ['{0:.3g}'.format(fexp.GetParameter(1))+'$\pm$'+'{0:.3g}'.format(fexp.GetParError(1))],
                        ['{0:.3g}'.format(fexp.GetParameter(2))+'$\pm$'+'{0:.3g}'.format(fexp.GetParError(2))]]
            the_table = plt.table(cellText=table_vals,
                                  rowLabels=row_labels,
                                  colLabels=col_labels,
                                  fontsize=8,
                                  colWidths=[0.3]*3,
                                  loc='upper right')
            #plt.title(tcut+' ; '+ str(avg_samples)+' samples/bin')
            plt.title(tcut+'\n'+ str(np.around(myprof.GetEntries()/myprof.GetNbinsX(),1))+' post-cut samples/bin')
            errormatrix=np.array([np.array(ely),np.array(ehy)])
            plt.errorbar(x,y,xerr=ex,yerr=errormatrix,markersize=4,fmt='ro')
            plt.plot(xfit,yfit,'b-')
            com = '((Tc-Ta)/log(an/cat)):(datime-'+str(t0)+')/3600.0'
            self.tree.Draw(com,tcut,'goff')
            cut_time = vtoa( self.tree.GetV2() , self.tree.GetSelectedRows() - 1 )
            cut_tau = vtoa( self.tree.GetV1() , self.tree.GetSelectedRows() - 1 )
            ptsel = self.superimpose.get() 
            if ptsel == 1 :
                plt.plot(norm_time_clean,lf_time,'go',markersize=1.5,zorder=-32)           
                plt.plot(cut_time,cut_tau,'m^',markersize=1.5,zorder=-32)           
            if time_choice == 'day': plt.xticks(rotation=25)
            if time_choice == 'hour': plt.xlabel('Hours') 
            plt.annotate('Median (assumes lognormality)', xy=(0.1, 0.95), xycoords='axes fraction',color='red',weight='bold')
            plt.ylabel('e$^{-}$ lifetime [$\mu$s]')
        elif llsel == str(5):#5 - Mode
            plt.close()
            self.tree.Draw('log((Tc-Ta)/log(an/cat)):datime','','goff')
            self.myhist=ROOT.TH2F()
            self.myhist=ROOT.TH2F('myhist','',nbinsX,0.0,(t1-t0)/3600.0,int(75*maxtau/100000.0),0.0,maxtau)
            print(self.myhist.ProjectionY('py_t',1,1).GetBinWidth(1))
            self.tree.Draw(com,tcut,'goff')
            x = []
            y = []
            ely = []
            ehy = []
            ex = []
            pxx = self.myhist.ProjectionX('pxx',1,2)
            modgraph=ROOT.TGraphAsymmErrors()
            for bin in range(1, self.myhist.GetNbinsX()+1 ) :
              if( myprof.GetBinEntries(bin)<=1 ): continue
              pyy = self.myhist.ProjectionY('pyy',bin,bin)
              mod = pyy.GetBinCenter(pyy.GetMaximumBin())
              pk = pyy.GetMaximum()
              thresh = pk - pyy.GetBinError(pyy.GetMaximumBin())
              lb = pyy.GetBinCenter(pyy.FindFirstBinAbove(thresh))
              ub = pyy.GetBinCenter(pyy.FindLastBinAbove(thresh))
              mod = 0.5*(lb+ub)
              erroryl = mod - lb 
              erroryh = ub - mod
              if erroryl <= pyy.GetBinWidth(bin)/2.0 : erroryl = pyy.GetBinWidth(bin)/2.0
              if erroryh <= pyy.GetBinWidth(bin)/2.0 : erroryh = pyy.GetBinWidth(bin)/2.0
              x.append(pxx.GetBinCenter(bin))
              y.append(mod)
              ely.append( erroryl )
              ehy.append( erroryh )
              #print(mod,ely[bin-1],ehy[bin-1])
              ex.append( myprof.GetBinWidth(bin)/2.0 )
              modgraph.SetPoint(modgraph.GetN(),x[-1],y[-1])
              modgraph.SetPointError(modgraph.GetN()-1,0.0,0.0,ely[-1],ehy[-1])
            #fexp = ROOT.TF1('fexp','expo + [2]',x[0],x[-1])
            modgraph.Fit(fexp,'NO')
            fpts = 500
            yfit = []
            xfit = []
            for pt in range(0,fpts):
              xfit.append( x[0] + pt*(x[-1] - x[0])/fpts )
              yfit.append( fexp.Eval( xfit[pt] ) )
            col_labels=['Fit']
            #row_labels=['$\chi^2$/ndf','Constant','Slope','Baseline']
            table_vals=[['{0:.3g}'.format(fexp.GetChisquare())+'/'+'{0:.3g}'.format(fexp.GetNDF())],
                        ['{0:.3g}'.format(fexp.GetParameter(0))+'$\pm$'+'{0:.3g}'.format(fexp.GetParError(0))],
                        ['{0:.3g}'.format(fexp.GetParameter(1))+'$\pm$'+'{0:.3g}'.format(fexp.GetParError(1))],
                        ['{0:.3g}'.format(fexp.GetParameter(2))+'$\pm$'+'{0:.3g}'.format(fexp.GetParError(2))]]
            the_table = plt.table(cellText=table_vals,
                                  rowLabels=row_labels,
                                  colLabels=col_labels,
                                  fontsize=8,
                                  colWidths=[0.3]*3,
                                  loc='upper right')
            #plt.title(tcut+' ; '+ str(avg_samples)+' samples/bin')
            plt.title(tcut+'\n'+ str(np.around(myprof.GetEntries()/myprof.GetNbinsX(),1))+' post-cut samples/bin')
            errormatrix=np.array([np.array(ely),np.array(ehy)])
            plt.errorbar(x,y,xerr=ex,yerr=errormatrix,markersize=4,fmt='ro')
            plt.plot(xfit,yfit,'b-')
            com = '((Tc-Ta)/log(an/cat)):(datime-'+str(t0)+')/3600.0'
            self.tree.Draw(com,tcut,'goff')
            cut_time = vtoa( self.tree.GetV2() , self.tree.GetSelectedRows() - 1 )
            cut_tau = vtoa( self.tree.GetV1() , self.tree.GetSelectedRows() - 1 )
            plt.annotate('Mode', xy=(0.1, 0.95), xycoords='axes fraction',color='red',weight='bold')
            ptsel = self.superimpose.get() 
            if ptsel == 1 :
                plt.plot(norm_time_clean,lf_time,'go',markersize=1.5,zorder=-32)           
                plt.plot(cut_time,cut_tau,'m^',markersize=1.5,zorder=-32)           
            if time_choice == 'day': plt.xticks(rotation=25)
            if time_choice == 'hour': plt.xlabel('Hours') 
            plt.ylabel('e$^{-}$ lifetime [$\mu$s]')

    def do_it(self) :
        print(self.analysis_option.get())

    def vtoa( self,buf, entries ):
        retarr = []
        for idx in range(0,entries-1) :
            retarr.append( buf[idx] )
        return retarr

    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.tree = ROOT.TTree('xpmdata','')
        self.myhist = ROOT.TH2F()
        self.myhist.SetName('myhist')


        # next two lines are for the texbox for entries
        self.dataFileInput = tk.Text( height=2, width=72) # text box( where user enters path)
        self.dataFileInput.insert(tk.END,os.getcwd() + os.sep + 'xpm_fitter_data' + os.sep + 'testData')
        self.dataFileInput.grid( row=2, column=0, columnspan=4,sticky=tk.W)

        self.analysis_option = tk.StringVar(master=self.parent)
        self.analysis_option.set('1 e- lifetime trend')
        self.analysis_pd_menu = tk.OptionMenu(self.parent,self.analysis_option,'1 e- lifetime trend','2 Cathode/Anode Trend','3 Power Monitors trend','4 Purity vs Cathode','5 Purity vs Power Monitors' )
        self.analysis_pd_menu.grid(row=3, column=0,columnspan=3,sticky=tk.W)
        
        self.plotting_option = tk.StringVar(master=self.parent)
        self.plotting_option.set('3 Scatter+Average')
        self.plotting_pd_menu = tk.OptionMenu(self.parent,self.plotting_option,'1 Average only','2 Scatter points only','3 Scatter+Average','4 Median (lognormal)','5 Mode')
        self.plotting_pd_menu.grid(row=4, column=0,columnspan=3,sticky=tk.W)
        
        self.opbLabel = tk.Label(height=1,width=30)
        self.opbLabel.config(text='Number of Samples to Average:')
        self.opbLabel.grid(row=5,column=0,sticky=tk.E)
        self.obsPerBin = tk.StringVar(self.parent)
        self.obsPerBin.set('10.0')  ### 33.0
        self.opb = tk.Spinbox(self.parent, increment=1.0, foreground='black', background='white', from_ = 1.0 , to = 1000000000.0 , width=24, textvariable = self.obsPerBin)
        self.opb.grid(row=5, column=1, sticky=tk.W)
        
        self.mdLabel = tk.Label(height=1,width=30)
        self.mdLabel.config(text='Minimum cathode-anode difference [mV]:')
        self.mdLabel.grid(row=6,column=0,sticky=tk.E)
        self.minacdiff = tk.StringVar(self.parent)
        self.minacdiff.set('0.2')  ### 33.0
        self.mdSpinBox = tk.Spinbox(self.parent, increment=0.01, foreground='black', background='white', from_ = 0.0 , to = 1000.0 , width=24, textvariable = self.minacdiff)
        self.mdSpinBox.grid(row=6, column=1, sticky=tk.W)

        self.binbyfibersave = tk.IntVar(value=0)
        self.bbfscheck = tk.Checkbutton( text='Bin by Fiber-save group?', variable = self.binbyfibersave, onvalue=1, offvalue=0 )
        self.bbfscheck.grid(row=7,column=0,sticky=tk.W)

        self.isrational = tk.IntVar(value=0)
        self.ircheck = tk.Checkbutton( text='Use rational function?', variable = self.isrational, onvalue=1, offvalue=0 )
        self.ircheck.grid(row=8,column=0,sticky=tk.W)

        self.superimpose = tk.IntVar(value=1)
        self.sicheck = tk.Checkbutton( text='Superimpose scatterplot?', variable = self.superimpose, onvalue=1, offvalue=0 )
        self.ircheck.grid(row=9,column=0,sticky=tk.W)

        self.fdLabel = tk.Label(height=1,width=30,font=('Arial',14))
        self.fdLabel.config(text='Time axis domain')
        self.fdLabel.grid(row=10,column=0,sticky=tk.E)
        self.fitdomain = tk.StringVar(self.parent)
        self.fitdomain.set(' ')  ### 33.0
        self.fd = tk.Entry(self.parent, font=('Arial',14),foreground='black', background='white', height=1,width=24,variable=self.fitdomain)
        self.fd.grid(row=10, column=1, sticky=tk.W)
        
        self.figure1 = Figure(figsize=(5, 4), dpi=100)
        self.canvas1 = FigureCanvasTkAgg(self.figure1, master=self.parent)
        self.plot_widget1 = self.canvas1.get_tk_widget()
        self.plot_widget1.grid(row=3, rowspan=7, column=2, columnspan=8)
        self.plt1 = self.figure1.add_subplot(111)
        plt.ion()
        self.toolbar = NavigationToolbar2Tk(self.canvas1, self.master, pack_toolbar=False )
        self.toolbar.update()
        self.toolbar.grid(row=11, column=2)
        self.canvas1.draw_idle()

        self.gobutton = tk.Button(text='Execute', command=self.do_it)
        self.gobutton.grid( row=12, column=1 )


