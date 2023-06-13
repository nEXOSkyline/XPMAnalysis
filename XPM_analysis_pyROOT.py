import matplotlib
matplotlib.use('TkAgg')
import sys, os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import dates as dt
import datetime
import ROOT

# Check version of python, works with two or three
version = sys.version_info.major

#hardcoded stuff
#path = 'C:/Users/exouser/Desktop/xpm_fitter_data/'    
#path = './'    
path = '/mnt/c/Users/exouser/Desktop/xpm_fitter_data/'    

tree = ROOT.TTree('xpmdata','')
myhist = ROOT.TH2F()
myhist.SetName('myhist')

def myinput( pt ):
	v = sys.version_info.major
	if v==2 :
		return raw_input( pt )
	else :
		return input( pt )
	

def vtoa( buf, entries ):
  retarr = []
  for idx in range(0,entries-1) :
    retarr.append( buf[idx] )
  return retarr

def rfc868(in_time):
    out_time = []
    offset_time = in_time-126144000-2208988800
    for i,e in enumerate(offset_time):
        out_time.append(datetime.datetime.fromtimestamp(e))
    return out_time

def tm2hrs(in_time): return (in_time-in_time[0])/3600.

def smooth(in_arr, window):
    buf_arr = in_arr
    for i in range(window-1):
        in_arr = np.row_stack((in_arr,np.roll(buf_arr,int(i-window/2))))
    avg_arr = np.mean(in_arr,axis=0)
    err_arr = np.std(in_arr,axis=0)
    avg_arr[0] = np.mean(buf_arr[:100])
    err_arr[0] = np.std(buf_arr[:100])
    avrg_arr = avg_arr[::window]
    errs_arr = err_arr[::window]
    return avrg_arr,errs_arr

def cut(in_arr,time_arr,t_min,t_max):
    return in_arr[np.argmin(np.absolute(time_arr-t_min)):np.argmin(np.absolute(time_arr-t_max))]

#Function to get the file
def file_sel(version):

    if version == 3:
        file_name = str(myinput('Enter filename (Without extension) '))
    else: 
        file_name = str(myinput('Enter filename (Without extension) '))
    
    return file_name 
    
#function to set the date time objects
def time_step(version):
    print ('For the time axis we will be using: hours, if that is okay press enter. If you want to use days instead of hours press the letter "d" and then enter.')

    if version == 3:
        selection = str(myinput(''))
    else: 
        selection = str(myinput(''))
    if selection == 'd':
        return 'day'
    else:
        return 'hour'

# check if the file exists, if not give user one more chance
file_name = file_sel(version)
try:
  os.listdir(path)
except OSError:
  path = './'
file_path = path + file_name + '.txt'
if not os.path.isfile(file_path):
    print('file is not found, choose one of the following files in the XPM fitter data folder')
    print(os.listdir(path))
    file_name = file_sel(version)
    
#load data from file
raw_input_file = np.genfromtxt(path+file_name+'.txt',delimiter=',')

raw_input_file_clean = raw_input_file[raw_input_file[:,5]>0]

#Set time step for X axis
time_choice = time_step(version)
    
#raw_input_file_clean = raw_input_file_clean[raw_input_file_clean[:,4]/raw_input_file_clean[:,5]>1]
if time_choice == 'day':
    norm_time_clean = rfc868(raw_input_file_clean[:,7])
    norm_time = rfc868(raw_input_file[:,7])
elif time_choice == 'hour':
    norm_time_clean = tm2hrs(raw_input_file_clean[:,7])
    norm_time = tm2hrs(raw_input_file[:,7])
else: sys.exit('Uncorrect time selection')



while(1):
    print('\nPlot Displayer:\n1 - Purity Trend (e- lifetime)\n2 - Cathode/Anode Trend\n3 - Power Monitors trend\n4 - Purity vs Cathode\n5 - Purity vs Power Monitors \n6 - Exit\n')
    sel = str(myinput('---- ')) 
    if sel == str(1):
        #plt.clf()
        tree.Reset()
        tree.ReadFile(path+file_name+'.txt','Tc:Ta:TcRise:TaRise:cat:an:offst:datime:IR:UV:chi2:nAvg:nTrig',',')
        tree.Draw('(Tc-Ta)/log(an/cat):datime','','goff')
        lf_time = vtoa( tree.GetV1() , tree.GetSelectedRows() - 1 )
        t0 = tree.GetV2()[0]
        t1 = tree.GetV2()[tree.GetEntries()-1]
        tau = vtoa( tree.GetV1() , tree.GetEntries() )
        maxtau = np.max( tau )
        bin_by_fibersave = False
        avg_samples = int(str(myinput('Number of Samples to Average: ')))
        binning_mode = (str(myinput('Bin by Fiber-Save group (Y/N)?')))
        if binning_mode.find('Y') >=0 or binning_mode.find('y') >=0 :
            ##################### binning by Fiber Save group ################## 
            bin_by_fibersave = True
            timebucket = []
            tlow = []
            tcent = []
            for entry in range(0,tree.GetSelectedRows()-1) :
                ind = entry - 1
                if ind < 0 :
                    ind = 0
                dtl = tree.GetV2()[entry] - tree.GetV2()[ind]
                dtr = tree.GetV2()[entry+1] - tree.GetV2()[entry]
                if  dtl < 1800.0 :
                    timebucket.append(tree.GetV2()[entry])
                elif len(timebucket) < avg_samples :
                    timebucket.append(tree.GetV2()[entry])
                else :
                    try :
                        tlow.append( np.array(timebucket).min()-1.0 )
                        tlow.append( np.array(timebucket).max()+1.0 )
                        tcent.append( np.array(timebucket).mean() )
                        #print(len(timebucket),(tree.GetV2()[entry-1]-t0)/3600.0)
                        timebucket = []
                        timebucket.append(tree.GetV2()[entry])
                    except ValueError :
                        pass
            tlow.append( np.array(timebucket).min() )
            tlow.append( np.array(timebucket).max() )
            tcent.append( np.array(timebucket).mean() )
            tlow = (np.array(tlow) - t0)/3600.0
            tcent = (np.array(tcent) - t0)/3600.0
            print('number of bins',len(tlow))
            ####################################################################
            myhist = ROOT.TH2F('myhist','',len(tlow)-1,tlow,200,0.0,maxtau)            
        else :
            nbinsX = int( tree.GetEntries()/avg_samples )
            myhist=ROOT.TH2F('myhist','',nbinsX,0.0,(t1-t0)/3600.0,int(maxtau/100.0),0.0,maxtau)
        
        tree.Draw('(Tc-Ta)/log(an/cat):(datime-'+str(t0)+')/3600.0','','goff')
        norm_time_clean = vtoa( tree.GetV2() , tree.GetSelectedRows() - 1 )
        #myhist=ROOT.TH2F()
        #maxtau = 50000.0
        com = '(Tc-Ta)/log(an/cat):(datime-'+str(t0)+')/3600.0>>myhist'
        mindiff = float(str(myinput('Minimum cathode-anode difference [mV]: ')))
        tcut = 'UV>30 && cat>0 && an>0 && an<(cat-'+str(mindiff)+')'
        print(tree.Draw(com,tcut,'goff'))
        #if bin_by_fibersave :
        #    myhist.Rebin(int(avg_samples/10.0))
        myprof = myhist.ProfileX()
        myprof.SetMarkerStyle(20)
        #myprof.Draw('e')
        x = []
        y = []
        e = []
        ex = []
        for bin in range(1,myprof.GetNbinsX()+1) :
            if( myprof.GetBinEntries(bin)<=1 ): continue
            x.append(myprof.GetBinCenter(bin))
            y.append(myprof.GetBinContent(bin))
            e.append(myprof.GetBinError(bin))
            ex.append(myprof.GetBinWidth(bin)/2.0)
        
        fitfunc = myinput( 'Fitting function: exponential+baseline or rational (E/R)?' )
        
        fexp = ROOT.TF1('fexp','expo + [2]',x[0],x[-1])
        row_labels=['$\chi^2$/ndf','Constant','Slope','Baseline']
        if fitfunc=='R' :
            fexp = ROOT.TF1('fexp','[0]/([1]+[2]*x)',x[0],x[-1])
            fexp.SetParameters(myprof.GetBinContent(1)*2,2.0,0.005) 
            row_labels=['$\chi^2$/ndf','A','B','C']
        try :
          mi,ma = str(myinput('Enter space-separated time bounds of fit domain, lower bound first, or just hit ENTER: ')).split()
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
        llsel = str(myinput('1 - Average only\n2 - Scatter points only\n3 - Scatter+Average\n4 - Median (lognormal)\n5 - Mode\n'))
        if llsel == str(1):
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
            plt.xlabel('Hours')
            plt.annotate('Average', xy=(0.1, 0.95), xycoords='axes fraction',color='red',weight='bold')
            plt.ylabel('e$^{-}$ lifetime [$\mu$s]')
            #plt.show(block=False)
        elif llsel == str(2):
            tree.Draw(com,tcut,'goff')
            cut_time = vtoa( tree.GetV2() , tree.GetSelectedRows() - 1 )
            cut_tau = vtoa( tree.GetV1() , tree.GetSelectedRows() - 1 )
            plt.close()
            plt.plot(norm_time_clean,lf_time,'go',markersize=1.5,zorder=-32)           
            plt.plot(cut_time,cut_tau,'m^',markersize=1.5,zorder=-32)           
            plt.xlabel('Hours')
            plt.ylabel('e$^{-}$ lifetime [$\mu$s]')
            plt.annotate('Scatter plot only', xy=(0.1, 0.95), xycoords='axes fraction',color='red',weight='bold')
            #plt.show(block=False)
        elif llsel == str(3):
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
            tree.Draw(com,tcut,'goff')
            cut_time = vtoa( tree.GetV2() , tree.GetSelectedRows() - 1 )
            cut_tau = vtoa( tree.GetV1() , tree.GetSelectedRows() - 1 )
            plt.plot(norm_time_clean,lf_time,'go',markersize=1.5,zorder=-32)           
            plt.plot(cut_time,cut_tau,'m^',markersize=1.5,zorder=-32)           
            #plt.errorbar(subsample_norm_time,avg,yerr=err/np.sqrt(avg_samples-1),fmt='ro')
            #plt.title('Raw Lifetime and Lifetime Averaged Over %i Samples'%avg_samples)
            if time_choice == 'day': plt.xticks(rotation=25)
            if time_choice == 'hour': plt.xlabel('Hours') 
            plt.ylabel('e$^{-}$ lifetime [$\mu$s]')
            plt.annotate('Average', xy=(0.1, 0.95), xycoords='axes fraction',color='red',weight='bold')
            plt.show(block=False)
        elif llsel == str(4):
            plt.close()
            tree.Draw('log((Tc-Ta)/log(an/cat)):datime','','goff')
            t0 = tree.GetV2()[0]
            t1 = tree.GetV2()[tree.GetEntries()-1]
            lnhist=ROOT.TH2F()
            tau = vtoa( tree.GetV1() , tree.GetEntries() )
            maxtau = np.max( tau )
            #maxtau = 50000.0
            lnhist=ROOT.TH2F('lnhist','',nbinsX,0.0,(t1-t0)/3600.0,int(np.exp(maxtau)/100.0),0.0,maxtau)
            com = 'log((Tc-Ta)/log(an/cat)):(datime-'+str(t0)+')/3600.0>>lnhist'
            tcut = 'UV>30 && cat>0 && an>0 && an<(cat-'+str(mindiff)+')'
            print(tree.Draw(com,tcut,'goff'))
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
            tree.Draw(com,tcut,'goff')
            cut_time = vtoa( tree.GetV2() , tree.GetSelectedRows() - 1 )
            cut_tau = vtoa( tree.GetV1() , tree.GetSelectedRows() - 1 )
            ptsel = str(myinput('Superimpose scatterplot (Y/N)?'))
            if ptsel == 'Y' :
                plt.plot(norm_time_clean,lf_time,'go',markersize=1.5,zorder=-32)           
                plt.plot(cut_time,cut_tau,'m^',markersize=1.5,zorder=-32)           
            if time_choice == 'day': plt.xticks(rotation=25)
            if time_choice == 'hour': plt.xlabel('Hours') 
            plt.annotate('Median (assumes lognormality)', xy=(0.1, 0.95), xycoords='axes fraction',color='red',weight='bold')
            plt.ylabel('e$^{-}$ lifetime [$\mu$s]')
        elif llsel == str(5):
            plt.close()
            tree.Draw('log((Tc-Ta)/log(an/cat)):datime','','goff')
            myhist=ROOT.TH2F()
            myhist=ROOT.TH2F('myhist','',nbinsX,0.0,(t1-t0)/3600.0,int(75*maxtau/100000.0),0.0,maxtau)
            print(myhist.ProjectionY('py_t',1,1).GetBinWidth(1))
            tree.Draw(com,tcut,'goff')
            x = []
            y = []
            ely = []
            ehy = []
            ex = []
            pxx = myhist.ProjectionX('pxx',1,2)
            modgraph=ROOT.TGraphAsymmErrors()
            for bin in range(1, myhist.GetNbinsX()+1 ) :
              if( myprof.GetBinEntries(bin)<=1 ): continue
              pyy = myhist.ProjectionY('pyy',bin,bin)
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
            tree.Draw(com,tcut,'goff')
            cut_time = vtoa( tree.GetV2() , tree.GetSelectedRows() - 1 )
            cut_tau = vtoa( tree.GetV1() , tree.GetSelectedRows() - 1 )
            plt.annotate('Mode', xy=(0.1, 0.95), xycoords='axes fraction',color='red',weight='bold')
            ptsel = str(myinput('Superimpose scatterplot (Y/N)?'))
            if ptsel == 'Y' :
                plt.plot(norm_time_clean,lf_time,'go',markersize=1.5,zorder=-32)           
                plt.plot(cut_time,cut_tau,'m^',markersize=1.5,zorder=-32)           
            if time_choice == 'day': plt.xticks(rotation=25)
            if time_choice == 'hour': plt.xlabel('Hours') 
            plt.ylabel('e$^{-}$ lifetime [$\mu$s]')
    elif sel == str(2):
        f, axarr = plt.subplots(2, sharex=True)
        axarr[0].plot(norm_time, raw_input_file[:,4],'co')
        axarr[0].set_ylabel('Cathode Amplitude [mV]')
        axarr[1].plot(norm_time, raw_input_file[:,5],'go')
        axarr[1].set_ylabel('Anode Amplitude [mV]')
        if time_choice == 'day': plt.xticks(rotation=25)
        if time_choice == 'hour': plt.xlabel('Hours')

    elif sel == str(3):
        plt.close()
        plt.plot(norm_time,raw_input_file[:,8],'b+',markersize=2,label='blue: IR')
        plt.plot(norm_time,raw_input_file[:,9],'m+',markersize=2,label='red: UV')
        plt.legend()
        plt.ylim(0,1000)
        if time_choice == 'day': plt.xticks(rotation=25)
        if time_choice == 'hour': plt.xlabel('Hours')
        plt.ylabel('Power Monitors [au]')

    elif sel == str(4):
        plt.ion()
        lf_time = (raw_input_file_clean[:,1]-raw_input_file_clean[:,0])/np.log(raw_input_file_clean[:,4]/raw_input_file_clean[:,5])
        plt.scatter(raw_input_file_clean[:,5],lf_time)
        plt.xlabel('Cathode Amplitude [mV]')
        plt.ylabel('e$^{-}$ lifetime [$\mu$s]')
        #plt.show(block=False)
        ct = str(myinput('1 - Input Time Cuts\nAny Keys to Go Back\n'))
        if ct == str(1):
          mi,ma = str(myinput('Enter min&Max: ')).split()
          mI = float(mi)
          Ma = float(ma)
          plt.clf()
          plt.scatter(cut(raw_input_file_clean[:,5],norm_time_clean,mI,Ma),cut(lf_time,norm_time_clean,mI,Ma),color='green')
          plt.xlabel('Cathode Amplitude [mV]')
          plt.ylabel('e$^{-}$ lifetime [$\mu$s]')
          #plt.show(block=False)

    elif sel == str(5):
        print('hello')
        plt.close()
        tree.Reset()
        tree.ReadFile(path+file_name+'.txt','Tc:Ta:TcRise:TaRise:cat:an:offst:datime:IR:UV:chi2:nAvg:nTrig',',')
        tree.Draw('UV:cat','cat>0','goff')
        cat_sig = vtoa( tree.GetV2() , tree.GetSelectedRows() - 1 )
        uv_sig = vtoa( tree.GetV1() , tree.GetSelectedRows() - 1 )
        plt.ion()
        #lf_time = (raw_input_file_clean[:,1]-raw_input_file_clean[:,0])/np.log(raw_input_file_clean[:,4]/raw_input_file_clean[:,5])
        #plt.scatter(raw_input_file_clean[:, 8],lf_time, color = 'blue', label='blue: IR')
        #plt.scatter(raw_input_file_clean[:, 9],lf_time, color = 'red', label='red: UV')
        plt.plot(cat_sig,uv_sig,'ro', markersize=1.5, label='red: UV')
        
        plt.xlabel('Cathode signal [mV]')
        plt.ylabel('UV signal [au]')
        plt.legend()
        #plt.show(block=False)

    elif sel == str(6): break

    else: print('Press 1-5 or 6 to exit')
    plt.show(block=False)
'''    elif sel == str(6):
        lf_time = (raw_input_file[:,1]-raw_input_file[:,0])/np.log(raw_input_file[:,4]/raw_input_file[:,5])
        avg,err = smooth(lf_time,2)
        buff = raw_input_file[:,5]
        buff = buff[::2]
        plt.scatter(buff,err)'''
        

