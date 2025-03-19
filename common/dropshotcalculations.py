from __future__ import print_function
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

sys.path.append('/sdf/group/lcls/ds/tools/smalldata_tools/pedplot')
from smalldata_tools.utilities import postRunTable
from smalldata_tools.utilities import postElogMsg

class calculations:
    def __init__(self, detslots, window, emailarg, pdfarg, minchisq=0.2):
        self.pdf=pdfarg
        self.email=emailarg
        self.window=window
        self.detslots=detslots
        self.mean={}
        self.standdev={}
        self.meanerr={}
        self.deletelist=[]
        self.msg = ''
        self.full_msg = ''
        self.minchisq=minchisq
    def chisquare(self, points):
        wtdMeanVar = 1./np.sum([1/v[1] for v in points])
        wtdMean = np.sum([v[0]/v[1] for v in points])*wtdMeanVar
        chisq = np.sum([(v[0]-wtdMean)**2/v[1] for v in points])
        dof = len(points)-1
        return chisq/dof
    def calculations(self):
        for detname in self.detslots:
            for value in self.detslots[detname]:
                if len(value)==0:
                    self.deletelist.append(detname)
                    break
            if detname in self.deletelist:
                continue
            for nstuff, stuff in enumerate(self.detslots[detname]):
                self.detslots[detname][nstuff]=np.array(self.detslots[detname][nstuff]) #Converts the nasty stuff into arrays for easy calculations
            self.mean[detname]=[]
            self.standdev[detname]=[]
            self.meanerr[detname]=[]
            for n, slot in enumerate(self.detslots[detname]): #Gettin' all the stuff for the chisquare calculation
                self.mean[detname].append(np.mean(slot))
                self.standdev[detname].append(np.std(slot))
                if self.standdev[detname][-1]==0 or np.isnan(self.standdev[detname][-1]):
                    self.deletelist.append(detname)
                    print("Standard deviation for ", detname, " was zero, so the detector was wiped")
                    self.mean.pop(detname)
                    self.meanerr.pop(detname)
                    self.standdev.pop(detname)
                    break
                self.meanerr[detname].append(np.std(slot)/(len(slot)**.5))
        for deleted_det in self.deletelist:
            self.detslots.pop(deleted_det)
    def chideterminer(self): #Uses values from self.calculations
            for detname in self.detslots:
                #if detname in self.deletelist:
                    #continue
                chisquares=[]
                dropshots=[]
                minchisq=931059 # a large number
                dropslot=90
                sets=[]
                maxpoint=[]
                for i in range(self.window):
                    maxpoint.append([self.mean[detname][i], self.standdev[detname][i]**2])
                chimax=self.chisquare(maxpoint)
                if chimax<1: # all points statistically agree with each other
                    print("No result for:",detname,'with chisquare %5.2f:'%chimax)
                    self.long_msg+="\n No results for "+str(detname)
                    continue
                for drop in range(self.window):
                    points=[]
                    for i in range(self.window):
                        points.append([self.mean[detname][i], self.standdev[detname][i]**2])
                    dropshots.append(drop)
                    points.pop(drop)
                    chisq=self.chisquare(points)
                    #print(chisq, drop, detname)
                    if chisq<minchisq:
                        minchisq=chisq
                        dropslot=drop
                    #print(set)
                if minchisq<self.minchisq and dropslot==int(self.window/2): 
                    print("Correct dropped shot for:",detname)
                    self.long_msg+="\n Correct dropped slot for "+str(detname)
                else:
                    print("Probably incorrect:",detname,'dropped shot in slot',dropslot,'with overall chisquare %5.2f'%chimax,'and minimum chisq',minchisq)
                    self.msg+="\n The dropped slot for "+str(detname)+" is "+str(dropslot)+"  out of "+str(self.window)+ " dropped shots, while the chimax was "+str(chimax)+" and the chiminimum was "+str(minchisq)
            #build a full messages, starting with the detector that seems to be incorrect
            self.long_msg=self.mgs+self.long_msg
    def graph(self,expstring):
        for detname in self.detslots:
            plt.clf()
            fname = expstring+'_'+str(detname)+'.pdf'
            fname = fname.replace(':','_')
            with PdfPages(fname) as pp:
                xcount=np.array([])
                xcount=np.arange(self.window)
                #print(mean[detname])
                #print(xcount)
                plt.errorbar(xcount, self.mean[detname], yerr=self.meanerr[detname], fmt="bo")
                pp.savefig(plt.gcf())

    def poster(self, experiment, run=None):
        if self.msg!="":
            if experiment in ['xpp','xcs','mfx','cxi','mec']:
                instrument_elog=f{"{experiment[:3].upper()} Instrument"}
            else
                instrument_elog=f{"{experiment[:3]}_Instrument"}
            #postElogMsg(exp=experiment,
            #            run=run,
            #            msg=self.long_msg,
            #            tag="OFF-BY-SOME",
            #            title="DETECTOR ALIGNMENT INFO")
            #postElogMsg(exp=instrument_elog,
            #            run=run,
            #            msg=self.long_msg,
            #            tag="OFF-BY-SOME",
            #            title="DETECTOR ALIGNMENT INFO")
            postElogMsg(exp='Controls',
                        #run=run,
                        msg=self.long_msg,
                        tag="OFF-BY-SOME",
                        title="DETECTOR ALIGNMENT INFO")
