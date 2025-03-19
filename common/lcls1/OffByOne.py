#!/usr/bin/env python
from __future__ import print_function
import psana
import numpy as np
import sys
import os
import logging
import argparse
from IPython import embed
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as mp
from collections import deque
fpath = os.path.dirname(os.path.abspath(__file__))
fpathup = "/".join(fpath.split("/")[:-1])
try:
    fpath = os.environ.get("MYDIR", fpathup).replace("/arp_scripts", "")
except:
    fpath = fpathup
sys.path.append(fpath)
from common.dropshotcalculations import calculations

class LCLS1_chi_group:

    def __init__(self, window): #Put variables here
        self.detslots={}
        self.msg="The PDFs are found in " + "/cpo/cporo"
        self.mean={}
        self.standdev={} #A bunch of empty dictionaries. Keywords are the detectors.
        self.meanerr={}
        self.maxdropshotcount=10
        self.window=window
        self.eventtime=deque("", window)
        self.eventcode=deque("", window)
        self.dets_in_keV = ['Epix10ka2M',
                            'Epix10ka',
                            'Epix10kaQuad',
                            'Epix100a',
                            'Jungfrau']
        self.areadets = [d for d in self.dets_in_keV]
        self.areadets.append('Opal1000')
        self.areadets.append('Zyla')
        self.areadets.append('Alvium')
        self.areadets.append('Rayonix')
        self.supported_dets = ['Wave8',
                               'EBeam',
                               'FEEGasDetEnergy',
                               'FEE-SPEC0',
                               'PhaseCavity',
                               'Pim', 'Ipm',
                               'BMMON', 'DIO']
        self.bigeventtimelist=[] #Used to store a bunch of event times

    def timestampgatherer(self):
        #Sets mode to small data mode, apparently because it's easier to scroll through
        ds = psana.DataSource(f"exp={args.experiment}:run={args.run}:smd") 
        ievr = 99
        for key in ds.env().configStore().keys():
            if key.src().__repr__().find("Evr") > 0:
                if int(key.src().__repr__()[-2]) < ievr:
                    ievr = int(key.src().__repr__()[-2])
                    evrName = (key.src().__repr__().split("(")[1])[:-1]
        if evrName is None:
            print("did not find EVR, this data is _weird_ ")
            return None
        self.evrdet = psana.Detector(evrName)
        for nevent,evt in enumerate(ds.events()):
            evtid=evt.get(psana.EventId) #...we get an event ID.
            if self.eventcode is None:
                continue
            seconds=int(evtid.time()[0])
            nanoseconds=evtid.time()[1] #Those event IDs get translated into timestamps for storage.
            fiducials=evtid.fiducials()
            self.eventcode.append(self.evrdet.eventCodes(evt))
            self.eventtime.append(psana.EventTime(int((seconds<<32)|nanoseconds), fiducials))
            if len(self.eventcode)<self.window:
                continue
            if 162 in self.eventcode[(window//2)]: #162 is the event code for a dropped shot. So if the middle of the dequeue is a dropped shot, it sticks the information in the dequeue into a list.
                self.bigeventtimelist.append(list(self.eventtime))
                continue
            if len(self.bigeventtimelist)==self.maxdropshotcount:
                print('Reached max drop shot count:',self.maxdropshotcount)
                break

    def detresultgetter(self):
        ds=psana.DataSource(f"exp={args.experiment}:run={args.run}:idx") #Switches the idx mode. More data or somethin'
        good_det_name=self.good_detectors() #Run the good_detectors function near the top. Returns a list.
        print('Supported detectors:',good_det_name) #Just a personal thing, where it spits out a list of detectors it's gonna go through.
        self.detlist=[psana.Detector(detname) for detname in good_det_name]
        run = next(ds.runs())
        for det in self.detlist:
            self.detslots[det.name]=[] #For every detector, it puts in a detector keyword with empty list "detname":[]
            for i in range(window):
                self.detslots[det.name].append([]) #Yoooooo, we got lists... IN LISTS. Perfect for slots
        #print(bigeventtimelist)
        for ndropshots, times in enumerate(self.bigeventtimelist): #Remember that bigeventtimelist are nested lists [[]]
            for islot,t  in enumerate(times): #Second layer of list
                event=run.event(t)
                if event==None:
                    print("None")
                    ipy.embed()
                    print('empty event')
                self.detresult(event, self.detlist, islot) #See function with corresponding name
        run.end()

    def calculations(self):
        for det in self.detlist:
            for nstuff, stuff in enumerate(self.detslots[det.name]):
                self.detslots[det.name][nstuff]=np.array(self.detslots[det.name][nstuff]) #Converts the nasty stuff into arrays for easy calculations

            self.mean[det.name]=[]
            self.standdev[det.name]=[]
            self.meanerr[det.name]=[]
            for n, slot in enumerate(self.detslots[det.name]): #Gettin' all the stuff for the chisquare calculations.
                self.mean[det.name].append(np.mean(slot))
                self.standdev[det.name].append(np.std(slot))
                self.meanerr[det.name].append(np.std(slot)/(len(slot)**.5))

    def detresult(self, event, detlist, islot): #Detresults translates the events into workable arrays and numbers
        for det in detlist:
            if det.name.dev in self.dets_in_keV:
                calib=det.calib(event)
                if calib is None:
                    continue
                #use 3 keV as threshold
                calib[calib<3.]=0
                result=np.sum(calib)
            elif det.name.dev in self.areadets:
                calib=det.calib(event)
                if calib is None:
                    continue
                result=np.sum(calib)
            elif "EBeam" == det.name.dev:
                ebeam=det.get(event)
                if ebeam is None:
                    continue
                result=ebeam.ebeamCharge()
            elif "FEEGasDetEnergy" == det.name.dev:
                feeg=det.get(event)
                if feeg is None:
                    continue
                result=feeg.f_11_ENRC()
            elif "PhaseCavity" == det.name.dev:
                phase=det.get(event)
                if phase is None:
                    continue
                result=phase.charge1()
            elif "BMMON" in det.name.dev or "DIO" in det.name.dev:
                bmmon=det.get(event)
                if bmmon is None:
                    continue
                result=bmmon.TotalIntensity()
            elif "Ipm" in det.name.dev or "Pim" in det.name.dev:
                val=det.channel(event)
                if val is None:
                    continue
                result=np.sum(val)
            elif "Wave8" in det.name.dev:
                wave=det.raw(event)
                if wave is None:
                    continue
                wavesum=wave[0].astype(float)
                for i in range(1,8):
                    wavesum+=wave[i]
                result=np.sum(wavesum)
            elif "FEE-SPEC0" in det.name.dev:
                spec_evt=det.get(event)
                #spec_list = [d for d in dir(spec_evt) if d[0]!='_']
                #print(spec_list)
                if spec_evt is None:
                    continue
                #spec = getattr(spec_evt,'hproj')()
                #speci = getattr(spec_evt,'integral')()
                #print('integral ',speci)
                #print('spec sum ',np.sum(spec.astype(float)))
                result=getattr(spec_evt,'integral')()
            else:
                print('Did not find detector:',det.name)
                continue
            self.detslots[det.name][islot].append(result)
            #print(self.detslots)
    def good_detectors(self):
        good_det=[]
        detnamelist=psana.DetNames("detectors")
        for detname in detnamelist:
            for supported_det in self.supported_dets:
                if supported_det in detname[0]:
                    good_det.append(detname[0])
            for areadet in self.areadets:
                if areadet in detname[0]:
                    good_det.append(detname[0])
        return good_det


# logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Define Args
parser = argparse.ArgumentParser()
parser.add_argument(
    "--run", help="run", type=str, default=os.environ.get("RUN_NUM", "")
)
parser.add_argument(
    "--experiment",
    help="experiment name",
    type=str,
    default=os.environ.get("EXPERIMENT", ""),
)
parser.add_argument("--postElog", help="post plot to elog", action="store_true")
parser.add_argument("--url", default="https://pswww.slac.stanford.edu/ws-auth/lgbk/")
parser.add_argument("--detnames", default="")
parser.add_argument("--pdf", help="Set to true if you want a series of graphs showing the means. Graphs are stored in the same directory this program is in, with the detector as the name.", action='store_true')
parser.add_argument("--email", help="Enter --email in order to get a useless email", action='store_true')
window=5

args = parser.parse_args()
logger.debug("Args to be used for data quality plots: {0}".format(args))
    
dropshotinfo=LCLS1_chi_group(window)
dropshotinfo.timestampgatherer()
ndrop = len(dropshotinfo.bigeventtimelist)
if ndrop>5:
    print('Found',ndrop,'dropped shots.')
else:
    print('Too few dropped shots:',ndrop)
    sys.exit(-1)
dropshotinfo.detresultgetter()

emailarg=args.email
pdfarg=args.pdf
calculationvar=calculations(dropshotinfo.detslots, window, emailarg, pdfarg)
calculationvar.calculations()
calculationvar.chideterminer()
#for debugging
print(dropshotinfo.detslots)
print(calculationvar.mean)

if args.pdf==True: #make thta holoviews/html!
    calculationvar.graph(expstring)
if args.email==True:
    calculationvar.mailer()
