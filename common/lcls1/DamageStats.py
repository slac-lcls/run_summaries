#!/usr/bin/env python
########################################################################
## load tools and packages required for data handling and plotting
########################################################################
import sys
import os
import argparse
import logging
import requests
from typing import Optional, Dict, List
import mimetypes
from pathlib import Path

import numpy as np
import panel as pn
import holoviews as hv
from holoviews import dim


hv.extension("bokeh")
pn.extension()

try:
    basestring
except NameError:
    basestring = str
fpath = os.path.dirname(os.path.abspath(__file__))
fpathup = "/".join(fpath.split("/")[:-1])
try:
    fpath = os.environ.get("MYDIR", fpathup).replace("/arp_scripts", "")
except:
    fpath = fpathup
sys.path.append(fpath)
sys.path.append('/sdf/group/lcls/ds/tools/smalldata_tools/pedplot')

from smalldata_tools.utilities import evtt2Rt
from smalldata_tools.utilities import postRunTable
from smalldata_tools.utilities import postElogMsg

#this should also be shared among the data quality plots and/or in code where it runs consistenly
def postDetectorDamageMsg(
    detectors: list,
    exp: str,
    run: int,
    damageVars: dict,
    *,
    tag: str = "SUMMARY_EVENT_DAMAGE",
    title: str = "EVENT DAMAGE INFO",
    smd_dir: Optional[str] = None,
    post_thresh: float = 0.1,
    save_elog = False,
) -> None:
    """Post detector event damage info for a specified run to the eLog.

    Parameters
    ----------
    detectors (list[str]) Names of detectors to report on.
    exp (str) Experiment name.
    run (int) Run number. Usually the current run.
    tag (str) Optional. Tag for the event damage summary posts.
    title (str) Optional. Title for event damage summary posts.
    smd_dir (str) Optional. Alternative directory for smalldata HDF5 files.
    post_thresh (float) Optional. Damage threshold (as a percentage)
        required to post a message to eLog. At least 1 detector must pass
        the threshold to post to eLog. Only detectors passing the threshold
        will be included.
    """
    table_header: str = (
        '<thead><tr><th colspan="3">' f"<center>{title}</center>" "</th></tr></thead>"
    )
    table_body: str = (
        "<tbody><tr>"
        "<td><b><center>Detector</center></b></td>"
        "<td><b><center>Missing/Damaged Events</center></b></td>"
        "<td><b><center>Missing/Damaged</center></b></td></tr>"
    )

    post_msg: bool = False

    runtable_damage_dict={}
    for det_name in detectors:
        try:
            damage_var: np.ndarray = damageVars[det_name]
            damage: int = len(damage_var[damage_var == 0])
        except Exception:
            continue
        dmg_frac: float = damage / len(damage_var)
        runtable_damage_dict[f"damage percent {det_name}"] = dmg_frac*100.
        runtable_damage_dict[f"N damage {det_name}"] = damage
        if dmg_frac > post_thresh:
            post_msg = True
            det_entry: str = (
                f"<tr><td><center>{det_name}</center></td>"
                f"<td><center>{damage}</center></td>"
                f"<td><center>{dmg_frac:.2%}</center></td></tr>"
            )
            table_body += det_entry
    table_body += "</tbody>"
    msg: str = f'<table border="1">{table_header}{table_body}</table>'
    if post_msg and save_elog:
        postElogMsg(exp=exp, msg=msg, tag=tag, title=title)
    print(runtable_damage_dict)
    postRunTable(runtable_damage_dict, exp, run)


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
parser.add_argument(
    "--directory",
    help="directory to read files from (def <exp>/hdf5/smalldata)",
    default=None,
)
parser.add_argument("--postElog", help="post plot to elog", action="store_true")
parser.add_argument("--url", default="https://pswww.slac.stanford.edu/ws-auth/lgbk/")
parser.add_argument("--detnames", default="")
args = parser.parse_args()
logger.debug("Args to be used for data quality plots: {0}".format(args))

##############################################
## Setup Global parameters and run numbers ###
##############################################
save_elog = args.postElog
expname = args.experiment
run = int(args.run)

if int(os.environ.get("RUN_NUM", "-1")) > 0:
    requests.post(
        os.environ["JID_UPDATE_COUNTERS"],
        json=[{"key": "<b>Damage Stats </b>", "value": "Started"}],
    )

######################################

######################################
S3DF_BASE = Path("/sdf/data/lcls/ds/")
hutch = expname[:3]
dirname = (
    f"{S3DF_BASE}/{hutch}/{expname}/hdf5/smalldata"
)
if args.directory: dirname = args.directory
fname = "%s/%s_Run%04d.h5" % (dirname, expname, run)
import tables
if os.path.isfile(fname):
    fh5 = tables.open_file(fname).root
else:
    print(fname,' done not exist, exit')
    sys.exit()
    
### Get data & define axis title&ranges.
eventTimeRaw = getattr(fh5, 'event_time').read()
eventTimeDim = hv.Dimension(("eventTimeR", "relative event time"))
eventTimeR = evtt2Rt(eventTimeRaw)

#now I need to get keys....
detNamesAll = [dmg for dmg in dir(getattr(fh5, 'damage')) if dmg[0]!='_']
detNamesAll = [dn for dn in detNamesAll if dn not in ['ControlData','evr2']]

if args.detnames != "":
    print('Looking at specified detectors only: ',args.detnames)
    detNames = args.detnames.split(',')
else:
    detNames = detNamesAll

plots = []
from holoviews.operation.timeseries import rolling

damageDims={}
damageVars={}
for detname in detNames:
    try:
        damageDims[detname]= hv.Dimension((f"damage/{detname}", f"{detname} Present"))
        damageVars[detname] = getattr(fh5, f"damage/{detname}").read()
    except Exception:
        continue

for detname in detNames:
    try:
        damageDim = damageDims[detname]
        damageVar = damageVars[detname]
        damagePlot = rolling(
            hv.Curve(damageVar, vdims=[damageDim], label=f"{detname}").opts(
                axiswise=True, color=hv.Palette("Spectral")
            ),
            rolling_window=10,
        )
        plots.append(damagePlot)
    except Exception:
        continue

multiDamagePlot = hv.Overlay(plots).opts(
    xlabel="Event (rolling average of 10)",
    ylabel="Present (Yes/No)",
    title="Missing/Damaged Data",
)

########################
# Tabs construction and finish
########################
tabs = pn.Tabs(multiDamagePlot)

if int(os.environ.get("RUN_NUM", "-1")) > 0:
    requests.post(
        os.environ["JID_UPDATE_COUNTERS"],
        json=[{"key": "<b>Damage Stats </b>", "value": "Done"}],
    )

postDetectorDamageMsg(
    detectors=detNames,
    exp=expname,
    run=run,
    damageVars=damageVars,
    title=f"EVENT DAMAGE INFO - r{run:04d}",
    smd_dir=args.directory,
    post_thresh=0.1,  # Percentage threshold to post to eLog
    save_elog=save_elog,
)

if save_elog:
    from summaries.summary_utils import prepareHtmlReport

    pageTitleFormat = "DamageStats/Run{run:04d}"
    prepareHtmlReport(tabs, expname, run, pageTitleFormat)

    if int(os.environ.get("RUN_NUM", "-1")) > 0:
        requests.post(
            os.environ["JID_UPDATE_COUNTERS"],
            json=[{"key": "<b>Damage Stats</b>", "value": "Posted"}],
        )

