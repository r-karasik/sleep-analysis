import csv
import numpy
from datetime import date, timedelta
import datetime


#universal parameters i need throughout
LAST_DAY='2016-01-11'

#possible reactions to the device 
rsp=['None', 'PreUseDisruption', 'Disruption', 'ReactGood','ReactBad', 'no react']

class childInfo(object):
    child_id = None
    
    def __init__(self):
        self.dates=dict() #key=date string,, value =dateInfo()
        self.initial=None  # initInfo object
        
          
class dateInfo(object):
    def __init__(self):
        self.time_for_bed=[]
        self.event1_time=[]
        self.device=None # devInfo() instance
        self.dev_time_ch=[] # from nightData table

        self.dev_time_fd=[] # from users_reply 
        self.dev_time_cu=[] # cumulative


class devInfo(object):
    def __init__(self):
        self.device_start_stop=[] #append by adding [start, end] times of company's device
        self.react=[] # all reactions for a given ch at given night
        self.shown_device_interval=[]# all interval for a given child at a given night

class initInfo(object):
    month_count=None
    night_count=None
    event3=None
    bed_sleep=None
    bed_dst_event=None



#********datetime functions**********
def cleanDT(dt):
    #cleans a string with datetime object +removes .? seconds
    return dt.split('.')[0]

def strToDT(dt):
    #takes string and return datetime format
    return datetime.datetime.strptime(dt, '%Y-%m-%d %H:%M:%S')
    

def dtToStr(dt):
    #converts data in datetime format to string

    return dt.strftime('%Y-%m-%d %H:%M:%S')


def smallTimeDiff(dt1, dt2):
    #returns true if difference between two dates is less than 3 min
    time1=strToDT(dt1)
    time2=strToDT(dt2)
    #print abs(time1 - time2).total_seconds() 
    return abs(time1 - time2).total_seconds() < 3*60

def isSameNight(dt1, dt2):
    #returns true if difference between two dates is less than 12 hr
    time1=strToDT(dt1)
    time2=strToDT(dt2)
    #print abs(time1 - time2).total_seconds() 
    return (time1 - time2).total_seconds() < 12*60*60

def isValid(dt):
    #determines if the date is valid
    return int(dt[:4])>2000

def dateDiff(dt1, dt2):
    #computes date difference in days
    date1=strToDT(dt1)
    date2=strToDT(dt2)
    return (date2-date1).days


##****start reading input from tables

class SchemaInfo(object):
    NIGHT_TABLE = "nightData.csv"

schema = SchemaInfo()

def parseNight(children):
    #read input from "night" Table

    for row in csv.DictReader(open("nightData.csv")):
        cid=int(row['child_id'])
        ch=children.get(cid)
        if ch is None:
            ch=children[cid]=childInfo()
            ch.child_id=cid
        #start reading other information from the same row
        
        t_asleep=cleanDT(row['event1_time'])
        t_bed=cleanDT(row['time_for_bed'])
        t_dst_event=cleanDT(row['dist_event_time'])
        #check if data is consistent
        if t_asleep and t_bed:
            if strToDT(t_asleep).date()==strToDT(t_bed).date():
                ch_date=strToDT(t_asleep).date().isoformat()
                if ch_date not in ch.dates:
                    ch.dates[ch_date]=dateInfo()
                    if isValid(t_bed):
                        ch.dates[ch_date].time_for_bed=[t_bed]
                    if isValid(t_asleep):
                        ch.dates[ch_date].event1_time=[t_asleep]
                    ch.dates[ch_date].device=devInfo()
                present=False
                for bt in ch.dates[ch_date].time_for_bed:
                    if smallTimeDiff(t_bed, bt):
                        present=True
                if not present:
                    if isValid(t_bed):
                        ch.dates[ch_date].time_for_bed.append(t_bed)
                present=False
                for ast in ch.dates[ch_date].event1_time:
                    if smallTimeDiff(t_asleep, ast):
                        present=True
                if not present:
                    if isValid(t_bed):
                        ch.dates[ch_date].event1_time.append(t_bed)
                if t_dst_event:
                    if isValid(t_dst_event):
                        ch.dates[ch_date].dev_time_ch.append(t_dst_event)
                        #can potentially be from wrong night
            
                l_start=cleanDT(row['device_ev1_time'])
                l_end=cleanDT(row['device_ev2_time'])
                react=row['reaction_to_device'][16:]
                if l_end:
                    ch.dates[ch_date].device.device_start_stop.append([l_start, l_end])
                    if react!='':
                        ch.dates[ch_date].device.react.append(react)
                    sh_interval=row['shown_interval']
                    if sh_interval:
                        ch.dates[ch_date].device.shown_device_interval.append(int(sh_interval))


def parseChildren(children):
    #parse info from table "children" with initial info
    for row in csv.DictReader(open("users.csv")):
        ch = children.get(int(row['id']))
        if not ch:
            # No reports from this child
            continue
        if ch.initial==None:
            ch.initial=initInfo()
            ch.initial.month_count = int(row['initial_dist_events_nightData_month'])
            ch.initial.night_count = int(row['initial_dist_events_night'])
    
#*****fixing data in class structure
def dev_timeFix(children):
    #upon initial read dev_times can be placed incorrectly if they are reported from a different date
    #potentially still have issues since can have the same date, moving to same date
    for ch in children.values():
        for da in ch.dates.values():
            for t_nt in da.dev_time_ch:
                if len(da.time_for_bed)>0:
                    c_date=da.time_for_bed[0]
                    c_date_str=strToDT(c_date).date().isoformat()
                    if not isSameNight(c_date, t_nt):
                        new_date=strToDT(t_nt).date().isoformat()
                        if new_date not in ch.dates:
                            ch.dates[new_date]=dateInfo()
                        ch.dates[new_date].dev_time_ch.append(t_nt)
                        ch.dates[c_date_str].dev_time_ch.remove(t_nt)
                        print 'inconsistency in dist event reporting', c_date, t_nt


#********Transfering information from tables into data structure
def dsEventsFromFeedback(children):
    #takes information about distuption events from users_reply tables and puts it into structure
    #works after info from "nightData" have been loaded
    for row in csv.DictReader(open("users_reply.csv")):
        cid=int(row['child_id'])
        ch=children.get(cid)
        if ch is None:
            ch=children[cid]=childInfo()
            ch.child_id=cid
        dst_event=row['dev_time_yesterday']
        if dst_event:
            if len(dst_event)<19:
                dst_event+=':00'
                print 'check if time has the right format?', dst_event
            #assume sleep with dist event started at most 8 hr ago, determine that date
            dst_event=cleanDT(dst_event)
            if isValid(dst_event):
                nt_date1=(strToDT(dst_event)-datetime.timedelta(hours=8)).date().isoformat()
                nt_date2=strToDT(dst_event).date().isoformat()
                #print 'rec date, dist event:', nt_date, dst_event
                if nt_date1 in ch.dates:
                    if len(ch.dates[nt_date1].time_for_bed)>0:
                        found=False
                        for bed_t in ch.dates[nt_date1].time_for_bed:
                            if isSameNight(bed_t, dst_event):
                                ch.dates[nt_date1].dev_time_fd.append(dst_event)
                                found=True
                        if not found:
                            if nt_date2 in ch.dates:
                                if len(ch.dates[nt_date2].time_for_bed)>0:
                                    found2=False
                                    for bed_t in ch.dates[nt_date2].time_for_bed:
                                        if isSameNight(bed_t, dst_event):
                                            ch.dates[nt_date2].dev_time_fd.append(dst_event)
                                            found2=True
                                    if not found2:
                                        ch.dates[nt_date1].dev_time_fd.append(dst_event)
                elif nt_date2 in ch.dates:
                    if len(ch.dates[nt_date2].time_for_bed)>0:
                        found2=False
                        for bed_t in ch.dates[nt_date2].time_for_bed:
                            if isSameNight(bed_t, dst_event):
                                ch.dates[nt_date2].dev_time_fd.append(dst_event)
                                found2=True
                        if not found2:
                            ch.dates[nt_date1]=dateInfo()
                            ch.dates[nt_date1].dev_time_fd.append(dst_event)
                else:
                    ch.dates[nt_date1]=dateInfo()
                    ch.dates[nt_date1].dev_time_fd.append(dst_event)
            
            


def mergeDistEv(children):
    #merges night_tr from users_reply and night Tables
    for ch in children.values():
        for d in ch.dates.values():
            if len(d.dev_time_ch)>0 and len(d.dev_time_fd)>0:
                #merge
                d.dev_time_cu.extend(d.dev_time_fd)
                for dst_event in d.dev_time_ch:
                    missing=True
                    for dsev_time in d.dev_time_fd:
                        if smallTimeDiff(dst_event, dsev_time):
                            missing=False
                    if missing:
                        #print 'adding missing dst_event', dst_event
                        d.dev_time_cu.append(dst_event)
                    
            else:
                d.dev_time_cu.extend(d.dev_time_ch)
                d.dev_time_cu.extend(d.dev_time_fd)


def parseData():
    #parses data from tables
    children=dict()
    parseNight(children)
    parseChildren(children)
    dev_timeFix(children)
    dsEventsFromFeedback(children)
    mergeDistEv(children)

    return children

#testing
#parseData()
