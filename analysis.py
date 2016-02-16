import csv
import numpy
from datetime import date, timedelta
import datetime
import pprint
import yaml
import sys

import parse
from parse import rsp
from parse import dateDiff
from parse import dtToStr
from parse import strToDT
from parse import cleanDT

import scipy.stats

#set up object that stores information for analysis
class analysisInfo(object):
    child_id=None
    start_date=None
    end_date=None
    std_time_for_bed=numpy.nan #std in time_for_bed
    std_event1_time=numpy.nan #std in event1 time
    std_time_to_event1=numpy.nan
    avg_time_to_event1=numpy.nan
    std_device_interval=numpy.nan
    std_device_start_time=numpy.nan
    std_offset_for_device_start=numpy.nan
    avg_offset_for_device_start=numpy.nan
    
    # fraction of days when reactgood rate is > rate_nr_tr
    reactgood_bigger_nt=None
    # pearson product-moment
    corr_dev_reactgood=None
    corr_dev_dev_total=None
    # spearman rank
    corr_spr_reactgood = None
    corr_spr_dev_total = None
    # spearman rank over sliding window
    corr_sprs_reactgood = None
    # works flag
    works_flag = None

    use_rate=None #use rate in 28 days prior to bigest improvement
    

    def __init__(self):
        self.use=[] # index of all days when app was used
        self.vibration_use=[] #index of all days vibrational device was used
        self.dist_event=[] #pair of numbers (days of occurance, number of dst_event)
        self.time_to_event1=[] # diff between time_for_bed and event1 time in minutes
        self.offset_for_device_start=[] #time diff between recommended device start time and actual time
        self.dev=dict() #for each reaction keep pair (days of occurance, number of occurance
        self.s_chart=dict() #(date as string)->chartRec
        # fields below have one element for each element of 'use'
        self.rate_dst_event=[]  #rate of distuption events computed in 7-day overlapping intervals
        self.biggest_improvement=[] #tuple with bigest improvement in rate of distuption events and the date when this improvement was seen , final rate
        self.devRate=dict() #key=reaction, value=rate of each reaction
        

class chartRec(object):
    """data for individual child sleep chart"""

    def __init__(self):
        self.date = None   # date as string
        # all number below are seconds since UTC midnight
        self.time_for_bed=[]  # black arrow down
        self.event1_time=[] # tuples (start, end)
        # end is derived from device event
        self.device_time=[] #event1_time+device_interval, random arrow
        self.device_vib=[] # tuples: (start, end, type)
        self.night_dsev_time=[] # up arrows



##**********output structure

def extract(children):
    #extracts info from children into analysisInfo() data structure
    #with most basic fields: child_id, start_date, end_date, initializes .dev
    analysis=dict()
    for ch in children.values():
        cid=ch.child_id
        if cid not in analysis:
            analysis[cid]=analysisInfo()
            analysis[cid].child_id=cid
        for d in ch.dates.values():
            if len(d.event1_time)>0:
                use_d=d.event1_time[0]
                if analysis[cid].start_date is None:
                    analysis[cid].start_date = use_d
                if analysis[cid].end_date is None:
                    analysis[cid].end_date = use_d
                analysis[cid].start_date=min(analysis[cid].start_date, use_d)
                analysis[cid].end_date=max(analysis[cid].end_date, use_d)
                for r in rsp:
                    analysis[cid].dev[r]=[]
        #print 'start date', analysis[cid].start_date
        #print 'end date', analysis[cid].end_date
    return analysis

def extend(children, analysis, start, end):
    #adds additional fields to analysisInfo()
    #such as time_for_bed (vector), event1_time (vector), dist_event, values for dev
    for ch in children.values():
        cid=ch.child_id
        for date, d in sorted(ch.dates.items()):
            t_bed=d.time_for_bed
            t_asleep=d.event1_time
            if len(t_bed)>0 or len(t_asleep)>0:
                if len(t_bed)>0:
                    analysis[cid].use.append(dateDiff(start, t_bed[0]))
                else:
                    analysis[cid].use.append(dateDiff(start, t_asleep[0]))
            #print "gg", analysis[cid].use
            dst_event=d.dev_time_cu
            num_tr=len(dst_event)
            if num_tr>0:
                if len(t_bed)>0:
                    analysis[cid].dist_event.append([dateDiff(start, t_bed[0]), num_tr])
                elif len(t_asleep)>0:
                    analysis[cid].dist_event.append([dateDiff(start, t_asleep[0]), num_tr])
                else:
                    dev_d=dtToStr((strToDT(dst_event[0])-datetime.timedelta(hours=8)))
                    #print '???', dst_event
                    analysis[cid].dist_event.append([dateDiff(start, dev_d), num_tr])
            #print 'device info:', d.device
            if d.device is not None:
                for r in d.device.react:
                    device_times=d.device.device_start_stop
                    #print device_times
                    num=len(device_times)
                    if num>0:
                        for times in device_times:
                            #print times
                            if len(t_bed)>0:
                                analysis[cid].dev[r].append([dateDiff(start, t_bed[0]),num ])
                            elif len(t_asleep)>0:
                                analysis[cid].dev[r].append([dateDiff(start, t_asleep[0]),num ])
                            else:
                                t_device=dtToStr((strToDT(device_times[0][0])-datetime.timedelta(hours=8)))
                                analysis[cid].dev[r].append([dateDiff(start, t_device),num ])
                if len(d.device.react)==0:
                    device_times=d.device.device_start_stop
                    if len(device_times)>0:
                        r='no react'
                        analysis[cid].dev[r]=[]
                        #print 'device no react', device_times
                        num=len(device_times)
                        if num>0:
                            for times in device_times:
                                #print times
                                if len(t_bed)>0:
                                    analysis[cid].dev[r].append([dateDiff(start, t_bed[0]),num ])
                                elif len(t_asleep)>0:
                                    analysis[cid].dev[r].append([dateDiff(start, t_asleep[0]),num ])
                                else:
                                    t_device=dtToStr((strToDT(device_times[0][0])-datetime.timedelta(hours=8)))
                                    analysis[cid].dev[r].append([dateDiff(start, t_device),num ])


def addSleepChart(children, analysis):
    #adds information about sleep chart to analysis
    for cid, crec in children.items():        
        cresult = analysis[cid]
        for date, dateinfo in crec.dates.items():
            chart_rec = chartRec()
            chart_rec.date = date
            cresult.s_chart[date] = chart_rec
            addSleepChartOneDate(date, dateinfo, chart_rec)
            
        if cid in []:  # 11, 57
            print
            yaml.dump_all([crec, cresult], sys.stdout)
            print

def addSleepChartOneDate(date, dateinfo, chart_rec):
    start = strToDT(date + ' 00:00:00')
    for time_str in dateinfo.time_for_bed:
        chart_rec.time_for_bed.append(
            (strToDT(time_str) - start).total_seconds()
        )
    event1_time=[]
    for time_str in dateinfo.event1_time:
        event1_time.append(
            (strToDT(time_str) - start).total_seconds()
        )
    #print "event1_time?", event1_time
    wakeup_times = []

    if dateinfo.device is not None:
        #print "checking", event1_time, chart_rec.time_for_bed
        for event1_time_sec, device_interval in zip(
                event1_time,
                dateinfo.device.shown_device_interval):
            if device_interval is None:
                # no prediction
                break
            assert event1_time_sec is not None
            chart_rec.device_time.append(
                event1_time_sec + device_interval * 60)

        for ldates, lreact in zip(
                dateinfo.device.device_start_stop,
                dateinfo.device.react):
            lstart_str, lstop_str = ldates
            lstart_sec = (strToDT(lstart_str) - start).total_seconds()
            chart_rec.device_vib.append([
                lstart_sec,
                (strToDT(lstop_str) - start).total_seconds(),
                lreact])
            wakeup_times.append(lstart_sec)

    for time_str in dateinfo.dev_time_cu:
        terr_time = (strToDT(time_str) - start).total_seconds()
        chart_rec.night_dsev_time.append(terr_time)
        wakeup_times.append(terr_time)

    wakeup_times.sort()

    for time_str in dateinfo.event1_time:
        asleep_start = (strToDT(time_str) - start).total_seconds()
        wakeup_times_after = [t for t in wakeup_times if t > asleep_start]
        if wakeup_times_after:
            asleep_end = wakeup_times_after[0]
        else:
            # no wakeup time. Assume 8 hours or end of day.
            asleep_end = min(asleep_start + 8 * 60 * 60,
                             24 * 60 * 60)
        chart_rec.event1_time.append([asleep_start, asleep_end])
    
def computeRateNtTr(analysis):
    #computes rate of distuption events in 7-day chunks with a sliding window
    for cid in analysis.keys():
        if len(analysis[cid].use)>=7:
            use=analysis[cid].use
            rate=[] #rate of distuption events
            nt_dict = dict(analysis[cid].dist_event)
            for i in range(len(use)):
                count=0
                for k in range(use[i] - 6, use[i] + 1):
                    count += nt_dict.get(k, 0)
                if i < 7:
                    rate.append(numpy.nan)
                else:
                    rate.append(1.0*count/7)

            analysis[cid].rate_dst_event = rate


def computeCompleteDevRate(analysis):
    #computes rate of device use in 7-day chunks with a sliding window
    for cid in analysis.keys():
        if len(analysis[cid].use)>=7:
            use=analysis[cid].use
            devInfo=dict()
            count=dict()
            rate=dict()
            for k in rsp:
                devInfo[k]=dict(analysis[cid].dev[k])
                rate[k]=[]
                for i in range(len(use)):
                    count[k]=0
                    for l in range(use[i] - 6, use[i] + 1):
                        count[k] += devInfo[k].get(l, 0)
                    if i < 7:
                        rate[k].append(numpy.nan)
                    else:
                        rate[k].append(1.0*count[k]/7)

            analysis[cid].complete_dev_rate =rate
            analysis[cid].dev_use_rate=numpy.zeros((len(rate[rsp[0]]),))
            for k in rsp:
                #print rate[k]
                analysis[cid].dev_use_rate+=numpy.array(rate[k])

def compBigImprov(analysis):
    #computes biggest improvement in rate of distuption events and the date for this improvement
    for cid in analysis.keys():
        current_max=None
        improve=[]
        if len(analysis[cid].use)>7:        #when for loop over dist event rate is done, select biggest improvement
            for r in analysis[cid].rate_dst_event:
                if not numpy.isnan(r):
                    if not current_max:
                        current_max=r
                    current_max=max(current_max, r)
                    improve.append(current_max-r)
                    #when for loop over dist event rate is done, select biggest improvement
            max_improve=max(improve)
            index_improve=improve.index(max_improve)+7 #since we skipped over first 7 days in rate_dst_event, which are nan
            date_improve=analysis[cid].use[index_improve]
            best_rate=analysis[cid].rate_dst_event[index_improve]
            analysis[cid].biggest_improvement=[max_improve, date_improve, best_rate]
            #print index_improve, analysis[cid].rate_dst_event
            #print 'improvement for id=', cid, analysis[cid].biggest_improvement



   

#prepare features for Machine Learning

def computeUseRate(analysis):
    #computes use rate for 28 days (if available) right before best improvement
    testPeriod=28 #days for improvement
    for ch in analysis.values():
        if len(ch.biggest_improvement)>0:
            bi_date=ch.biggest_improvement[1]
            start_d=bi_date-testPeriod
            end_d=bi_date+1
            count=0 #count number of days in ch.use between start_d and end_d
            for d in ch.use:
                if d>start_d and d<end_d:
                    count+=1
            period=min(28, bi_date-ch.use[0]+1)
            ch.use_rate=1.0*count/period
            #print 'best date', bi_date, 'all days', ch.use
            #print "use rate", ch.use_rate

def computeDevRate(analysis):
    #computes use of device rate for 28 days (if available) right before best improvement
    testPeriod=28 #days for improvement
    #print '*******looking at device us'
    count=dict()
    
    for ch in analysis.values():
        if len(ch.biggest_improvement)>0:
            for k in rsp:
                count[k]=0
            bi_date=ch.biggest_improvement[1]
            start_d=bi_date-testPeriod
            end_d=bi_date+1
            for k in ch.dev.keys():
                for data in ch.dev[k]:
                    if data[0]>start_d and data[0]<end_d:
                        count[k]+=data[1]
            for k in count.keys():
                ch.devRate[k]=1.0*count[k]/testPeriod
            #print ch.child_id, count 
            #print ch.devRate

def computeTimeToEvent1(analysis):
    #computes time a child needs to event 1 on each day
    #uses earliest time_for_bed and event1_time
    for cid in analysis.keys():
        for date in sorted(analysis[cid].s_chart.keys()):
            sc=analysis[cid].s_chart[date]
            data_exists=False
            if len(sc.time_for_bed)>0:
                if len(sc.event1_time)>0:
                    data_exists=True
                    #print sc.time_for_bed[0], sc.event1_time[0][0]
                    time_diff=sc.event1_time[0][0]-sc.time_for_bed[0]
                    if time_diff<0:
                        print 'for child_id=', cid, ' date=', date, 'mistaken time_for_bed and event1_time'
                    analysis[cid].time_to_event1.append(time_diff)
            if not data_exists:
                analysis[cid].time_to_event1.append(numpy.nan)
        #print "time to event 1", cid, analysis[cid].time_to_event1
        if numpy.nan in  analysis[cid].time_to_event1:
            event1_times=[t for t in  analysis[cid].time_to_event1 if not numpy.isnan(t)]
        else:
            event1_times= analysis[cid].time_to_event1
        analysis[cid].avg_time_to_event1=numpy.mean( event1_times)
        analysis[cid].std_time_to_event1=numpy.std( event1_times)
        #print 'time to event 1 stats:', analysis[cid].std_time_to_event1,  analysis[cid].avg_time_to_event1


def computeOffsetDeviceStartTime(analysis):
    #computes time difference between recommened device start time and actual device start time
     for cid in analysis.keys():
        for date in sorted(analysis[cid].s_chart.keys()):
            sc=analysis[cid].s_chart[date]
            data_exists=False
            if len(sc.device_time)>0:
                if len(sc.device_vib)>0:
                    data_exists=True
                    #print sc.time_for_bed[0], sc.event1_time[0][0]
                    time_diff=sc.device_vib[0][0]-sc.device_time[0]
                    analysis[cid].offset_for_device_start.append(time_diff)
                    #if time_diff>2000.:
                        #print "large of set:", cid, date, time_diff
            if not data_exists:
                analysis[cid].offset_for_device_start.append(numpy.nan)
        #print "device offset", cid, analysis[cid].offset_for_device_start
        if numpy.nan in  analysis[cid].offset_for_device_start:
            offset_for_device_start=[t for t in  analysis[cid].offset_for_device_start if not numpy.isnan(t)]
        else:
            offset_for_device_start= analysis[cid].offset_for_device_start
        analysis[cid].avg_offset_for_device_start=numpy.mean(offset_for_device_start)
        analysis[cid].std_offset_for_device_start=numpy.std(offset_for_device_start)
        #print 'std_offset_for_device_start stats:', analysis[cid].std_offset_for_device_start,  analysis[cid].avg_offset_for_device_start

def computeVibUse(analysis):
    #computes dates for vibration use
    
    for ch in analysis.values():
        days=set()
        for r in ch.dev.values():
            for el in r:
                days.add(el[0])
        #print 'days', days
        ch.vibration_use=list(days)

def spearman_safe(a, b):
    assert len(a) == len(b)
    if len(a) < 2: 
        return 0
    rho, p = scipy.stats.spearmanr(a, b)
    if numpy.isnan(rho):
        return 0
    return float(rho)

def computeCorrCoef(analysis):
    #computes corrCoef between dist_ev_rate and "reactgood" reaction
    print "correlation between night tr and devRate for ReactGood"
    cnt = 0
    for ch in analysis.values():
        if len(ch.use)<=20: continue
        #if ch.start_date<='2015-08-31': continue

        #print ch.complete_dev_rate['ReactGood']
        dist_event=[i for i in ch.rate_dst_event if not numpy.isnan(i)]
        reactgood_rec=[i for i in ch.complete_dev_rate['ReactGood'] if not numpy.isnan(i)]
        dev_used=[i for i in ch.dev_use_rate if not numpy.isnan(i)]
        assert len(dist_event) == len(reactgood_rec)
        assert len(dist_event) == len(dev_used)


        ch.reactgood_bigger_nt = sum(
            a > b for (a, b) in zip(reactgood_rec, dist_event)) * 1.0 / len(ch.use)

        ch.corr_spr_reactgood=spearman_safe(dist_event, reactgood_rec)
        ch.corr_spr_dev_total=spearman_safe(dist_event,dev_used)

        ch.corr_dev_reactgood=float(numpy.corrcoef(dist_event, reactgood_rec)[0][1])
        ch.corr_dev_dev_total=float(numpy.corrcoef(dist_event, dev_used)[0][1])

        sprs_move_vals = []
        # we slide window over every day in use
        for idx1, day1 in enumerate(ch.use):
            idx2 = idx1
            full_window = False
            while idx2 < len(ch.use):
                if ch.use[idx2] > (day1 + 21):
                    full_window = True
                    break
                idx2+=1
            if (not full_window) and (idx1 != 0):
                # we got to the end, and this is not a first record
                break            
            win_dist_event=dist_event[idx1:idx2]
            if len(win_dist_event) < 2:
                continue
            win_reactgood_rec=reactgood_rec[idx1:idx2]
            reactgood_bigger_cnt=sum(a > b for (a, b) in zip(win_reactgood_rec, win_dist_event))
            reactgood_bigger = reactgood_bigger_cnt > (len(win_reactgood_rec) / 2.0)
            sprs_move_vals.append((spearman_safe(win_dist_event, win_reactgood_rec),
                                   idx1, idx2, reactgood_bigger))
     
        # get value and borders for 'best' window
        ch.corr_sprs_reactgood, idx1, idx2, reactgood_bigger = min(sprs_move_vals)

        reactgood_bigger = ch.reactgood_bigger_nt > 0.5
        ch.works_flag = reactgood_bigger and ((ch.corr_sprs_reactgood + ch.corr_spr_reactgood) < 0)        
        #plt.figure()
        #plt.plot(dist_event, reactgood_rec, '*')

        #print ch.child_id, ch.corr_sprs_reactgood, sprs_move_vals
        #cnt += 1
        #if cnt > 10: fail        


def addFeatureForML(children, analysis):
    #adds features needed for machine learning, such as std_time_for_bed
    for cid in analysis.keys():
        time_for_beds=[]
        intervals=[]
        start_t=[]
        event1_times=[]
        for d in children[cid].dates.values():
            if len(d.time_for_bed)>0:
                tm=strToDT(d.time_for_bed[0]).time()
                time_for_bed=tm.second+tm.minute*60+tm.hour*3600
                #print time_for_bed
                time_for_beds.append(time_for_bed)
            if len(d.event1_time)>0:
                tm=strToDT(d.event1_time[0]).time()
                event1_time=tm.second+tm.minute*60+tm.hour*3600
                event1_times.append(event1_time)
            if hasattr(d.device, 'shown_device_interval'):
                if len(d.device.shown_device_interval)>0:
                    intervals.append(d.device.shown_device_interval[0])
            if hasattr(d.device,  'device_start_stop'):
                if len(d.device.device_start_stop)>0:
                    stm=strToDT(d.device.device_start_stop[0][0]).time()
                    stms=stm.second+stm.minute*60+stm.hour*3600
                    start_t.append(stms)
                    
        time_for_beds=numpy.array(time_for_beds)
        if len(time_for_beds)>0:
            analysis[cid].std_time_for_bed=numpy.std(time_for_beds)
        if len(event1_times)>0:
            analysis[cid].std_event1_time=numpy.std(event1_times)
            #print 'test event1 time', cid, analysis[cid].std_event1_time
        intervals=numpy.array(intervals)
        if len(intervals)>1:
            analysis[cid].std_device_interval=numpy.std(intervals)
        #print 'sdt intervals', analysis[cid].std_device_interval
        if len(start_t)>1:
            #print start_t
            analysis[cid].std_device_start_time=numpy.std(start_t)
    computeUseRate(analysis)
    computeDevRate(analysis)
    computeTimeToEvent1(analysis)   
    computeOffsetDeviceStartTime(analysis)
    computeVibUse(analysis)

##********analysis
def findStartDate(analysis):
    #find the strt date in the table
    start=None
    for cho in analysis.values():
        #print '****', cho.start_date
        if cho.start_date!=None:
            if not start:
                start=cho.start_date
            #print "33", start, cho.start_date
            start=min(start, cho.start_date)

    return start

def findEndDate(analysis):
    #find the end date in the table
    end=None
    for cho in analysis.values():
        #print '****', cho.start_date
        if cho.end_date!=None:
            if not end:
                end=cho.end_date
            #print "33", start, cho.start_date
            end=max(end, cho.end_date)

    return end

def userReport(analysis):
    #prints report about users' sussess using the device
    count_u=0
    count_0=0
    num_day_use=20
    for cid in analysis.keys():
        if len(analysis[cid].use)>num_day_use:
            count_u+=1
            if len(analysis[cid].biggest_improvement)>1:
                #print analysis[cid].biggest_improvement
                if analysis[cid].biggest_improvement[2]==0.0:
                    count_0+=1
            
    print '***Report***'
    print 'total number of users:', count_u, 'number users who got 0 rate:', count_0    
     
def searchMissingUsers(analysis):
    #search for missing users: users without start date (i.e. no time_for_bed of event1_time)
    #prints a lot of data
    missing_users=[]
    
    for cid in analysis.keys():
        if analysis[cid].start_date==None:
            if len(analysis[cid].use)==0:
                missing_users.append(cid)
    print "Missing users' report for", len(missing_users), 'users'
    count_early_users=0
    for row in csv.DictReader(open("nightData.csv")):
        child_id=int(row['child_id'])
        if child_id in missing_users:
            dt=cleanDT(row['created_at'])
            if dt<'2015-09-01':
                count_early_users+=0
            for k in row.keys():
                if row[k]!='':
                    print k, row[k],
            print
    print 'early users', count_early_users

def processData():
    #takes data from database and process it to produce analysisInfo with features that can be used for machine learning

    children=parse.parseData()
    analysis = extract(children)

    startT=findStartDate(analysis)
    endT=findEndDate(analysis)

    print "date diff", dateDiff(startT, endT)
    extend(children, analysis, startT, endT)
    addSleepChart(children, analysis)
    computeRateNtTr(analysis)
    compBigImprov(analysis)
    computeCompleteDevRate(analysis)
    userReport(analysis)
    #searchMissingUsers(analysis)
    addFeatureForML(children, analysis)
    computeCorrCoef(analysis)
    return children, analysis, startT, endT


#******Testing
#children, analysis, startT, endT=processData()
