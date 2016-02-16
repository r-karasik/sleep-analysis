import matplotlib.pyplot as plt
import matplotlib
import numpy
import psycopg2
import sys
import yaml
import json

#import my code from files
import parse
import analysis



#variable to be used throughout
#color codes for device reaction
COLOR_CODES = {
    'None': 'yellow',
    'ReactGood': 'green',
    'ReactBad': 'blue',
    'Disruption': 'red',
    'PreUseDisruption': 'orange'
}


#**************************
#plotting
#*******************************



#*****************************
#plots with matplotlib and for the webapp that show dist event rates and device use rate
def plotRate(children, analysis):
    #plots rate of distuption events for users who used device for more than 21 days
    #plots done with matplotlib and saved as .jpg's in folder "figures2"
    count=0
    for cid in analysis.keys():
        if len(analysis[cid].use)>20:
            fig, ax=plt.subplots()
            ax.set_title('rate of distuption events for the last 7 days')
            ax.set_xlabel('number of days')
            ax.set_ylabel('dist event rate')
            rate=analysis[cid].rate_dst_event
            use=analysis[cid].use
            #print rate
            #print rate, use
            ax.plot(use, rate, 'bs-')
            init=children[cid].initial
            value=init.month_count*init.night_count/30.0
            ax.plot([min(use) - 1],value, 'ro')
            ax.set_xlim((-2, max(use)+1))
            ax.set_ylim((-0.2, max([value] + rate)*1.2))
            fig.savefig('figures2/dev_rate'+str(cid)+'.jpg')
            plt.close(fig)
            #plt.show()
            count+=1
            #if count>10:
            #    break


def plotAllRates(children, analysis):
    #plots rate of distuption events for users who used device for more than 21 days
    #plots the rate how often device was used with best possible reaction
    #figures are created with matplotlibe and saved as .jpg files in folder "figures2"
    count=0
    for cid in analysis.keys():
        if len(analysis[cid].use)>20:
            fig, ax=plt.subplots()
            ax.set_title('rate of events for the last 7 days')
            ax.set_xlabel('number of days')
            ax.set_ylabel('rate')
            rate=analysis[cid].rate_dst_event
            use=analysis[cid].use
            reactgood_rate=analysis[cid].complete_dev_rate['ReactGood']
            dev_rate=analysis[cid].dev_use_rate
            #print rate
            #print rate, use
            ax.plot(use, rate, 'bs-')
            ax.plot(use, dev_rate, 'k,-')
            ax.plot(use, reactgood_rate, 'go-')
            init=children[cid].initial
            value=init.month_count*init.night_count/30.0
            ax.plot([min(use) - 1],value, 'ro')
            ax.set_xlim((-2, max(use)+1))
            ax.set_ylim((-0.2, max([value] + rate)*1.2))
            #fig.savefig('figures2/dev_rate'+str(cid)+'.jpg')
            #plt.close(fig)
            plt.show()
            count+=1
            if count>10:
                break


#******************************
#plots for sleep chart

def plotSleepChart(crec, cresult, f_name):
    #plots a slep chart for one child and shows event 3 time, event 1 time, time for distuption events, time when device was used, etc.
    fig, ax = plt.subplots()

    #yaml.dump_all([crec, cresult], sys.stdout)

    dates = list()
    for date, chart_rec in sorted(cresult.s_chart.items()):
        n = len(dates)
        dates.append(date)
        y_pos = len(cresult.s_chart) - n - 1

        for tm in chart_rec.time_for_bed:
            ax.arrow(tm / 3600.0, y_pos - 0.1, 0, 0.8, width=0.01, 
                     length_includes_head=True, fc='k', ec='k')

        for tm in chart_rec.night_dsev_time:
            ax.arrow(tm / 3600.0, y_pos - 0.1 + 1, 0, -0.8, width=0.01, 
                     length_includes_head=True, 
                     fc='red', ec='red')

        for tm in chart_rec.device_time:
            x_val = tm / 3600.0
            ax.plot([x_val, x_val], [y_pos + 0.1, y_pos + 0.9],
                    color='magenta')

        x_vals = []
        for tm1, tm2 in chart_rec.event1_time:
            x_vals.append([tm1 / 3600.0,
                           (tm2 - tm1) / 3600.0]) 
        if x_vals:
            ax.broken_barh(x_vals, [y_pos + 0.3, 0.5], facecolors='blue')
        
        x_vals = []
        color_vals = []
        for tm1, tm2, ttype in chart_rec.device_vib:
            duration = (tm2 - tm1) / 3600.0
            x_vals.append([tm1 / 3600.0,
                           max(0.02, duration)]) 
            color_vals.append(COLOR_CODES.get(ttype, 'magenta')) 
        if x_vals:
            ax.broken_barh(x_vals, [y_pos + 0.1, 0.8], 
                           facecolors=color_vals, edgecolors=color_vals)

        
    ax.set_xlim(0,24)
    ax.set_xlabel('time (hours)')
    ax.set_xticks(range(0, 25, 4))

    rangeLabel=range( len(dates))
    rangeLabel.reverse()
    rangeLabelTicks=numpy.array(rangeLabel)
    rangeLabel=numpy.array(rangeLabel)+0.5
    #print rangeLabel
    ax.set_ylim(0, len(dates))
    ax.set_yticks(rangeLabelTicks, minor=True)
    ax.set_yticks(rangeLabel, minor=False)
    ax.set_yticklabels(dates)

    ax.grid(True)
    ax.yaxis.grid(False, which='major')
    ax.yaxis.grid(True, which='minor')

    if f_name is None:
        plt.show()
    else:
        fig.savefig(f_name)
    plt.close(fig)

def plotAllSleepCharts(children, analysis):
    #plots all sleep charts
    for cid, cresult in sorted(analysis.items()):
        crec = children[cid]
        plotSleepChart(crec, cresult, None)
        # if cid == 63:
        #     plotSleepChart(crec, cresult, None)
        #     fail

####***************************
#exploratory plots
#*****************************


#plots for CDF
def computeAndPlotCDF(analysis, n, end):
    #computes cdf for number of days people used the app
    userCount=numpy.zeros((n,1))
    for ch in analysis.values():
        index=len(ch.use)
        userCount[index]+=1
    cdf=numpy.zeros((n,1))
    for i in range(1,n):
        cdf[i]=sum(userCount[:i])
    diff=numpy.zeros((n,1))
    for i in range(1,n):
        diff[i]=cdf[i]-cdf[i-1]
    print 'all users',  numpy.transpose(diff[:50])
    print 'max', max(cdf), sum(userCount)
    plt.figure()
    plt.bar(range(n), cdf*1.0/sum(userCount), align='center')
    plt.title('cdf: number of people who stopped using the app')
    plt.xlabel('days')
    plt.ylabel('% of all users')
    plt.figure()

    plt.plot(diff, 'g*--')
    plt.title('number of people lost')
    plt.ylabel('number of people')
    plt.xlabel('days')
    #assume users are not lost until "cutOff" days after last use
    cutOff=10
    userCount_old=numpy.zeros((n,1))
    new_userCount=numpy.zeros((cutOff,1))
    for ch in analysis.values():
        if ch.start_date!=None:
            #print ch.end_date, end
            if dateDiff(ch.start_date, end)>=cutOff:
                if ch.start_date>'2015-08-31':
                    index=len(ch.use)
                    userCount_old[index]+=1
            else:
                userCount_old[-1]+=1
                index=len(ch.use)
                if index>10:
                    print ch.start_date
                new_userCount[index]+=1
        elif len(ch.use)>0:
            print 'issue with data, id=', ch.child_id
        else:
            index=len(ch.use)
            userCount_old[index]+=1
    print 'number of all users', len(analysis), 'number of old users', sum(userCount_old), sum(userCount)
    cdf_old=numpy.zeros((n,1))
    cdf_new=numpy.zeros((cutOff,1))
    for i in range(1,n):
        cdf_old[i]=sum(userCount_old[:i])
    for i in range(1, cutOff):
        cdf_new[i]=sum(new_userCount[:i])
    diff_old=numpy.zeros((n,1))
    for i in range(1,n):
        diff_old[i]=cdf_old[i]-cdf_old[i-1]
    diff_new=numpy.zeros((cutOff,1))
    for i in range(1, cutOff):
        diff_new[i]=cdf_new[i]-cdf_new[i-1]
    print 'new users', userCount_old[-1], 'old diff', numpy.transpose(diff_old[:50])
    print 'new diff', numpy.transpose(diff_new)
    plt.figure()
    plt.bar(range(n), cdf_old/sum(userCount_old), align='center')
    plt.title('cdf: number of people who stopped using device for >'+str(cutOff))
    plt.xlabel('days')
    plt.ylabel('number of people')

    plt.figure()
    #plt.plot(range(n), 1-cdf*1.0/sum(userCount_old),'b-', ls='steps',linewidth=4.0)
    plt.plot(range(n), 1-cdf_old/sum(userCount_old),'b-',ls='steps', linewidth=4.0)
    plt.title('Remaining Users')
    plt.xlabel('days')
    plt.xlim((0,50))
    plt.ylabel('fraction of all users')
    #matplotlib.rcParams.update({'font.size': 22})

    plt.figure()
    plt.plot(diff_old, 'g*--')
    plt.title('number of people lost for >'+str(cutOff))
    plt.ylabel('number of people')
    plt.xlabel('days')
    plt.show()
   

def computeAndPlotCDFforVibUse(analysis, n, end):
    #computes cdf for number of days people used vibration on the device

    earlyDate='2015-08-31' #assume data from users is good after this date

    #assume users are not lost until cutOff number of days after last use
    cutOff=10
    userCount_old=numpy.zeros((n,1))
    new_userCount=numpy.zeros((cutOff,1))
    for ch in analysis.values():
        if ch.start_date!=None and ch.end_date>earlyDate:
            #print ch.end_date, end
            if dateDiff(ch.start_date, end)>=cutOff:
                #if ch.start_date>'2015-09-01':
                index=len(ch.vibration_use)
                userCount_old[index]+=1
            else:
                index=len(ch.vibration_use)
                if index>cutOff:
                    print "problem with vibration use logic", ch.start_date
                new_userCount[index]+=1
        elif len(ch.vibration_use)>0:
            print 'issue with data, id=', ch.child_id
        else:
            index=len(ch.vibration_use)
            userCount_old[index]+=1
    print  'number of old users using vibration', sum(userCount_old)
    cdf_old=numpy.zeros((n,1))
    cdf_new=numpy.zeros((cutOff,1))
    for i in range(1,n):
        cdf_old[i]=sum(userCount_old[:i])
    for i in range(1, cutOff):
        cdf_new[i]=sum(new_userCount[:i])
    diff_old=numpy.zeros((n,1))
    for i in range(1,n):
        diff_old[i]=cdf_old[i]-cdf_old[i-1]
    diff_new=numpy.zeros((cutOff,1))
    for i in range(1, cutOff):
        diff_new[i]=cdf_new[i]-cdf_new[i-1]
    print 'new users', userCount_old[-1], 'old diff', numpy.transpose(diff_old[:50])
    print 'new diff', numpy.transpose(diff_new)
    plt.figure()
    plt.bar(range(n), cdf_old/sum(userCount_old), align='center')
    plt.title('cdf: number of people who stopped using the device')
    plt.xlabel('days')
    plt.ylabel('% of all users')
    plt.figure()

    plt.plot(diff_old, 'g*--')
    plt.title('number of people lost for >'+str(cutOff))
    plt.ylabel('number of people')
    plt.xlabel('days')
    plt.show()

#******************************************
#feature relations, correlations?
def plotUseRateVsBestRate(analysis):
    #shows plot of use rate vs best dist event rate from biggest improvement
    use_rate=[]
    best_rate=[]
    for ch in analysis.values():
        if len(ch.use)>20:
            if len(ch.biggest_improvement):
                use_rate.append(ch.use_rate)
                best_rate.append(ch.biggest_improvement[2])
    plt.figure()
    plt.title('best rate vs use freq in 4 weeks prior to best rate')
    plt.xlabel('best rate')
    plt.ylabel('use freq in 4 weeks')
    plt.plot(best_rate, use_rate, 'bo')
    plt.show()


def plotUseRateVsImprovRate(analysis):
    #shows plot of use rate vs biggest improvement diff rate for distuption events
    use_rate=[]
    improv_diff=[]
    best_rate=[]
    for ch in analysis.values():
        if len(ch.use)>20:
            if len(ch.biggest_improvement):
                use_rate.append(ch.use_rate)
                improv_diff.append(ch.biggest_improvement[0])
                best_rate.append(ch.biggest_improvement[2])


    imdiffLp6=[]
    imdiffBp6=[]
    brateLp6=[]
    brateBp6=[]
    for i in range(len( use_rate)):
        if use_rate[i]<.6:
            imdiffLp6.append(improv_diff[i])
            brateLp6.append(best_rate[i])
        else:
            imdiffBp6.append(improv_diff[i])
            brateBp6.append(best_rate[i])

    print '******Report 2******'
    print 'for users with use rate <.6, improvement stats: mean', numpy.mean(imdiffLp6), 'std', numpy.std(imdiffLp6)
    print 'for users with use rate <.6, best rate stats: mean', numpy.mean(brateLp6), 'std', numpy.std(brateLp6)
    print 'for users with use rate >=.6, improvement stats: mean', numpy.mean(imdiffBp6), 'std', numpy.std(imdiffBp6)
    print 'for users with use rate >=.6, best rate stats: mean', numpy.mean(brateBp6), 'std', numpy.std(brateBp6)

    plt.figure()
    plt.title('improvement difference vs use freq in 4 weeks prior to best rate')
    plt.xlabel('improvement difference')
    plt.ylabel('use freq in 4 weeks')
    plt.plot(improv_diff, use_rate, 'bo')
    plt.show()

def plotDevRateVsImprov(analysis):
    #plots dev rate vs biggest improvement in dist event rate
    improvement=[]
    reactgood_rate=[]
    for ch in analysis.values():
        if len(ch.biggest_improvement)>0:
            improvement.append(ch.biggest_improvement[0])
            reactgood_rate.append(ch.devRate['ReactGood'])
    plt.figure()
    plt.plot(reactgood_rate, improvement, 'gs')
    plt.title('sussessful device use rate vs improvement')
    plt.ylabel('improvement in distuption events')
    plt.xlabel('rate of correct reaction')
    plt.show()



def plotCorrVsOffset(analysis):
    #plots correlation between dist event to reactgood and offset for start of company's device
    correlation=[]
    std_offset=[]
    mean_offset=[]
    for ch in analysis.values():
        if len(ch.use)>0:
            reactgood= ch.corr_dev_reactgood
            std_off=ch.std_offset_for_device_start
            mean=ch.avg_offset_for_device_start
            if reactgood!=None and not numpy.isnan(std_off):
                correlation.append(reactgood)
                std_offset.append(std_off)
                mean_offset.append(mean)
    plt.figure()
    plt.plot(correlation, std_offset, 'b*')
    plt.figure()
    plt.plot(correlation, mean_offset, 'r*')
    plt.figure()
    plt.plot(mean_offset, std_offset, 'g*')
    plt.show()         


#*************************************
#prepare data and plots for webapp
  

def colorCode(r):
    #input=codes for reaction to the device
    #output=color for the plot for each reaction
    if r=='None':
        c='y'
    elif r=='ReactGood':
        c='g'
    elif r=='ReactBad':
        c='b'
    elif r=='Disruption':
        c='r'
    elif r=='PreUseDisruption':
        c='k'
    else:
        c='c'
    #print c
    return c


##*************plotting
def plotUserReport(crec, cresult, f_name): 
    #matplotlib plotting, figures saved in files
    X_use=cresult.use
    Y_use=[0]*len(X_use)
    X_nt=[]
    Y_nt=[]
    for el in cresult.dist_event:
        X_nt.append(el[0])
        Y_nt.append(el[1])
    X=dict()
    Y=dict()
    for r in rsp:
        X[r]=[]
        Y[r]=[]
    for r in cresult.dev.keys():
        for el in cresult.dev[r]:
            X[r].append(el[0])
            Y[r].append(el[1])
        k=rsp.index(r)
        Y[r]=numpy.array(Y[r])+(k+1)*0.1
    #if i>100 and i<120:
    fig=plt.figure()
    ax=fig.add_subplot(111)
    #ax.text(180,3.55, 'id '+str(cid))
    ax.text(150, 4.25, 'history: %d, %d' % (
        crec.initial.month_count, 
        crec.initial.night_count) )
    ax.set_title('user'+str(cresult.child_id))
    ax.set_xlabel('days')
    ax.set_ylabel('events')
    ax.plot(X_use, Y_use, 'cs')
    if len(Y_nt)>0:
        ax.plot(X_nt, Y_nt, 'k*')
    for r in X.keys():
        #print colorCode(r)
        if len(X[r])>0:
            ax.plot(X[r], Y[r], colorCode(r)+'o')
    ax.set_ylim((-1,5))
    ax.set_xlim((0,200))
    fig.savefig(f_name)
    plt.close(fig)



#generate webapp data
def prepJsPoints(x, y):
    xypairs = numpy.transpose(numpy.vstack((x, y)))
    good = xypairs[~numpy.isnan(xypairs).any(1), :]
    return tuple(map(tuple, good))

def makeRateChartData(crec, cresult):
    #for each user, collect data in dict that will be put into database that will be used for plots on webapp
    if len(cresult.use)<20:
        return {}
    init=crec.initial
    init_rate=init.month_count*init.night_count/30.0
    init_day=min(cresult.use) - 1
    return dict(
        init_rate=init_rate,
        init_day=init_day,
        rate_disr=prepJsPoints(cresult.use,
                               cresult.rate_dst_event),
        rate_worked=prepJsPoints(cresult.use,
                                 cresult.complete_dev_rate['ReactGood']),
        rate_use=prepJsPoints(cresult.use,
                              cresult.dev_use_rate),
    )

def nan_to_str(fmt, val):
    #works with measures of "success" to the treatment so that can be properly displayed on the webapp i.e., displaying invalid values such as None or numpy nan
    try:
        if val is None or numpy.isnan(val):
            return None
        return fmt % val
    except:
        print repr(val)
        raise
        
    return float(val)

def makeSummaryData(cresult):
    #generates summary outcomes (various ways to measure how well the treatment worked) to be displayed on the webapp
    return [
        ('correlation G', 
         nan_to_str("%.3f", cresult.corr_dev_reactgood)),
        ('correlation G (SPR)', 
         nan_to_str("%.3f", cresult.corr_spr_reactgood)),
        ('correlation B', 
         nan_to_str("%.3f", cresult.corr_dev_dev_total)),
        ('correlation B (SPR)', 
         nan_to_str("%.3f", cresult.corr_spr_dev_total)),
        ('biggest improvement', nan_to_str("%.2f", cresult.biggest_improvement[0])),
        ('best rate', nan_to_str("%.2f", cresult.biggest_improvement[2]))
    ]

def mkTable():
    #generates tables that contain values used in plots on the webapp
    con=psycopg2.connect(dbname='testdb')
    cur=con.cursor()
    cur.execute('DROP TABLE IF EXISTS report_data')
    cur.execute('CREATE TABLE report_data ('
                '  child_id INT PRIMARY KEY, '
                '  description VARCHAR(255), '
                '  user_report VARCHAR(255), '
                '  sleep_chart VARCHAR(255), '
                '  rate_chart_data JSON, '
                '  summary_data JSON'
                ')')
    return con

SELECTED_IDS = [151, 174, 207, 238, 258, 281, 306, 331, 360, 390]
assert len(SELECTED_IDS) == 10

def plotData(analysis):
    #generates data for ploting on the webapp
    i=0
    con=mkTable()
    cur=con.cursor()
    for cid, cresult in sorted(analysis.items()):
        crec = children[cid]

        user_report='figures2/user%04d.jpg' % cid
        #plotUserReport(crec, cresult, user_report)

        sleep_chart = 'figures2/schart%04d.jpg' % cid
        #plotSleepChart(crec, cresult, sleep_chart)
        #if cid == 63:
        #    plotSleepChart(crec, cresult, None)
        #    fail

        if cid not in SELECTED_IDS: continue

        rc_data=makeRateChartData(crec, cresult)
        if not rc_data:
            continue


        i+=1
        if (i % 10) == 0: print >>sys.stderr, cid,

        sum_data=makeSummaryData(cresult)

        try:
            cur.execute(
                'INSERT INTO report_data(child_id, description, user_report, '
                'sleep_chart, rate_chart_data, summary_data) VALUES (%s, %s, %s, %s, %s, %s)',
                (cid, 'child'+str(cid), user_report, sleep_chart,
                 json.dumps(rc_data), json.dumps(sum_data)))
        except:
            yaml.dump_all([crec, cresult], sys.stdout)
            raise

    assert i >= 10
    print >>sys.stderr
    con.commit()


#testing
children, analysis, start, end=analysis.processData()
#plotAllSleepCharts(children, analysis)
#plotUseRateVsBestRate(analysis)
#plotUseRateVsImprovRate(analysis)
#plotDevRateVsImprov(analysis)
#plotCorrVsOffset(analysis)
plotData(analysis)
