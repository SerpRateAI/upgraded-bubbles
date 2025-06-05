import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import glob
#import datetime
#from scipy import signal
from obspy import read, UTCDateTime
from obspy.signal.cross_correlation import correlation_detector

#from multiprocessing import Pool
import sys
day = str(sys.argv[1])
#year_month = '2020-01'

#path = '/home/karina/New_bubbles/'
path = '/fp/projects01/ec332/data/'
   
do_plot = False
if do_plot:
#    savefigto = '/Figures/'
    savefigto = '/fp/homes01/u01/ec-karinal/New_bubbles/Figures/'
    lts = 10
    fts = 15
    #date_form = mdates.DateFormatter("%y-%m-%d %H:%M:%S.%f")
    date_form = mdates.DateFormatter("%H:%M:%S.%f")
 
# Template 
templ_label = ['Tremor1', 'Tremor2', 'Tremor3', 'Tremor4', 'ChirpA1', 'ChirpA2', 'ChirpB', 'Earthquake']
templ1 = [UTCDateTime('2019-07-18 19:49:45.4'), UTCDateTime('2019-06-02 01:18:07'), UTCDateTime('2019-08-13 04:59:07.5'), UTCDateTime('2019-12-02 05:33:20.3'),  UTCDateTime('2019-10-10 08:06:06'), UTCDateTime('2019-10-10 08:06:12'), UTCDateTime('2019-09-16 02:04:33'), UTCDateTime('2019-11-12 03:29:06.5')]
templ2 = [UTCDateTime('2019-07-18 19:49:57'), UTCDateTime('2019-06-02 01:18:20'), UTCDateTime('2019-08-13 04:59:13'), UTCDateTime('2019-12-02 05:33:27'), UTCDateTime('2019-10-10 08:06:10.5'), UTCDateTime('2019-10-10 08:06:15.5'), UTCDateTime('2019-09-16 02:04:35'), UTCDateTime('2019-11-12 03:29:13')]

templ_label = ['Tremor1', 'Tremor2', 'Tremor3', 'Tremor4', 'ChirpA1', 'ChirpA2', 'ChirpB', 'Earthquake']
templ1 = [UTCDateTime('2019-07-18 19:49:45.4'), UTCDateTime('2019-06-02 01:18:07'), UTCDateTime('2019-08-13 04:59:07.5'), UTCDateTime('2019-12-02 05:33:20.3'), UTCDateTime('2019-10-10 08:06:06'), UTCDateTime('2019-10-10 08:06:12'), UTCDateTime('2019-09-16 02:04:33'), UTCDateTime('2019-11-12 03:29:06.5')]
templ2 = [UTCDateTime('2019-07-18 19:49:49'), UTCDateTime('2019-06-02 01:18:11.5'), UTCDateTime('2019-08-13 04:59:11'), UTCDateTime('2019-12-02 05:33:22.5'), UTCDateTime('2019-10-10 08:06:10.5'), UTCDateTime('2019-10-10 08:06:15.5'), UTCDateTime('2019-09-16 02:04:35'), UTCDateTime('2019-11-12 03:29:10.5')]

# read and process function
def read_pros(filename, flip=False):
    tr = read(filename)
    trStat = tr[0].stats
    statName = ('.').join((trStat.network,trStat.station,trStat.location,trStat.channel))
    dateName = str(trStat.starttime).split('T')[0]
    tr.merge(fill_value='latest')
    if flip: tr[0].data = -1 * tr[0].data
    tr.detrend('demean')
    #tr.filter('highpass', freq=15, corners=4, zerophase=True)
    tr.filter('bandpass', freqmin=10,freqmax=50, corners=4, zerophase=True)

    return tr, statName, dateName

def convert_trace_to_pascals(trace):
    """
    converts hydrophone data to pascals
    % convert digital counts to Pa
    % Q330S+ has 419430 counts/volt
    % hydrophone sensitivity is -165 dB, re: 1V/uPa
    % or p in uPa = V/10^-16.5
    """
    trace.data = trace.data/419430 # converts to volts
    trace.data = trace.data/(10**(-165/20.)) # converts to microPascals
    trace.data = trace.data/1e6 # converts to pascals
    # return trace
    return trace.data

def convert_stream_to_pascals(stream):
    """
    wrapper function to convert an entire stream to pascals
    """
    for n, tr in enumerate(stream):
        # print(type(tr))
        stream[n].data = convert_trace_to_pascals(tr)
        
    return stream

def make_template(bubbles, templ_label, channels, do_plot=False):
    no_stats = len(channels)
    templates = []
    for (j, bub1), bub2, tempL in zip(enumerate(bubbles[0]), bubbles[1], templ_label):
    #    print('template ' + str(j))
        # First station
        day = str(bub1).split('T')[0]
        fileName = path + 'hydrophones/' + day + '/7F.'+channels[0]+'.GDH.mseed'
        tr, statName, dateName = read_pros(fileName, flip=False)
        
        template = tr.select(id=statName).slice(bub1,bub2)
        
        # Save and plot templates from day 1
        if do_plot:
            fig, axs = plt.subplots(no_stats,figsize=(20, 20))
            ax = axs[0]
            ax.grid(which='both', axis='both', linestyle='--', linewidth=0.4) 
            ax.plot(template[0].times("matplotlib"),template[0].data)
            #axs.set_xlabel('Velocity [m/s]', fontsize=fts)
            ax.set_title(tempL + ': ' + dateName + '\n' + statName, fontsize=fts)
            ax.xaxis.set_major_formatter(date_form)
            ax.tick_params(labelsize=lts, rotation=10)
        
        for i, cl in enumerate(channels[1:]):
#            if j==int(len(bubbles[0])-1): # Adjust the last template to remove other bubbles
#                bub1-= 0.02
#                bub2-=0.02
            
            if cl=='B00.04': flip=True
            else: flip=False
            
            fileName =  path + 'hydrophones/' + day + '/7F.'+cl+'.GDH.mseed'
            
            tr, statName, dateName = read_pros(fileName, flip=flip)
            templateP = tr.select(id=statName).slice(bub1,bub2)
            template+= templateP
            
            if do_plot:
                ax = axs[i+1]
                ax.grid(which='both', axis='both', linestyle='--', linewidth=0.4) 
                ax.plot(templateP[0].times("matplotlib"),templateP[0].data)
                #axs.set_xlabel('Velocity [m/s]', fontsize=fts)
                ax.set_title(statName, fontsize=fts)
                ax.xaxis.set_major_formatter(date_form)
                ax.tick_params(labelsize=lts, rotation=10)
            
            del templateP
            
        if do_plot:
            plt.tight_layout()
            plt.savefig(savefigto+str(j)+tempL+"_templates.png", orientation='landscape')
        
        templates.append(template)
        del template
        
    return templates

def match_day(day, templates, channels):
    no_stats = len(channels)
    # Match with the templates
    fileName = path + 'hydrophones/' + day + '/7F.'+channels[0]+'.GDH.mseed'
    tr, statName, dateName = read_pros(fileName, flip=False)
    
    data = tr
    for cl in channels:
        fileName = path + 'hydrophones/' + day + '/7F.'+cl+'.GDH.mseed'
        
        if cl=='B00.04': flip=True
        else: flip=False
            
        tr, statName, dateName = read_pros(fileName, flip=flip)
        data+= tr
        
        height = 0.6
        distance = 1.1
    data.merge(fill_value='latest')
    
    print(dateName)
    detections, sims = correlation_detector(
                                    stream=data,
                                    templates=templates,
                                    heights=height,
                                    distance=distance,
                                    plot=None
                                    # , similarity_func=tm.simf
                                    )   
    del data, sims
    
    
#    df = pd.DataFrame(detections)
#    print('detected {n} events on day {y}.{d}'.format(n=df.shape[0], d=dateName))
##     print('elapsed time:', datetime.now()-start)
##    print('elapsed time:',useful_variables.elapsed_time())
#    df.to_csv(data_writing_location + '.csv', index=False)
    return pd.DataFrame(detections)

#channels = ['A00.03','A00.04','A00.05','A00.06','B00.01', 'B00.02','B00.03','B00.04','B00.05','B00.06']
print('number of templates: {n}'.format(n=len(templ_label)))

channels = ['A00.03','A00.04','A00.05','A00.06']
templates = make_template([templ1,templ2], templ_label, channels, do_plot=do_plot)
dfA = match_day(day, templates, channels)
dfA['Channel'] = 'A'

channels = ['B00.01', 'B00.02','B00.03','B00.04','B00.05','B00.06']
templates = make_template([templ1,templ2], templ_label, channels, do_plot=do_plot)
dfB = match_day(day, templates, channels)
dfB['Channel'] = 'B'

df = pd.concat([dfA, dfB], ignore_index=True)
del dfA, dfB

if len(df)==0:
    print('No match')
else:
    print(np.unique(df.template_id, return_counts=True))
    temp_typ = np.full_like(df.template_id, 'nan', dtype='U10')
    for i, tl in enumerate(templ_label):
        temp_typ[df.template_id==i] = tl 
    print(np.unique(temp_typ, return_counts=True))
    df['template_type'] = temp_typ

    data_writing_location = path + 'event_detections/' + day
    df.to_csv(data_writing_location + '.csv', index=False)