import numpy as np
import pandas as pd
import glob


detected_events = pd.read_csv('events7greater.csv')

event_typ = np.full_like(detected_events.starttime, 'No class', dtype='U10') 
plotfolder = '/home/karina/New_bubbles/geq7_plots/'
types = ['Tremor', 'Earthquake', 'Bubble', 'Other', 'Chirps']

for typ in types:
    
    events = np.sort(glob.glob(plotfolder + '/' + typ + '/*'))
    
    for ev in events:
        evSt = ev.split('/')[-1][:-4]
        
        event_typ[detected_events.starttime==evSt] = typ

#%%
print(np.unique(event_typ, return_counts=True))

detected_events['Type'] = event_typ

detected_events.to_csv('events7greater_wclass.csv')

#%% Template list
temps = np.sort(glob.glob(plotfolder + '/' + 'Template' + '/*'))
    
for ev in temps:
    evSt = ev.split('/')[-1].split('_')[0]
    tp = ev.split('/')[-1].split('_')[1]
    
    print(tp)
    print(detected_events.starttime[detected_events.starttime==evSt])
    print(detected_events.endtime[detected_events.starttime==evSt])
    
    
    
#%%
bubs = pd.read_csv('/home/karina/New_bubbles/detections/all.csv')
print(np.unique(bubs.template_id, return_counts=True))

temp_typ = np.full_like(bubs.template_id, 'nan', dtype='U10')

templ_label = ['Tremor1', 'Tremor2', 'Tremor3', 'Tremor4', 'ChirpA1', 'ChirpA2', 'ChirpB', 'Earthquake']
for i, tl in enumerate(templ_label):
    
    temp_typ[bubs.template_id==i] = tl 

print(np.unique(temp_typ, return_counts=True))
bubs['template_type'] = temp_typ