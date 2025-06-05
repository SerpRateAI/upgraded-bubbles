import numpy as np
import os
import pickle
import glob
from obspy import read, UTCDateTime
from matplotlib import dates as mdates
import matplotlib.pyplot as plt
plt.rcParams["date.converter"] = "concise"

from scatseisnet import ScatteringNetwork

path = '/fp/projects01/ec332/data/'
dirpath_save = path + "Scattering_networks/"
#dirpath_save = "seisscatt_test/"

# Create directory to save the results
#os.makedirs(dirpath_save, exist_ok=True)

do_plot = False
if do_plot:
    savefigto = '/fp/homes01/u01/ec-karinal/New_bubbles/seisscatt_test/'
    lts = 10
    fts = 15
    date_form = mdates.DateFormatter("%H:%M:%S.%f")

# Scattering network Parameters 
# Select parameters
segment_duration_seconds = 20.0
sampling_rate_hertz = 200.0
samples_per_segment = int(segment_duration_seconds * sampling_rate_hertz)
bank_keyword_arguments = (
    {"octaves": 6, "resolution": 4, "quality": 2},
    {"octaves": 6, "resolution": 1, "quality": 2},
)

# Create scattering network
network = ScatteringNetwork(
    *bank_keyword_arguments,
    bins=samples_per_segment,
    sampling_rate=sampling_rate_hertz,
)

print(network)

# Extract segment length (from any layer)
segment_duration_seconds = network.bins / network.sampling_rate
overlap = 0.5

# Save the scattering network with Pickle
filepath_save = os.path.join(dirpath_save, "scattering_network.pickle")
with open(filepath_save, "wb") as file_save:
    pickle.dump(network, file_save, protocol=pickle.HIGHEST_PROTOCOL)

if do_plot:
    # Visualize the filter banks
    # Loop over network layers
    for bank in network.banks:
    
        # Create axes (left for temporal, right for spectral domain)
        fig, ax = plt.subplots(1, 2, sharey=True)
    
        # Show each wavelet
        for wavelet, spectrum, ratio in zip(
            bank.wavelets, bank.spectra, bank.ratios
        ):
    
            # Time domain
            ax[0].plot(bank.times, wavelet.real + ratio, "C0")
    
            # Spectral domain (log of amplitude)
            ax[1].plot(bank.frequencies, np.log(np.abs(spectrum) + 1) + ratio, "C0")
    
        # Limit view to three times the temporal width of largest wavelet
        width_max = 3 * bank.widths.max()
    
        # Labels
        ax[0].set_ylabel("Octaves (base 2 log)")
        ax[0].set_xlabel("Time (seconds)")
        ax[0].set_xlim(-width_max, width_max)
        ax[0].grid()
        ax[1].set_xscale("log")
        ax[1].set_xlabel("Frequency (Hz)")
        ax[1].grid()

#%% Read, change endtime, preprocess and 
days_all = np.sort(glob.glob(path + 'hydrophones/*/'))

channels = ['A00.03','A00.04','A00.05','A00.06']
#channels = ['B00.01', 'B00.02','B00.03','B00.04','B00.05','B00.06']

# Start with day 1
fileName = days_all[0] + '/7F.'+ channels[0] +'.GDH.mseed'
stream = read(fileName)
for cl in channels[1:]:
        fileName = days_all[0] + '/7F.'+ cl +'.GDH.mseed'
        # Read each day's data
        stream += read(fileName)
trStat = stream[0].stats
print(trStat.starttime)

# Loop over every 3 days
for idx in np.arange(1,len(days_all),3):
#for idx in np.arange(1,6,3):
    #idx = np.where(days_all==path+day+'/')[0][0]
    print(idx)
    days = days_all[idx:idx+2]
    # Read each day's data
    for day in days:
        #fileName = '/fp/projects01/ec332/data/hydrophones/' + days[0] + '/7F.B00.01.GDH.mseed'
        for cl in channels:
            fileName = day + '/7F.'+ cl +'.GDH.mseed'
            # Read each day's data
            stream += read(fileName)
    
    #print('Raw')
    #print(stream)
    #Combine
    stream.sort()       # Sort traces by id and start time
    stream.merge(method=1, fill_value='interpolate')

    #Slice
    #print('Old endtime: ' + str(stream[0].stats.endtime))
    new_entime = UTCDateTime(str(stream[0].stats.endtime).split('T')[0]+'T00:00:00.1')
    streamNew = stream.copy()
    streamNew = streamNew.slice(endtime=new_entime)
    #print('New endtime: ' + str(new_entime))


    # NEED TO ADD A CHECK HERE THAT THE NEW STREAMS ARE THE SAME LENGTH!!!!! (SEE LAST SLURM)
    
    #print('New')
    print(streamNew)
    
    stream = stream.slice(starttime=new_entime)
    #print('Next')
    #print(stream)
    trStat = stream[0].stats
    print(trStat.starttime)
    
    if do_plot:
        print('Raw')
        streamNew.plot(rasterized=True)

    # Process
    sampling_rate0 = streamNew[0].stats.sampling_rate
    streamNew.decimate(factor=int(sampling_rate0/sampling_rate_hertz), strict_length=False, no_filter=True)
    
    streamNew.detrend('demean')
    #stream.filter('highpass', freq=10, corners=4, zerophase=True)
    streamNew.filter('bandpass', freqmin=10,freqmax=50, corners=4, zerophase=True)

    print('Processed')
    #print(streamNew)
    if do_plot:
        streamNew.plot(rasterized=True)

    # Scattering network
    # Gather list for timestamps and segments
    timestamps = list()
    segments = list()
    # Collect data and timestamps
    for traces in streamNew.slide(segment_duration_seconds, segment_duration_seconds * overlap):
        #print(traces)
        timestamps.append(mdates.num2date(traces[0].times(type="matplotlib")[0]))
        segments.append(np.array([trace.data[:-1] for trace in traces]))

    print('Segmented')
    ## Scattering transformation
    scattering_coefficients = network.transform(segments, reduce_type=np.median)
    
    # Save
    np.savez(dirpath_save+"scattering_coefficientsA_day"+str(idx)+".npz",
    order_1=scattering_coefficients[0],
    order_2=scattering_coefficients[1],
    times=timestamps)
             
    if do_plot:
        ##Observe result from a single channel
        # Extract the first channel
        channel_id = 0
        trace = streamNew[channel_id]
        order_1 = np.log10(scattering_coefficients[0][:, channel_id, :].squeeze())
        center_frequencies = network.banks[0].centers
        
        # Create figure and axes
        fig, ax = plt.subplots(2, sharex=True, dpi=300)
        # Plot the waveform
        ax[0].plot(trace.times("matplotlib"), trace.data, rasterized=True, lw=0.5)
        # First-order scattering coefficients
        ax[1].pcolormesh(timestamps, center_frequencies, order_1.T, rasterized=True)
        # Axes labels
        ax[1].set_yscale("log")
        ax[0].set_ylabel("Counts")
        ax[1].set_ylabel("Frequency (Hz)")
        # Show
        plt.tight_layout()
        plt.savefig(savefigto+"scattering_coefficients0.png", orientation='landscape')

    del streamNew, traces
    
