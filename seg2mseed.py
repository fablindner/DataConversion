from obspy import read, Stream, UTCDateTime
import numpy as np
import sys
from preprocessing import correct_timeshift

# file directory
fd = "/media/TOSHIBA EXT/flindner/G7Jul01/SEG2/"

# filenames
stn = sys.argv[1]
fs = sys.argv[2]
# read station and channel names
stn, chn = np.loadtxt(stn, skiprows=1, usecols=(1,2), unpack=True, dtype=bytes).astype(str)
# number of traces
ntr = len(chn)
# files to process
fls = np.loadtxt(fs, dtype=bytes).astype(str)
# number of files
nfls = len(fls)

# read first file
st = read(fd + fls[0])
if len(st) != ntr:
    raise ValueError("different number of traces in %s and %s" % (st, stn))
# assign station and channel number
for tr in range(ntr):
    st[tr].stats.station = stn[tr]
    st[tr].stats.channel = chn[tr]

# read following files, assign station and channel, and add to first stream
for f in range(1, nfls):
    print(fls[f])
    new = read(fd + fls[f])
    for tr in range(ntr):
        new[tr].stats.station = stn[tr]
        new[tr].stats.channel = chn[tr]
    st += new
# sort for station and channel
st.sort(keys=["station", "channel"])

# empty stream for time corrected streams
tcor_st = Stream()
# unique station-channel pair
stations = list(np.unique(stn))
channels = list(np.unique(chn))
channels = channels[::-1]
# select traces for one station-channel pair, shift in time, merge and add to container
for station in stations:
    for channel in channels:
        stream = st.select(station=station, channel=channel)    
        stream.sort(keys=["starttime"])
        cor_stream = correct_timeshift(stream)
        tcor_st += cor_stream

print(tcor_st.__str__(extended=True))
# get date and time from stream
time = tcor_st[0].stats.starttime + 10*60
y = time.year
m = time.month
d = time.day
h = time.hour

t1_trim = UTCDateTime(y, m, d, h, 0)
t2_trim = UTCDateTime(y, m, d, h + 1, 0)
tcor_st.trim(t1_trim, t2_trim)
for tr in tcor_st:
    tr.stats.network = "4D"
    fname = "/media/TOSHIBA EXT/flindner/G7Jul01/MSEED/" \
            + "4D." + tr.stats.station + ".." + tr.stats.channel + "." \
            + "%04d%02d%02d%02d" % (y, m, d, h) + "_fromSEG2"
    tr.write(fname, format="MSEED")