from obspy import read, Stream, UTCDateTime
import numpy as np
import numpy.ma as ma
import sys
from preprocessing import obtain_timeshifts, correct_timeshifts

try:
    # filenames
    ff = sys.argv[1]
    stn = sys.argv[2]
    fs = sys.argv[3]
except:
    print("USAGE: python seg2mseed.py <fileformat> <stationlist> <filelist>")

# directory containing files to convert
fd = "/scratch/flindner/G7Jul01/%s2/" % ff

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
if ff == "SEG":
    # assign station and channel number
    for tr in range(ntr):
        st[tr].stats.station = stn[tr]
        st[tr].stats.channel = chn[tr]

# read following files, assign station and channel (if format=SEG2), and add to first stream
print("read files ...")
for f in range(1, nfls):
    new = read(fd + fls[f])
    if ff == "SEG":
        for tr in range(ntr):
            new[tr].stats.station = stn[tr]
            new[tr].stats.channel = chn[tr]
    st += new
# sort for station and channel
st.sort(keys=["station", "channel"])
print("... done!")

# empty stream for time corrected streams
tcor_st = Stream()
# unique station-channel pair
stations = list(np.unique(stn))
channels = list(np.unique(chn))
# reverse ordering
channels = channels[::-1]
# select traces for one station-channel pair, shift in time, merge and add to container
run = 0
print("calculate time shifts and merge files ...")
for station in stations:
    for channel in channels:
        stream = st.select(station=station, channel=channel)    
        stream.sort(keys=["starttime"])
        # test continuity of stream
        test = stream.copy()
        test.merge(method=1)
        if not ma.is_masked(test[0].data):
            if run == 0:
                # copy stream
                dtst = stream.copy()
                # calculate time shifts
                dts, tp = obtain_timeshifts(dtst)
            print("... ", station, channel, " ...")
            # correct for time shifts and add merged trace to container
            cor_stream = correct_timeshifts(stream, dts, tp)
            tcor_st += cor_stream
            run += 1
        else:
            print("Data gaps detected")
print("... done!")

# get date and time from stream
time = tcor_st[0].stats.starttime + 10*60
y = time.year
m = time.month
d = time.day
h = time.hour
delta = tcor_st[0].stats.delta

# save traces as mseed files
print("write traces ...")
t1_trim = UTCDateTime(y, m, d, h, 0)
if h < 23:
    t2_trim = UTCDateTime(y, m, d, h + 1, 0)
elif h == 23:
    t2_trim = UTCDateTime(y, m, d + 1, 0, 0)
tcor_st.trim(t1_trim, t2_trim - delta)
for tr in tcor_st:
    tr.stats.network = "4D"

    fname = "/scratch/flindner/G7Jul01/MSEED/" \
            + "4D." + tr.stats.station + ".." + tr.stats.channel + "." \
            + "%04d%02d%02d%02d" % (y, m, d, h) + "_from%s2" % ff
    tr.write(fname, format="MSEED")
    print(tr)
