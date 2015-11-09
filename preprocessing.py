from obspy import read, Stream, Trace
import numpy as np


def correct_timeshift(st):
    '''function to merge overlapping time series associated with time shifts'''
    '''st: obspy stream containing one station channel pair             '''

    # check if all traces are from same station/channel
    for tr in st:
        if tr.stats.station == st[0].stats.station:
            pass
        else:
            raise NameError("Traces in st recorded by different stations")
        if tr.stats.channel == st[0].stats.channel:
            pass
        else:
            raise NameError("Traces in st recorded on different channels")
    
    # empty stream for time corrected traces; add first trace
    cor = Stream() + st[0]
    st = st[1:]
    
    # add next trace
    for tr in st:
        cor += tr
        # helper stream - overlap
        helper = cor.copy()
        t1 = helper[1].stats.starttime
        t2 = helper[0].stats.endtime
        helper.trim(t1, t2)
        delta = helper[0].stats.delta
        npts = helper[0].stats.npts
        # even number of samples
        if (npts % 2) != 0: 
            helper.trim(t1, t2-delta)
            delta = helper[0].stats.delta
        
        # cross-correlation in order to obtain timeshift
        xcorr = np.correlate(helper[0].data, helper[1].data, "same")
        xcorr = Trace(data=xcorr)
        xcorr.stats.delta = delta
        xcorr.resample(1000)
        ndelta = xcorr.stats.delta
        nnpts = xcorr.stats.npts
        dt = (np.argmax(xcorr.data) - nnpts/2.) * ndelta
        dt = round(dt, 3)
    
        # add estimated time shift to second trace
        cor[1].stats.starttime += dt
    
        # test if time series are perfectly aligned
        check = True
        count = 1
        while check:
            start = cor[1].stats.starttime
            end = cor[0].stats.endtime
            test = cor.slice(start, end)
            if test[0].stats.starttime == test[1].stats.starttime:
                eq_st = True
            eq_ar = np.array_equal(test[0].data, test[1].data)
            if eq_ar == True and eq_st == True:
                check = False  
            # if not perfectly aligned, shift by one sample and try again
            else:
                if count == 1:
                    cor[1].stats.starttime -= 100*delta
                elif count > 1:
                    cor[1].stats.starttime += delta
            count += 1
        # if aligned, merge!
        cor.merge(method=1)

    return cor