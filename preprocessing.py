from obspy import read, Stream, Trace
import numpy as np
import warnings

def obtain_timeshifts(st):
    '''function to compute time shift of overlapping files              '''
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
    cor = Stream()
    # "take_previous" determines if first trace is subject to syncronization to previous mseed file
    take_previous = False
    # check if previous mseed file is available
    try:
        stn = st[0].stats.station
        chn = st[0].stats.channel
        stime = st[0].stats.starttime
        prev = read("/scratch/flindner/G7Jul01/MSEED/4D.%s..%s.%04d%02d%02d%02d*" \
                % (stn, chn, stime.year, stime.month, stime.day, stime.hour), starttime=stime - 5)
        cor += prev[0]
        # available -> synchronize
        take_previous = True
    except:
        pass
    # in case of no sync to previous file, start with timeshift between first and second trace
    if take_previous == False:
        cor += st[0]
        st = st[1:]
        warnings.warn("Synchronization to previous MSEED file failed!")
    # empty list for dts
    dts = []
    
    # add next trace
    for i, tr in enumerate(st):
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
        count = 0
        while check:
            start = cor[1].stats.starttime
            end = cor[0].stats.endtime
            test = cor.slice(start, end)
            eq_st = False
            if round(test[0].stats.starttime.timestamp, 4) == round(test[1].stats.starttime.timestamp, 4):
                eq_st = True
            eq_ar = np.array_equal(test[0].data, test[1].data)
            if eq_ar == True and eq_st == True:
                # if perfectly aligned, quit testing loop
                check = False  
                if count == 0:
                    dts.append(dt) 
            # if not perfectly aligned, undo time shift obtained by xcorr and shift manually
            else:
                if count == 0:
                    # undo time shift obtained from xcorr
                    cor[1].stats.starttime -= dt
                    mshift = int((end - start) / delta) + 1000
                    cor[1].stats.starttime -= mshift*delta
                elif count > 0:
                    cor[1].stats.starttime += delta
                cor[1].stats.starttime.microsecond = np.round(cor[1].stats.starttime.microsecond, 3)
            count += 1
        # in case of take_previous == True, remove previous file from stream (was just needed for time sync)
        if i==0 and take_previous:
            for tr in cor.select(network="4D"):
                cor.remove(tr)
        # time shift obtained from shifting
        if count > 1:
            dt = (-mshift + count - 2) * delta
            dts.append(dt) 

        cor.merge(method=1)
    # check number of written times shifts
    dts = np.asarray(dts)
    if len(dts) != len(st):
        raise ValueError("Number of time shifts is not correct")
    return dts, take_previous




def correct_timeshifts(st, dts, take_previous):
    '''function to merge overlapping time series associated with time shifts'''
    '''st: obspy stream containing one station channel pair                 '''
    '''dts: array containing dts                                            '''
    '''take_previous: use file of previous hour to synchronize              '''

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
    
    # check if nuber of dts is correct
    if take_previous:
        length = len(st)
    else:
        length = len(st) - 1
    if len(dts) != length:
        raise ValueError("number of dts is not correct")

    # empty stream for time corrected traces; add first trace
    cor = Stream() + st[0]
    st = st[1:]
    if take_previous:
        cor[0].stats.starttime += dts[0]
        dts = dts[1:]
    
    # add next trace, apply timeshift, check alignment, possibly correct time shift and merge
    for it, tr in enumerate(st):
        cor += tr
        dt = dts[it]
        dt = round(dt, 3)
        dt = "%.3f" % dt
        dt = float(dt)
        cor[1].stats.starttime += dt
        delta = cor[0].stats.delta

        count = 0
        check = True
        while check:
            # test stream - overlap
            start = cor[1].stats.starttime
            end = cor[0].stats.endtime
            test = cor.slice(start, end)
    
            # test if time series are perfectly aligned
            eq_st = False
            if round(test[0].stats.starttime.timestamp, 4) == round(test[1].stats.starttime.timestamp, 4):
                eq_st = True
            eq_ar = np.array_equal(test[0].data, test[1].data)
            if eq_ar == True and eq_st == True:
                # if aligned, merge!
                cor.merge(method=1)
                check = False
            else:
                if count == 0:
                    cor[1].stats.starttime += delta
                elif count == 1:
                    cor[1].stats.starttime -= 2*delta
                else:
                    raise ValueError("Overlapping samples do not match: %s - %s" % (start, end))
            count += 1
    return cor
