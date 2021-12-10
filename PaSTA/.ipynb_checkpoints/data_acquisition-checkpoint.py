import obspy
from obspy import read
from obspy import UTCDateTime, read
from obspy.clients.fdsn import Client
import numpy as np
import obspy.geodetics.base
from obspy.geodetics.base import gps2dist_azimuth
import obspy.taup
from obspy.taup import TauPyModel
import matplotlib.pyplot as plt
import os
import shutil
import csv



def SNR(trZ,trR,thresZ,thresR,fmin,fmax):
    '''
    Function that takes in stack files (SAC format) in the vertical and radial components, applies a bandpass filter, and evaluates if
    their signal to noise ratios (SNR) meet a threshold. Returns a 1 if they do, a 0 if they don't.

    Parameters
    ----------
    trZ: Obspy Trace object
        Seismogram data stack corresponding to the z (vertical) component.
    trR: Obspy Trace object
        Seismogram data stack corresponding to the r (horizontal) component.
    thresZ: int or float
        Minimum value that the signal to noise ratio of the z stack can be and still have the function return a 1.
    thresR: int or float
        Minimum value that the signal to noise ratio of the r stack can be and still have the function return a 1.
    fmin: int or float
        Minimum cutoff frequency used when applying the bandpass filter before calculating the SNR.
    fmax: int or float
        Miaximum cutoff frequency used when applying the bandpass filter before calculating the SNR.
    
    Returns
    -------
    snr_indicator: int
        Value of 1 if SNR meets inputted thresholds, value of 0 otherwise.
    '''
    trZ_f = trZ.copy()
    trR_f = trR.copy()
    trZ_f.filter(type="bandpass",freqmin=fmin,freqmax=fmax,zerophase=True)
    trR_f.filter(type="bandpass",freqmin=fmin,freqmax=fmax,zerophase=True)

    ref = trZ_f.stats.starttime
    Nstart = ref + 7.5
    Nend = ref + 25
    Pstart = ref + 30
    Pend = ref + 37.5
    Psstart = ref + 45
    Psend = ref +52.5

    trZ_noise = trZ_f.copy()
    trZ_signalP = trZ_f.copy()
    trZ_signalPs = trZ_f.copy()
    trR_noise = trR_f.copy()
    trR_signalP = trR_f.copy()
    trR_signalPs = trR_f.copy()

    trZ_noise.trim(starttime=Nstart,endtime=Nend)
    trZ_signalP.trim(starttime=Pstart,endtime=Pend)
    trZ_signalPs.trim(starttime=Psstart,endtime=Psend)
    trR_noise.trim(starttime=Nstart,endtime=Nend)
    trR_signalP.trim(starttime=Pstart,endtime=Pend)
    trR_signalPs.trim(starttime=Psstart,endtime=Psend)

    nZ = trZ_noise.max()
    pZ = trZ_signalP.max()
    psZ = trZ_signalPs.max()
    SNR_Z = 10*np.log10((pZ/nZ)**2)
    SNS_Z = 10*np.log10((pZ/psZ)**2)

    nR = trR_noise.max()
    pR = trR_signalP.max()
    psR = trR_signalPs.max()
    SNR_R = 10*np.log10((pR/nR)**2)
    SNS_R = 10*np.log10((pR/psR)**2)

    if (SNR_Z>=thresZ)&(SNR_R>=thresR):
        return 1
    else:
        return 0
    
    
    
def acquire_data(start_date, end_date, min_magnitude, min_distance, max_distance, max_event_depth, oned_model, before, after, latitude_min, latitude_max, longitude_min, longitude_max, save_directory, cat_client = "IRIS", wf_client = "IRIS",network = "XP",sta = "*",location = "*",cha = "BH?"):
    '''
    Function that locates seismogram stations within a certain region, obtains data from them from events that happened within certain
    specifications, and then creates directories to store that data in. Each event is assigned a code with its corresponding
    subdirectory having that name. All these subdirectories go into the specified save directory. Each event subdirectory has a csv
    file called 'evinfo' which contains a table where each row corresponds to a station and the columns are station name, event year,
    month, day, hour, minute, second, event latitude, longitude, depth (in kilometers), magnitude, great-circle distance from event to
    station in degrees, event back-azimuth (the direction of the event with respect to the station). Also in these subdirectories are 
    additional subdirectories corresponding to each station that had data for that event. These contain the SAC files with the data for 
    that station for that event (named STACK_R, STACK_T, STACK_Z for the radial, transverse, and vertical components, respectively).

    Parameters
    ----------
    start_date: str
        A string in the format of "year-month-day" or "0000-00-00". This specifies the beginning of the range of time from which we
        want data/seismograms.
    end_date: str
        A string in the format of "year-month-day" or "0000-00-00". This specifies the end of the range of time from which we
        want data/seismograms.
    min_magnitude: int or float
        The minimum magnitude of events we want to look at.
    min_distance: int or float
        The minimum distance in degrees an event must be to the specified latitude and longitude for its data to be collected.
    max_distance: int or float
        The maximum distance in degrees an event must be to the specified latitude and longitude for its data to be collected.
    max_event_depth: int or float
        The maximum depth in kilometers of the events being collected by this function.
    oned_model: str
        The 1-D model used to estimate where the start time for each data should be. The start time for each SAC file will be 'before'
        seconds before the expected P wave arrival estimated by this station.
    before: int or float
        The number of seconds before the estimated P wave arrival for an event which defines the beginning of the time window for the
        stack files.
    after: int or float
        The number of seconds after the estimated P wave arrival for an event which defines the end of the time window for the stack
        files. It is recommended that this be at least several hundred seconds if the S wave arrival is to be captured.
    latitude_min: float
        The minimum above which the latitude of all the stations collected must be.
    latitude_max: float
        The maximum below which the latitude of all the stations collected must be.
    longitude_min: float
        The minimum above which the longitude of all the stations collected must be.
    longitude_max: float
        The maximum below which the longitude of all the stations collected must be.
    save_directory: str
        The directory to which all the data directories and files will be saved to.
    cat_client: str
        The source catalog for the events.
    wf_client: str
        The source catalog for the waveforms, likely the same as the cat_client.
    network: str
        The code for the stations whose data we want.
    sta: str
        The specific code for the stations, use "*" to leave this value open.
    location: str
        Location code for the desired network of stations, use "*" to leave this value open.
    cha: str
        The code for the desired channels.

    '''
    
    
    #####GET STATIONS###################
    model = TauPyModel(model=oned_model)
    client = Client(wf_client)

    print("** Getting station data from "+wf_client+" **")


    starttime = UTCDateTime(start_date)
    endtime = UTCDateTime(end_date)


    print(cha)

    inv = client.get_stations(network = network, station = sta, level="response",  channel = cha,
          minlongitude=longitude_min, maxlongitude=longitude_max, minlatitude=latitude_min, maxlatitude=latitude_max, 
          starttime=starttime, endtime=endtime)

    print(starttime,endtime)
    stime = starttime
    etime = endtime

    station_list = []
    lat_list = []
    lon_list = []
    no = 0
    for net in inv:
        for station in net:
            if ((station.latitude > latitude_min) & (station.latitude < latitude_max)):
                if ((station.longitude > longitude_min) & (station.longitude < longitude_max)):
                    no += 1
                    station_info = [no,f'{station.latitude:7.4f}',f'{station.longitude:8.4f}',station.elevation/1000.0,station.code]
                    station_list.append(station_info)
                    lat_list.append(station.latitude)
                    lon_list.append(station.longitude)
                    print(station_info)                

    with open(save_directory + 'all_sta_cord', "w") as output:
        writer = csv.writer(output, delimiter=' ', lineterminator='\n')
        writer.writerows(station_list)

    
    
    
    
    
    
    #####GET EVENTS###################
    slat = np.mean(lat_list)
    slon = np.mean(lon_list)
    radm = 6378137.0 #Radius of the Earth in meters
    pie = np.pi
    radian = radm*pie/180.0

    evclient = Client(cat_client)
    print("** Getting event data from "+cat_client+" **")
    cat = evclient.get_events(starttime=stime,
          endtime=etime,minmagnitude=min_magnitude, maxdepth=max_event_depth,latitude=slat,longitude=slon,
          minradius=min_distance,maxradius=max_distance,orderby='time-asc')
    print(cat)
    
    
    
    
    
    #####GET WAVEFORMS###################
    ev_list = []
    count = 0
    print("Grab some popcorn, this is going to take a while!")
    for event in cat[0:-1]:
        count += 1
        t = UTCDateTime(event.origins[0].time)
        ev_name = str(t.year) + str(t.julday).zfill(3) + str(t.hour).zfill(2) + str(t.minute).zfill(2)
        #print("Event No." + str(count) + " - " + ev_name)
        ev_dir = save_directory + ev_name + '/'
        if os.path.exists(ev_dir):
            shutil.rmtree(ev_dir)
        os.makedirs(ev_dir)
        ev_lat = event.origins[0].latitude
        ev_lon = event.origins[0].longitude
        ev_dep = event.origins[0].depth/1000.
        mag = event.magnitudes[0].mag

        rawdata = 0
        rows = []
        for net in inv:
            for sta in net:
                slat = sta.latitude
                slon = sta.longitude
                aaa = gps2dist_azimuth(ev_lat,ev_lon,slat,slon)
                epi = aaa[0]/radian
                baz = aaa[2]

                arrP = model.get_travel_times(source_depth_in_km=ev_dep, \
                                           distance_in_degree=epi,phase_list=["P"])
                ptime = arrP[0].time
                if sta.end_date is None:
                    sta.end_date = UTCDateTime(2099, 1, 1, 0, 0, 0)
                if ((sta.start_date+24*3600 < (t+ptime-before)) & (sta.end_date-24*3600 > (t+ptime+after))):
                    try:
                        st = client.get_waveforms(network=net.code, station=sta.code,
                                                  location=location, channel=cha,
                                                  starttime=t+ptime-before, endtime=t+ptime+after,attach_response=True)
                    except:
                        #print(sta.code + ' somehow does not have data for event ' + ev_name)
                        continue
                    try:
                        st.rotate(method="->ZNE", inventory=inv)
                        st.rotate(method="NE->RT",back_azimuth = baz)
                    except:
                        #print('Channels of ' + sta.code + ' cannot be rotated, probably due to unaligned time series.')
                        continue
                    rawdata += 1
                    st.detrend('demean')
                    st.detrend('linear')
                    st.filter(type="bandpass",freqmin=0.03,freqmax=4.9,zerophase=True)
                    st.taper(0.05,type='Hamming',max_length=7.5)
                    try:
                        st[2].write(ev_dir+'dum'+st[2].stats.channel[2]+'.sac',format='SAC')
                        st[1].write(ev_dir+'dum'+st[1].stats.channel[2]+'.sac',format='SAC')
                        st[0].write(ev_dir+'dum'+st[0].stats.channel[2]+'.sac',format='SAC')
                    except:
                        #print(sta.code + ' does not have 3 required channels for ' + ev_name)
                        continue
                    number_of_points = (before + after)*40
                    if (st[2].stats.npts < number_of_points) | (st[1].stats.npts < number_of_points) | (st[0].stats.npts < number_of_points):
                        #print(sta.code + ' does not have complete record for ' + ev_name)
                        continue

                    stZ = read(ev_dir+'dumZ.sac')
                    stR = read(ev_dir+'dumR.sac')
                    stT = read(ev_dir+'dumT.sac')
                    os.remove(ev_dir+'dumZ.sac')
                    os.remove(ev_dir+'dumR.sac')
                    os.remove(ev_dir+'dumT.sac')
                    trZ = stZ[0]
                    trR = stR[0]
                    trT = stT[0]
                    f1 = SNR(trZ,trR,5,4,0.03,0.3)
                    f2 = SNR(trZ,trR,5,4,0.03,0.6)
                    f3 = SNR(trZ,trR,5,4,0.03,1.0)
                    trZ.stats.sac.user1=f1
                    trZ.stats.sac.user2=f2
                    trZ.stats.sac.user3=f3
                    trR.stats.sac.user1=f1
                    trR.stats.sac.user2=f2
                    trR.stats.sac.user3=f3
                    trT.stats.sac.user1=f1
                    trT.stats.sac.user2=f2
                    trT.stats.sac.user3=f3
                    trZ.stats.sac.user0=arrP[0].ray_param/6371
                    trR.stats.sac.user0=arrP[0].ray_param/6371
                    trT.stats.sac.user0=arrP[0].ray_param/6371

                    if f1+f2+f3>0:
                        sta_dir = ev_dir + sta.code +'/'
                        os.makedirs(sta_dir)
                        trZ.write(sta_dir + 'STACK_Z',format='SAC')
                        trR.write(sta_dir + 'STACK_R',format='SAC')
                        trT.write(sta_dir + 'STACK_T',format='SAC')
                        row = [sta.code,str(t.year),str(t.month).zfill(2),str(t.day).zfill(2),
                               str(t.hour).zfill(2),str(t.minute).zfill(2),
                               str(t.second).zfill(2)+'.'+str(int(t.microsecond/100)).zfill(3),
                              f'{ev_lat:7.4f}',f'{ev_lon:8.4f}',ev_dep,mag,f'{epi:7.4f}',f'{baz:7.4f}']
                        rows.append(row)

        with open(ev_dir + 'evinfo', "w") as output:
            writer = csv.writer(output, delimiter=' ', lineterminator='\n')
            writer.writerows(rows)

        qcdata = len(rows)
        #print('Number of stations: ' + str(rawdata) + ' Passed QC: ' + str(qcdata))
        ev_sta = [ev_name,str(rawdata),str(qcdata)]
        ev_list.append(ev_sta)

        if qcdata < 8:
            shutil.rmtree(ev_dir)

    with open(save_directory + 'ev_list', "w") as output:
        writer = csv.writer(output, delimiter=' ', lineterminator='\n')
        writer.writerows(ev_list)

