import numpy as np
import pandas as pd
import os
import obspy
from obspy.taup import TauPyModel
from obspy import UTCDateTime

class Event():
    def __init__(self,event_path,oned_model):#event path will specify directory corresponding to that event's data,
        #oned model will be a string wiht the model name OR a dictionary with keys "interface depth" and "velocities", the values will be lists; the length of "velocities" list will be one more than the length of "interface depth" list
        '''
        __init__ function for Event class. Sets initial attributes.
        
        Parameters
        ----------
        event_path: string
            Path to subdirectory with event info and data.
        oned_model: string or dictionary
            1-D model to be used to determine predicted travel times. Either a string that specifies a known model
            or a dictionary with keys "interface depth" and "velocities", the values will be lists; the length of
            "velocities" list will be one more than the length of "interface depth" list.
        '''
        
        event_df = pd.read_csv(event_path + "/evinfo",sep="\s",engine='python',header=None)
        event_df.index=event_df[0]
        self.event_df=event_df
        self.station_list = list(event_df[0])
        self.depth = event_df.iloc[0,9]
        self.longitude = event_df.iloc[0,8]
        self.latitude = event_df.iloc[0,7]
        self.year = event_df.iloc[0,1]
        self.month = event_df.iloc[0,2]
        self.day = event_df.iloc[0,3]
        self.hour = event_df.iloc[0,4]
        self.minute = event_df.iloc[0,5]
        self.second = event_df.iloc[0,6]
        self.event_time = UTCDateTime(str(self.year) + "-" +  str(self.month) + "-" + str(self.day) + "T" +
                                      str(self.hour) + ":"+ str(self.minute) + ":" + str(self.second))
        self.magnitude = event_df.iloc[0,10]
        self.angle_to_stations = event_df.iloc[:,11] #pd dataframe
        self.distance_to_stations = event_df.iloc[:,11] * (np.pi/180) * 6,378.1370 #pd dataframe
        self.oned_model = oned_model
        self.event_path = event_path
        
    def identify_estimated_arrival(self,sequence): #Would be cool to implement a more sophisticated algorithm here!
        '''
        Method that takes in the seismogram data of one component and outputs the index of the array corresponding to 
        the arrival of either the P phase (if it is the veritcal component) or the S phase (if it is the radial or
        traverse component).
        
        Parameters
        ----------
        sequence: ndarray
            Array of numbers representing displacements of seismogram over time.
            
        Returns
        -------
        arrival_index: int
            Index of inputted array corresponding to arrival of desired phase (either P or S).
        '''
        
        absolute_max = np.max(np.absolute(sequence))
        arrival_index = np.argwhere(np.absolute(sequence)==absolute_max)[0][0]
        
        return arrival_index
    
    def get_observed_arrival_times_P(self):
        '''
        Method that returns the time difference between the time of the event and the time of the inital P arrival at
        each station. This method determines the station with the earliest arrival in the vertical component, makes
        that the "master" station, then cross-correlates with the data from the vertical components of the other
        stations to find their respective arrival times for the P phase. Stores final dataframe in attribute as well as
        returns it.
            
        Returns
        -------
        observed_arrival_times_P: DataFrame
            DataFrame with station codes and their observed P phase arrival times.
        '''
        
        earliest_arrival = None
        master_station = None
        master_index = None
        master_sampling_rate = None
        z_stacks=pd.DataFrame(dtype=float,columns=self.station_list)
        start_time_list = []
        
        for index,station in enumerate(self.station_list):
            stack = obspy.read(self.event_path + '/' + station + '/STACK_Z')
            vertical_data = stack[0].data
            z_stacks.loc[:,station] = vertical_data
            start_time = stack[0].stats.starttime
            start_time_list.append(start_time)
            sampling_rate = stack[0].stats.sampling_rate
            estimated_arrival_time = (self.identify_estimated_arrival(vertical_data)/sampling_rate) + start_time
            if earliest_arrival is None:
                earliest_arrival = estimated_arrival_time
                master_station = station
                master_index = index
                master_sampling_rate = sampling_rate
            elif estimated_arrival_time < earliest_arrival:
                earliest_arrival = estimated_arrival_time
                master_station = station
                master_index = index
                master_sampling_rate = sampling_rate
                
        arrival_time_array = np.zeros(len(self.station_list))
        arrival_time_array[master_index] = earliest_arrival - self.event_time
                
        if np.any(np.isnan(test_df.loc[:,master_station])):
            raise Exception("stacks are not the same length")
        
        for index,station in enumerate(self.station_list):
            if index != master_index:
                if np.any(np.isnan(test_df.loc[:,station])):
                    raise Exception("stacks are not the same length")
                
                master_data = z_stacks.loc[:,master_station]
                current_data = z_stacks.loc[:,station]
                station_corr = scipy.signal.correlate(master_data,current_data)
                max_corr_ind = np.argwhere(np.max(station_corr)==station_corr)[0][0]
                center = master_data.size-1
                data_offset = (center - max_corr_ind)/master_sampling_rate
                start_time_offset = start_time_list[index] - start_time_list[master_index]
                offset = data_offset + start_time_offset
                arrival_time_array[index] = earliest_arrival - self.event_time + offset
                
        arrival_time_list = list(arrival_time_array)
        self.observed_arrival_times_P = pd.DataFrame(data = {'Times':arrival_time_list}, index = self.station_list)
        
        return self.observed_arrival_times_P #this will be a dataframe with the observed P arrival times per station that can be indexed with the station code or by the regular index
    
    def get_observed_arrival_times_S():
        '''
        Method that returns the time difference between the time of the event and the time of the inital S arrival at
        each station. This method determines the station with the earliest arrival in the vertical component, makes
        that the "master" station, then cross-correlates with the data from the vertical components of the other
        stations to find their respective arrival times for the S phase. Stores final dataframe in attribute as well as
        returns it.
            
        Returns
        -------
        observed_arrival_times_S: DataFrame
            DataFrame with station codes and their observed S phase arrival times.
        '''
        earliest_arrival = None
        master_station = None
        master_index = None
        master_sampling_rate = None
        r_stacks=pd.DataFrame(dtype=float,columns=self.station_list)
        t_stacks=pd.DataFrame(dtype=float,columns=self.station_list)
        p_stacks=pd.DataFrame(dtype=float,columns=self.station_list)
        start_time_list = []
        
        for index,station in enumerate(self.station_list):
            stack_r = obspy.read(self.event_path + '/' + station + '/STACK_R')
            stack_t = obspy.read(self.event_path + '/' + station + '/STACK_T')
            radial_data = stack_r[0].data
            transverse_data = stack_r[0].data
            r_stacks.loc[:,station] = radial_data
            t_stacks.loc[:,station] = transverse_data
            if stack_r[0].stats.starttime != stack_t[0].stats.starttime:
                raise Exception("horizontal components have different start times")
            start_time = stack_r[0].stats.starttime
            start_time_list.append(start_time)
            if stack_r[0].stats.sampling_rate != stack_t[0].stats.sampling_rate:
                raise Exception("horizontal components have sampling rates")
            sampling_rate = stack_r[0].stats.sampling_rate
            radial_variance = np.var(radial_data)
            transverse_variance = np.var(transverse_data)
            if radial_variance > transverse_variance:
                clearest_data = radial_data
            else:
                clearest_data = transverse_data
            p_stacks.loc[:,station] = clearest_data
            estimated_arrival_time = (self.identify_estimated_arrival(clearest_data)/sampling_rate) + start_time
            if earliest_arrival is None:
                earliest_arrival = estimated_arrival_time
                master_station = station
                master_index = index
                master_sampling_rate = sampling_rate
            elif estimated_arrival_time < earliest_arrival:
                earliest_arrival = estimated_arrival_time
                master_station = station
                master_index = index
                master_sampling_rate = sampling_rate
                
        arrival_time_array = np.zeros(len(self.station_list))
        arrival_time_array[master_index] = earliest_arrival - self.event_time
                
        if np.any(np.isnan(test_df.loc[:,master_station])):
            raise Exception("stacks are not the same length")
        
        for index,station in enumerate(self.station_list):
            if index != master_index:
                if np.any(np.isnan(test_df.loc[:,station])):
                    raise Exception("stacks are not the same length")
                
                master_data = p_stacks.loc[:,master_station]
                current_data = p_stacks.loc[:,station]
                station_corr = scipy.signal.correlate(master_data,current_data)
                max_corr_ind = np.argwhere(np.max(station_corr)==station_corr)[0][0]
                center = master_data.size-1
                data_offset = (center - max_corr_ind)/master_sampling_rate
                start_time_offset = start_time_list[index] - start_time_list[master_index]
                offset = data_offset + start_time_offset
                arrival_time_array[index] = earliest_arrival - self.event_time + offset
                
        arrival_time_list = list(arrival_time_array)
        self.observed_arrival_times_S = pd.DataFrame(data = {'Times':arrival_time_list}, index = self.station_list)
        
        return self.observed_arrival_times_S #this will be a dataframe with the observed S arrival times per station that can be indexed with the station code or by the regular index
    
    
    def get_predicted_arrival_times_P(self):
        '''
        Method that returns the predicted arrival times for the P phase for an event based on the 1-D model specified
        for that Event object. Stores final dataframe in attribute as well as returns it.
            
        Returns
        -------
        predicted_arrival_times_P: DataFrame
            DataFrame with station codes and their predicted P phase arrival times.
        '''
        if isinstance(self.oned_model, str):
            model = TauPyModel(model=self.oned_model)
        else:
            #Figure out how to construct custom model form Obspy documentation!!!!!
            #my_model = ...
            #model = TauPyModel(model=my_model)
    
        arrival_time_list = []
    
        for station in self.station_list:
            arrival = model.get_travel_times(source_depth_in_km=self.depth,
                                              distance_in_degree=self.angle_to_stations[station], phase_list=["P"])
            arrival_time_list.append(arrival[0].time)
            
        self.predicted_arrival_times_P = pd.DataFrame(data = {'Times':arrival_time_list}, index = self.station_list)
        
        return self.predicted_arrival_times_P
        
    def predicted_arrival_times_S(self):
        '''
        Method that returns the predicted arrival times for the S phase for an event based on the 1-D model specified
        for that Event object. Stores final dataframe in attribute as well as returns it.
            
        Returns
        -------
        predicted_arrival_times_S: DataFrame
            DataFrame with station codes and their predicted S phase arrival times.
        '''
        if isinstance(self.oned_model, str):
            model = TauPyModel(model=self.oned_model)
        else:
            #Figure out how to construct custom model form Obspy documentation!!!!!
            #my_model = ...
            #model = TauPyModel(model=my_model)
    
        arrival_time_list = []
    
        for station in self.station_list:
            arrival = model.get_travel_times(source_depth_in_km=self.depth,
                                              distance_in_degree=self.angle_to_stations[station], phase_list=["S"])
            arrival_time_list.append(arrival[0].time)
            
        self.predicted_arrival_times_S = pd.DataFrame(data = {'Times':arrival_time_list}, index = self.station_list)
        
        return self.predicted_arrival_times_S
    
    def get_residuals(self):
        '''
        Calculate and output residuals for P and S.
        '''
        
    
    
def read_in_events(event_list, oned_model):#Reads in event path list and outputs a list of Event objects
    '''
    This function takes a list of event names (which should represent directories) and outputs the corresponding Event objects.
        
    Parameters
    ----------
    event_list: list of Strings
        List of the code names representing each event we are dealing with.
    oned_model: String that represents the name of an established 1-D velocity model (i.e. "AK135") or a dictionary with keys   
        "interface depth" and "velocities", the values will be lists; the length of "velocities" list will be one more than the
        length of "interface depth" list.

            
    Returns
    -------
    list_of_events: list of Event objects
        List of event objects that correspond to the code names inputted into the function.
    '''
    list_of_events = []
    for event_name in event_list:
        new_event = Event(event_name,oned_model)
        list_of_events.append(new_event)
        
    return list_of_events
    
    