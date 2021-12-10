import numpy as np
import pandas as pd
import scipy.signal
import obspy
from obspy.taup import TauPyModel
from obspy import UTCDateTime
import os


class Event():
    def __init__(self,event_path,oned_model):#event path will specify directory corresponding to that event's data,
        '''
        __init__ function for Event class. Sets initial attributes.
        
        Parameters
        ----------
        event_path: string
            Path to subdirectory with event info and data.
        oned_model: string or dictionary
            1-D model to be used to determine predicted travel times. Either a string that specifies a known model
            or a dictionary with keys "depths", "P velocities", and "S velocities", the values will be lists, the "depths" list must
            start with zero). Each element of the "velocities" list will be a list of two elements: P wave velocity then S wave
            velocity. For the customizable model to work, the velocities all the way down to the center of the Earth must be specified.
            Example: {"depths": [0,20,20,35,...], "P velocities": [5.8,5.8,6.5,6.5,...], "S velocities": [3.46,3.46,3.85,3.85,...]}
        '''
        
        event_df = pd.read_csv(event_path + "/evinfo",sep="\s",engine='python',header=None)
        if np.any(event_df.iloc[0,:]=='''"'''):
            col_to_destroy=np.argwhere(event_df.iloc[0,:]=='''"''')
            event_df=event_df.drop(columns=col_to_destroy[0][0])
        event_df.index=event_df[0]
        self.event_df=event_df
        self.station_list = list(event_df[0])
        self.depth = event_df.iloc[0,9]
        if isinstance(self.depth, str):
            self.depth = np.float64(self.depth.strip('''" '''))
        self.longitude = event_df.iloc[0,8]
        if isinstance(self.longitude, str):
            self.longitude = np.float64(self.longitude.strip('''" '''))
        self.latitude = event_df.iloc[0,7]
        if isinstance(self.latitude, str):
            self.latitude = np.float64(self.latitude.strip('''" '''))
        self.year = event_df.iloc[0,1]
        if isinstance(self.year, str):
            self.year = np.int(self.year.strip('''" '''))
        self.month = event_df.iloc[0,2]
        if isinstance(self.month, str):
            self.month = np.int(self.month.strip('''" '''))
        self.day = event_df.iloc[0,3]
        if isinstance(self.day, str):
            self.day = np.int(self.day.strip('''" '''))
        self.hour = event_df.iloc[0,4]
        if isinstance(self.hour, str):
            self.hour = np.int(self.hour.strip('''" '''))
        self.minute = event_df.iloc[0,5]
        if isinstance(self.minute, str):
            self.minute = np.int(self.minute.strip('''" '''))
        self.second = event_df.iloc[0,6]
        if isinstance(self.second, str):
            self.second = np.float64(self.second.strip('''" '''))
        self.event_time = UTCDateTime(str(self.year) + "-" +  str(self.month) + "-" + str(self.day) + "T" +
                                      str(self.hour) + ":" + str(self.minute) + ":" + str(self.second))
        self.magnitude = np.float64(event_df.iloc[0,10])
        self.angle_to_stations = event_df.iloc[:,11] #pd dataframe
        self.distance_to_stations = self.angle_to_stations * (np.pi/180) * 6,378.1370 #pd dataframe
        self.oned_model = oned_model
        self.event_path = event_path
        
    
    def identify_estimated_arrival(self,sequence,sampling_rate,phase,pred_time): 
        '''
        Method that takes in the seismogram data of one component and outputs the index of the array corresponding to 
        the arrival of the P phase.
        
        Parameters
        ----------
        sequence: ndarray 
            Array of numbers representing displacements of seismogram over time.
        sampling_rate: int
            Number of samples taken per second
        phase: str
            Phase to be identified in the seismogram
        pred_time: float
            The number of seconds after the start time of the file that the oned_model predicts the phase will arrive.
            
        Returns
        -------
        arrival_index: int
            Index of inputted array corresponding to arrival of desired phase (either P or S).
        '''
        if phase == 'P':
            absolute_max = np.max(np.absolute(sequence[int(np.rint((pred_time-10)*sampling_rate)):int(np.rint((pred_time+10)*sampling_rate))]))
            arrival_index = np.argwhere(np.absolute(sequence)==absolute_max)[0][0]
            
        if phase == 'S':
            absolute_max = np.max(np.absolute(sequence[int(np.rint((pred_time-20)*sampling_rate)):int(np.rint((pred_time+20)*sampling_rate))]))
            arrival_index = np.argwhere(np.absolute(sequence)==absolute_max)[0][0]
        
        return arrival_index
    
    
    
    def get_observed_arrival_times_P(self,estimate_arrival=None):
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
        estimate_arrival: function
            User can input a function to estimate the observed arrival times instead of using the algorithm implemented below.
            This 'estimated_arrival' function would take in a list of vertical component stack files each corresponding to a
            station that detected this event. Should output a list of the P observed arrival times. If None, then the method below is used.
        '''
        
        if estimate_arrival is None:
        
            earliest_arrival = None
            master_station = None
            master_index = None
            master_sampling_rate = None
            z_stacks=pd.DataFrame(dtype=float,columns=self.station_list)
            start_time_list = []
            pred_times = self.get_predicted_arrival_times_P()

            for index,station in enumerate(self.station_list):
                stack = obspy.read(self.event_path + '/' + station + '/STACK_Z')
                vertical_data = stack[0].data
                z_stacks.loc[:,station] = vertical_data
                start_time = stack[0].stats.starttime
                start_time_list.append(start_time)
                sampling_rate = stack[0].stats.sampling_rate
                pred_time = np.float64(pred_times.loc[station]) - (start_time - self.event_time)
                estimated_arrival_time = start_time + (self.identify_estimated_arrival(vertical_data,int(sampling_rate),'P',pred_time)/sampling_rate)
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

            if np.any(np.isnan(z_stacks)):
                raise Exception("stacks are not the same length")

            for index,station in enumerate(self.station_list):
                if index != master_index:
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
            
        else:
            stacks=[]
            for index,station in enumerate(self.station_list):
                stack = obspy.read(self.event_path + '/' + station + '/STACK_Z')
                stacks.append(stack)
            arrival_time_list = estimate_arrival(stacks)
            self.observed_arrival_times_P = pd.DataFrame(data = {'Times':arrival_time_list}, index = self.station_list)
            
        
        
        return self.observed_arrival_times_P #this will be a dataframe with the observed P arrival times per station that can be indexed with the station code or by the regular index
    
    def get_observed_arrival_times_S(self,estimate_arrival=None):
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
        estimate_arrival: function
            User can input a function to estimate the observed arrival times instead of using the algorithm implemented below.
            This 'estimated_arrival' function would take in two lists, representing the two horizontal components, of stack files
            each corresponding to a station that detected this event. Should output a list of the S observed arrival times. If
            None, then the method below is used.
        '''
        
        if estimate_arrival is None:
        
            earliest_arrival = None
            master_station = None
            master_index = None
            master_sampling_rate = None
            r_stacks=pd.DataFrame(dtype=float,columns=self.station_list)
            t_stacks=pd.DataFrame(dtype=float,columns=self.station_list)
            p_stacks=pd.DataFrame(dtype=float,columns=self.station_list)
            start_time_list = []
            pred_times = self.get_predicted_arrival_times_S()

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
                #if radial_variance > transverse_variance:
                    #clearest_data = radial_data
                #else:
                    #clearest_data = transverse_data
                clearest_data = radial_data
                p_stacks.loc[:,station] = clearest_data
                pred_time = np.float64(pred_times.loc[station]) - (start_time - self.event_time)
                estimated_arrival_time =  start_time + (self.identify_estimated_arrival(clearest_data,int(sampling_rate),'S',pred_time)/sampling_rate)
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

            if np.any(np.isnan(r_stacks)):
                raise Exception("stacks are not the same length")

            for index,station in enumerate(self.station_list):
                if index != master_index:

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
            
        else:
            stacks1=[]
            stacks2=[]
            for index,station in enumerate(self.station_list):
                stack1 = obspy.read(self.event_path + '/' + station + '/STACK_R')
                stack2 = obspy.read(self.event_path + '/' + station + '/STACK_T')
                stacks1.append(stack1)
                stacks2.append(stack2)
            arrival_time_list = estimate_arrival(stacks1,stacks2)
        
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
            with open('oned_model.nd', 'w') as f:
                lines = []
                for index, depth in enumerate(self.oned_model.get("depths")):
                    lines.append(str(depth) + " " + str(self.oned_model.get("P velocities")[index]) + " " + str(self.oned_model.get("S velocities")[index]) + "\n")
                f.writelines(lines)
            TauPyModel.taup_create.build_taup_model('oned_model.nd',output_folder='.')
            model = TauPyModel(model='oned_model.npz')
    
        arrival_time_list = []
    
        for station in self.station_list:
            arrival = model.get_travel_times(source_depth_in_km=self.depth,
                                              distance_in_degree=self.angle_to_stations[station], phase_list=["P","S"])
            
       
            arrival_time_list.append(arrival[0].time)
            
        self.predicted_arrival_times_P = pd.DataFrame(data = {'Times':arrival_time_list}, index = self.station_list)
        
        return self.predicted_arrival_times_P
        
    def get_predicted_arrival_times_S(self):
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
            with open('oned_model.nd', 'w') as f:
                lines = []
                for index, depth in enumerate(self.oned_model.get("depths")):
                    lines.append(str(depth) + " " + str(self.oned_model.get("P velocities")[index]) + " " + str(self.oned_model.get("S velocities")[index]) + "\n")
                f.writelines(lines)
            TauPyModel.taup_create.build_taup_model('oned_model.nd',output_folder='.')
            model = TauPyModel(model='oned_model.npz')
            
        arrival_time_list = []
    
        for station in self.station_list:
            arrival = model.get_travel_times(source_depth_in_km=self.depth,
                                              distance_in_degree=self.angle_to_stations[station], phase_list=["S"])
            
            arrival_time_list.append(arrival[0].time)
            
        self.predicted_arrival_times_S = pd.DataFrame(data = {'Times':arrival_time_list}, index = self.station_list)
        
        return self.predicted_arrival_times_S
    
    def get_residuals(self):
        '''
        Method that returns the residuals for for the P and S wave arrival times for the stations corresponding to the
        events.
            
        Returns
        -------
        p_residuals: DataFrame
            DataFrame with that is the difference in seconds between the observed arrival of the P phase at different stations for
            different events and the predicted arrival times.
        s_residuals: DataFrame
            DataFrame with that is the difference in seconds between the observed arrival of the S phase at different stations for
            different events and the predicted arrival times.
        '''
        obs_p = self.get_observed_arrival_times_P()
        obs_s = self.get_observed_arrival_times_S()
        pred_p = self.get_predicted_arrival_times_P()
        pred_s = self.get_predicted_arrival_times_S()
        
        p_residuals = obs_p - pred_p
        s_residuals = obs_s - pred_s
        
        return p_residuals,s_residuals
    
    def get_stacks(self,station):
        '''
        Method that returns the a list of the sac files corresponding to the specified station for this event.
            
        Returns
        -------
        stack_list: list of sac objects
            Stack files corresponding to the specified station for this event in the component order z,r,t
        '''
        z=obspy.read(self.event_path + '/' + station + '/STACK_Z')
        r=obspy.read(self.event_path + '/' + station + '/STACK_R')
        t=obspy.read(self.event_path + '/' + station + '/STACK_T')
        stack_list=[z,r,t]
        return stack_list
    
        
    
    
    
def read_in_events(events_directory_path, oned_model):
    '''
    This function takes the path to the directory where the event subdirectories are stored and outputs the corresponding Event
    objects.
        
    Parameters
    ----------
    events_directory_path: str
        String indicating the path to the directory where the event subdirectories are stored.
    oned_model: String that represents the name of an established 1-D velocity model (i.e. "AK135") or a dictionary with keys   
        "interface depth" and "velocities", the values will be lists; the length of "velocities" list will be one more than the
        length of "interface depth" list.

            
    Returns
    -------
    list_of_events: list of Event objects
        List of event objects that correspond to the code names inputted into the function.
    '''
    
    event_codes=[]
    for code in next(os.walk(events_directory_path))[1]:
        if code[:3] == '20':
            event_codes.append(code)
    
    
    list_of_events = []
    for event_name in event_codes:
        new_event = Event(event_name,oned_model)
        list_of_events.append(new_event)
        
    return list_of_events


