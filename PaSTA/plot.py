import numpy as np
import pandas as pd
import os
import obspy
import matplotlib.pyplot as plt
import scipy.odr

    

def fit_line(x,y):
    '''
    This function applies total least squares regression to the inputted x and y data. Returns list with the slope and the
    intercept of said line.
        
    Parameters
    ----------
    x: ndarray
        x data to fit to a line
    y: ndarray
        y data to fit to a line
            
    Returns
    -------
    beta: ndarray
        An array of length two where the first term is the slope of the fitted line and the second term is the y-intercept.
    '''
    def f(B, x):
        return B[0]*x + B[1]

    linear = scipy.odr.Model(f)
    sx=np.std(x)
    sy=np.std(y)
    my_data = scipy.odr.Data(x,y)
    #my_data = scipy.odr.Data(x, y, wd=1./(sx**2), we=1./(sy**2))
    my_odr = scipy.odr.ODR(my_data, linear, beta0=[1., 0])
    my_output = my_odr.run()
    
    return my_output.beta
    

def histogram_residuals(events,phase,limits=None,**kwargs):
    '''
    This function creates a histogram of the time residuals for the inputted list of event objects for the phase specified (either 
    P or S).
        
    Parameters
    ----------
    events: list of Event objects
        A list of Event objects whose time residuals will be extracted and plotted.
    phase: str
        String specifying the desired phase to be plotted (either 'S' or 'P').
    limits: ndarray, optional
        Numpy array of length two that specifies the minimum and maximum range of residuals that should be plotted. 
            
    Returns
    -------
    fig: Figure
        The figure object corresponding to the histogram.
    ax: Axes
        The axes object corresponding to the histogram.
    '''

    residuals=np.asarray([])
    for event in events:
        p_resid,s_resid = event.get_residuals()
        if phase=='P':
            residuals=np.concatenate((residuals,np.squeeze(p_resid)))
        elif phase=='S':
            residuals=np.concatenate((residuals,np.squeeze(s_resid)))
        else:
            raise Exception("Invalid phase entered")
    
    
    
    fig,ax = plt.subplots(figsize=(13,13));
    
    plt.rcParams.update({
        "font.size": 25})
    if limits is not None:
        mask=np.all([residuals < limits[1], residuals > limits[0]],axis=0)
        ax.hist(residuals[mask],**kwargs);
    else:
        ax.hist(residuals,**kwargs);
    if phase=='P':
        ax.set_title("P wave travel time residuals");
    else:
        ax.set_title("S wave travel time residuals");
        
    return fig,ax
    
    
    
def plot_residuals(events,limits=None,line=False,**kwargs):
    '''
    This function creates a plot of the S wave time residuals versus corresponding the P wave time residuals. If the line argument
    is set to True, it will also fit and plot a line to these points with orthogonal distance regression.
        
    Parameters
    ----------
    events: list of Event objects
        A list of Event objects whose time residuals will be extracted and plotted.
    limits: ndarray, optional
        A 2x2 numpy array that specifies the min and max range for the P phase then the S phase for what residuals should be
        included in the fit (this is to exclude outliers that result from bad cross correlation of time pickings). First row is
        the limits for the P residuals, the second row is the limits for the S residuals.
    line: boolean, optional
        If this value is set to True, a line will be fitted and plotted over these points with orthogonal distance regression.
            
    Returns
    -------
    fig: Figure
        The figure object corresponding to the plot.
    ax: Axes
        The axes object corresponding to the plot.
    '''
    
    p_residuals=np.asarray([])
    s_residuals=np.asarray([])
    for event in events:
        p_resid,s_resid = event.get_residuals()
        p_residuals=np.concatenate((p_residuals,np.squeeze(p_resid)))
        s_residuals=np.concatenate((s_residuals,np.squeeze(s_resid)))
        
    
    fig,ax = plt.subplots(figsize=(13,13));
    
    plt.rcParams.update({
        "font.size": 25})
       
    ax.tick_params(length=15,width=1,labelsize=25,direction='inout');
    ax.set_xlabel("P wave travel time residuals");
    ax.set_ylabel("S wave travel time residuals");
    if limits is None:
        ax.scatter(p_residuals,s_residuals,**kwargs);
        if line:
            beta = fit_line(p_residuals,s_residuals)
            y = p_residuals*beta[0] + beta[1]
            ax.plot(p_residuals,y);
    else:
        mask1 = np.all([p_residuals < limits[0][1], p_residuals > limits[0][0]],axis=0)
        mask2 = np.all([s_residuals < limits[1][1], s_residuals > limits[1][0]],axis=0)
        ax.scatter(p_residuals[np.all([mask1,mask2],axis=0)],s_residuals[np.all([mask1,mask2],axis=0)],**kwargs);
        if line:
            beta = fit_line(p_residuals[np.all([mask1,mask2],axis=0)],s_residuals[np.all([mask1,mask2],axis=0)])
            y = p_residuals[np.all([mask1,mask2],axis=0)]*beta[0] + beta[1]
            ax.plot(p_residuals[np.all([mask1,mask2],axis=0)],y);
    
        

    
    return fig,ax



def fit_residuals(events,limits=None):
    '''
    This function applies orthogonal least squares regression to the corresponding P and S wave time residuals for the inputted
    list of Event objects. Returns the list of the slope and the intercept of said line.
        
    Parameters
    ----------
    events: list of Event objects
        A list of Event objects whose time residuals will be fitted to a line with orthogonal least squares regression.
    limits: ndarray, optional
        A 2x2 numpy array that specifies the min and max range for the P phase then the S phase for what residuals should be
        included in the fit (this is to exclude outliers that result from bad cross correlation of time pickings). First row is
        the limits for the P residuals, the second row is the limits for the S residuals.
            
    Returns
    -------
    beta: ndarray
        An array of length two where the first term is the slope of the fitted line and the second term is the y-intercept.
    '''

    residuals=np.asarray([])
    for event in events:
        p_resid,s_resid = event.get_residuals()
        p_residuals=np.concatenate((residuals,np.squeeze(p_resid)))
        s_residuals=np.concatenate((residuals,np.squeeze(s_resid)))
        
    if limits is None:
        beta = fit_line(p_residuals,s_residuals)
    else:
        mask1 = np.all([p_residuals < limits[0][1], p_residuals > limits[0][0]],axis=0)
        mask2 = np.all([s_residuals < limits[1][1], s_residuals > limits[1][0]],axis=0)
        beta = fit_line(p_residuals[np.all([mask1,mask2],axis=0)],s_residuals[np.all([mask1,mask2],axis=0)])
    return beta

def plot_tx(events,phase,time_reduction_fraction=0,limits=None,**kwargs):
    '''
    This function creates a time curve for the arrival times of the specified phase (either P or S). If time_reduction_fraction is
        given a non-zero value, it will plot reduced time (T - X*time_reduction_fraction) vs distance.
        
    Parameters
    ----------
    events: list of Event objects
        A list of Event objects whose arrival times and distances will be plotted.
    phase: str
        The phase whose time curve is to be plotted. Either 'P' or 'S'.
    time reduction_fraction: float, optional
        The fraction of the distance to be subtracted from T to get the reduced time.
    limits: ndarray, optional
        Numpy array of length two that specifies the minimum and maximum range of times to be plotted.

    Returns
    -------
    fig: Figure
        The figure object corresponding to the plot.
    ax: Axes
        The axes object corresponding to the plot.
    '''
    fig,ax = plt.subplots(figsize=(13,13));

    plt.rcParams.update({"font.size": 25})

    ax.tick_params(length=15,width=1,labelsize=25,direction='inout');
    if time_reduction_fraction == 0:
        ax.set_ylabel("Time (s)");
    else:
        ax.set_ylabel("Reduced Time (s), T - X/%d" % (int(1/time_reduction_fraction)));
        ax.set_xlabel("Distance (km)");
        
        
        
    travel_times=np.asarray([])
    distances=np.asarray([])
    for event in events:
        if phase=='P':
            observed_arrivals = np.squeeze(event.get_observed_arrival_times_P())
            X = event.distance_to_stations[0]
            travel_times=np.concatenate((travel_times,observed_arrivals))
            distances=np.concatenate((distances,X))
        elif phase=='S':
            observed_arrivals = np.squeeze(event.get_observed_arrival_times_S())
            X = event.distance_to_stations[0]
            travel_times=np.concatenate((travel_times,observed_arrivals))
            distances=np.concatenate((distances,X))
        else:
            raise Exception("Invalid phase entered")
        

        
    
    if phase == 'P':
        ax.set_title("P wave time curve");
        if limits is None:
            ax.scatter(distances,travel_times - (distances*time_reduction_fraction),**kwargs);
        else:
            mask=np.all([travel_times < limits[1], travel_times > limits[0]],axis=0)
            ax.scatter(distances[mask],travel_times[mask] - (distances[mask]*time_reduction_fraction),**kwargs);
    elif phase == 'S':
        ax.set_title("S wave time curve");             
        if limits is None:
            ax.scatter(distances,travel_times - (distances*time_reduction_fraction),**kwargs);
        else:
            mask=np.all([travel_times < limits[1], travel_times > limits[0]],axis=0)
            ax.scatter(distances[mask],travel_times[mask] - (distances[mask]*time_reduction_fraction),**kwargs);
        
    return fig,ax


def plot_stacks(event,station):
    '''
    Method that plots the wave forms corresponding to the specified station for this event.
            
    Returns
    -------
    fig: Figure
        The figure object corresponding to the plots.
    (ax1,ax2,ax3): tuple of Axes objects
        Axes objects corresponding to the plots.
    '''
    fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(30,30));
    
    plt.rcParams.update({"font.size": 25})
    
    ax1.tick_params(length=15,width=1,labelsize=25,direction='inout');
    ax2.tick_params(length=15,width=1,labelsize=25,direction='inout');
    ax3.tick_params(length=15,width=1,labelsize=25,direction='inout');
    
    z=obspy.read(event.event_path + '/' + station + '/STACK_Z')
    r=obspy.read(event.event_path + '/' + station + '/STACK_R')
    t=obspy.read(event.event_path + '/' + station + '/STACK_T')
    
    ax1.plot(np.linspace(0,z[0].stats.npts/z[0].stats.sampling_rate,z[0].stats.npts),z[0].data);
    ax1.set_title("Vertical Component Waveform for Station {}".format(station));
    ax2.plot(np.linspace(0,r[0].stats.npts/r[0].stats.sampling_rate,r[0].stats.npts),r[0].data);
    ax2.set_title("Radial Component Waveform for Station {}".format(station));
    ax3.plot(np.linspace(0,t[0].stats.npts/t[0].stats.sampling_rate,t[0].stats.npts),t[0].data);
    ax3.set_title("Transverse Component Waveform for Station {}".format(station));
    
    
    return fig,(ax1,ax2,ax3)