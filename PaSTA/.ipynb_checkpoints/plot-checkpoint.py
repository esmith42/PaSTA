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
    

def histogram_residuals(events,phase,**kwargs):
    '''
    This function creates a histogram of the time residuals for the inputted list of event objects for the phase specified (either 
    P or S).
        
    Parameters
    ----------
    events: list of Event objects
        A list of Event objects whose time residuals will be extracted and plotted.
    phase: str
        String specifying the desired phase to be plotted (either 'S' or 'P').
            
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
        
    ax.hist(residuals,**kwargs);
    if phase=='P':
        ax.set_title("P wave travel time residuals");
    else:
        ax.set_title("S wave travel time residuals");
        
    return fig,ax
    
    
    
def plot_residuals(events,line=True,**kwargs):
    '''
    This function creates a plot of the S wave time residuals versus corresponding the P wave time residuals. If the line argument
    is set to True, it will also fit and plot a line to these points with orthogonal distance regression.
        
    Parameters
    ----------
    events: list of Event objects
        A list of Event objects whose time residuals will be extracted and plotted.
    line: boolean
        If this value is set to True, a line will be fitted and plotted over these points with orthogonal distance regression.
            
    Returns
    -------
    fig: Figure
        The figure object corresponding to the plot.
    ax: Axes
        The axes object corresponding to the plot.
    '''
    
    residuals=np.asarray([])
    for event in events:
        p_resid,s_resid = event.get_residuals()
        p_residuals=np.concatenate((residuals,np.squeeze(p_resid)))
        s_residuals=np.concatenate((residuals,np.squeeze(s_resid)))
        
    
    fig,ax = plt.subplots(figsize=(13,13));
    
    plt.rcParams.update({
        "font.size": 25})
       
    ax.tick_params(length=15,width=1,labelsize=25,direction='inout');
    ax.set_xlabel("P wave travel time residuals");
    ax.set_ylabel("S wave travel time residuals");    
    ax.scatter(p_residuals,s_residuals,**kwargs);
    if line:
        beta = fit_line(p_residuals,s_residuals)
        y = p_residuals*beta[0] + beta[1]
        ax.plot(p_residuals,y);

    
    return fig,ax



def fit_residuals(events):
    '''
    This function applies orthogonal least squares regression to the corresponding P and S wave time residuals for the inputted
    list of Event objects. Returns the list of the slope and the intercept of said line.
        
    Parameters
    ----------
    events: list of Event objects
        A list of Event objects whose time residuals will be fitted to a line with orthogonal least squares regression.
            
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
        
    beta = fit_line(p_residuals,s_residuals)
    return beta

def plot_tx(events,phase,time_reduction_fraction=0,**kwargs):
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
    if phase == 'P':
        ax.set_title("P wave time curve");
        for event in events:
            observed_arrivals = event.get_observed_arrival_times_P()                
            ax.scatter(self.distance_to_stations,observed_arrivals - (self.distance_to_stations*time_reduction_fraction),**kwargs);
    elif phase == 'S':
        ax.set_title("S wave time curve");
        for event in events:
            observed_arrivals = event.get_observed_arrival_times_S()                
            ax.scatter(self.distance_to_stations,observed_arrivals - (self.distance_to_stations*time_reduction_fraction),**kwargs);
        
    return fig,ax