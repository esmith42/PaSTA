# PaSTA
This package was built to help with the analysis of arrival time of P and S waves from earthquakes. It is my the final project for my Scientific Computing in Astrophysics class in Fall 2021.

This package was specifically built to perform the type of analysis carried out in [Menke et al's paper](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016GL070918), specifically, where they plot the S and P wave time residuals against each other and measure the slope of the best-fit line. This package also has a couple extra features (TX curve and waveform plotting) that should be useful to body wave seismologists.

The acronym stands for P And S wave Time Analysis (PASTA!).

This package would not have been possible without the help of Professor Maureen Long and her student Yantao Luo of Yale University. Thank you to both of them!


# Installation
To install, download the zip file from the main github page, unzip the file, cd into the directory, and do 'pip install .'

# Data Acquisition
The acquire_data() function communicates with IRIS's (Incorporated Research Institutions for Seismology) public data stores and acquires data from stations (seismometers) that fall within the specifications of the user (such as location) corresponding to earthquake events that also fall within the user's specifications (depth range, magnitude, etc.) It uses a known 1-D Earth velocity model specified by the user (such as ak-135) to predict when the initial P wave arrival will be for each station, and extracts the seismogram's data starting 30 seconds prior to that. The user then specifies at what time after that P wave arrival estimate the data range should be cutoff. While this time range could in theory include the surface wave arrivals, this package was made for body wave seismologists (sorry) and therefore is only built to deal with the P and S phases. For any non seismologists, the gifs below demosntrate the way that P and S waves are propagate, respectively.

![Alt Text](https://i.stack.imgur.com/wji2Z.gif)

Besides extracting the data, the function also takes the data and stores it in the user's computer in the specified directory. Each event is assigned a code with its corresponding subdirectory having that name. All these subdirectories go into the specified save directory. Each event subdirectory has a csv file called 'evinfo' which contains a table where each row corresponds to a station and the columns are station name, event year, month, day, hour, minute, second, event latitude, longitude, depth (in kilometers), magnitude, great-circle distance from event to station in degrees, event back-azimuth (the direction of the event with respect to the station). Also in these subdirectories are  additional subdirectories corresponding to each station that had data for that event. These contain the SAC files with the data for that station for that event (named STACK_R, STACK_T, STACK_Z for the radial, transverse, and vertical components, respectively). For non-seismologists, SAC files are what seismologists use to store seismogram data (they are to seismologists what FITS files are to astronomers). To help any non-seismologists understand what it means to have different components of a seismometer, see the gifs below for a helpful visual aid.

![Alt Text](https://opengeology.org/textbook/wp-content/uploads/2018/07/09.5-seismograph_vert.gif)


![Alt Text](https://i.makeagif.com/media/12-16-2015/e81Yff.gif)


Here is an example of code one write run to get data from the central Appalachian region (see the documentation for it for the meaning of each of the parameters):

```
acquire_data("2015-08-01","2019-08-01", 5.5, 30, 90, 750, 'ak135', 30, 800, 41.68, 42.18, -73.6, -71.4,'./',"IRIS","IRIS","XP","*","*","BH?")
```

**Note:** The bulk of the code for this part of the package was written by Yantao Luo of Yale University. It was restructured and generalized to be made compatible with this package.
**Warning:** This function can take a loooooonnnngg time to run (like 1-2 hours!)

# Data Analysis


# Plotting

# Example Data
Since the data acquisition script takes a long time to run (owing to its need to communicate with the IRIS server), I have provided a folder with example data of five events measured by stations in the central Appalachian region between 2017 and 2019.
