# PaSTA
This package was built to help with the analysis of arrival time of P and S waves from earthquakes. It is my the final project for my Scientific Computing in Astrophysics class in Fall 2021.

This package was specifically built to perform the type of analysis carried out in [Menke et al's paper](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016GL070918), specifically, where they plot the S and P wave time residuals against each other and measure the slope of the best-fit line. This package also has a couple extra features (TX curve and waveform plotting) that should be useful to body wave seismologists.

The acronym stands for P And S wave Time Analysis (PASTA!).

This package would not have been possible without the help of Professor Maureen Long and her student Yantao Luo of Yale University. Thank you to both of them!


# Installation
To install, download the zip file from the main github page, unzip the file, cd into the directory, and do 'pip install .'

# Data Acquisition
The acquire_data() method communicates with IRIS's (Incorporated Research Institutions for Seismology) public data stores and acquires data from stations (seismometers) that fall within the specifications of the user (such as location) corresponding to earthquake events that also fall within the user's specifications (depth range, magnitude, etc.) It uses a known 1-D Earth velocity model specified by the user (such as ak-135) to predict when the initial P wave arrival will be for each station, and extracts the seismometer's data starting 30 seconds prior to that. The user then specifies at what time after that P wave arrival estimate the data range should be cutoff. While this time range could in theory include the surface wave arrivals, this package was made for body wave seismologists (sorry) and therefore is only built to deal with the P and S phases. For any non seismologists, the gifs below demosntrate the way that P and S waves are propagate, respectively.

[](https://web.ics.purdue.edu/~braile/edumod/waves/Pwave_files/image001.gif)
[](https://web.ics.purdue.edu/~braile/edumod/waves/Swave_files/image001.gif)

# Data Analysis

# Plotting

# Example Data
Since the data acquisition script takes a long time to run (owing to its need to communicate with the IRIS server), I have provided a folder with example data of five events measured by stations in the central Appalachian region between 2017 and 2019.
