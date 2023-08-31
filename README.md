# Magstar

## Description

A GUI based Least Square Deconvolution and Mean Longitudinal Magnetic Field Analysis Tool. The LSDpy computation is done using LSDpy code from folsomcp and magnetic field determination with the help of package specpolflow.
It has the following functionality:

1. Convertion of VALD line file downloaded from VALD3 in short format, to LSDpy supported input Mask format.
2. Computation of LSD profiles and viewing the profiles, for a single file or multiple files.
3. Finding the Mean Longitudinal Magnetic Field (Bz) for the LSD profiles computed.
4. Mask file functionalities like: Calculating number of lines flagged in the Mask file, Setting all flags to zero in the mask file.

Work on the following features are ongoing:

1. Display the Bz variability with phase of star.
2. Automatic flagging of mask file with te help of spectral file and vald line file.

## Requirements

Following packages are required:

1. Tkinter
2. Mendeleev
3. re
4. pandas
5. os
6. shutil
7. time
8. subprocess
9. PIL
10. sys
11. astropy

## Setup

Download the folder LSDpy above. Also download specpolflow from here: [Specpolflow]([url](https://github.com/folsomcp/specpolFlow))
To run the GUI, go to the folder LSDpy you just downloaded and open terminal there and type: 

`python3 magstar.py`


## Operations:

### For finding LSD Profile

1. Select the Spectra files using the 'Choose Spectra' button. You can select many files or a single one too. You will get the message in the result box as: 'n Spectra Files Choosen', where n is the nu,ber of files you chose.
2. If you have downloaded the VALD line file choose the VALD line file by pressing 'Choose VALD File' and then convert it into lsdpy supported mask file by pressing 'Convert VALD to MASK'. Then you need to flag the lines you need in your lsd computation. If you have MASK file with lines flagged already, you can skip this step.
3. Then choose the Mask file you wish to use by clicking button 'Choose Mask File'.
4. Click 'Compute LSD Profile(s)' button. And it will store the LSD Profile text file in working directory in the dataout folder. It will also store 
