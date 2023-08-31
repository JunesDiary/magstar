import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from mendeleev import element
import re
import pandas as pd
import os
from IPython.display import clear_output
import shutil
import time
import subprocess
from subprocess import Popen, PIPE
from PIL import Image, ImageTk
import sys
import astropy.units as u
from os import listdir
from os.path import isfile, join


#FUNCTIONS DEFINITIONS

#FLAG COUNTER
#not dependant on directory
def elem():
    
    try:
        specname = element(int(str(word[1])[0:2])).symbol
        specstate = int(str(word[1])[-1]) + 1
        return str(specname) + ' ' + str(specstate)
    except:  
        specname = element(int(str(word[2])[0:1])).symbol
        specstate = int(str(word[2])[-1]) + 1
        return str(specname) + ' ' + str(specstate)
    

def main_flag_counter(valdfile): 
    i = 0
    global word
    flagcount = 0
    #print('Element Wavelength Lande Excit')
    #print('======= ========== ===== =====')
    data = open(valdfile, 'r')
    for line in data:
        i = i + 1
        if i == 1:
            continue
        else:
            word = [x.strip() for x in line.split('  ')]  
            if (int(word[-1])):
                flagcount = flagcount + 1
                #try:
                #    print(elem(), '   ',(word[0]),'',   (word[-2]), (word[-3]))
                #except:
                #    print(elem(), '   ',(word[0]),'',   (word[-2]), (word[-3]))
    data.close()
    return 'Lines Flagged: ' + str(flagcount) + ' == Total Lines: ' + str(i) + ' == Lines Used: '+ str(flagcount*100/i)[:5] + '%' + '\n======================================================'

#======================================================================================
#SET ALL FLAGS ZERO
#not dependant on directory
def set_zero_main(maskfile): 
    i = 0
    count = -1
    
    zero = []
    zero.append('0')
    
    with open(maskfile, 'r') as f:
        lines = f.readlines()    
        f.close()
    os.remove(os.getcwd() + "/mask.txt")
    with open('mask.txt', 'w') as f:
        for line in lines:
            count = count + 1
            
            if count != 0:
                word = [x.strip() for x in line.split('  ')]
               
                f.write(word[0] + '  ')
                f.write(word[1] + '  ')
                f.write(word[2] + '  ')
                f.write(word[3] + '  ')
                f.write(word[4] + '  ')
                f.write(word[5] + '  ')
                f.write(word[6] + '  ')
                f.write(word[7] + '  ')
                f.write(word[8] + '  ')
                f.write('0' + '\n')
                           
        f.close()
    return 
#===================================================================================
#VALD TO MASK CONVERTION
#dependant on directory: till LSDpy
#output labelled
def ion_state(name, ion_no):
    at_no = element(name[0:3]).atomic_number
    return "{:.2f}".format(float(at_no + (float(ion_no)-1)/100))

def vald_to_mask_main(valdfile): 
    
    data = open(valdfile, 'r')
    for line in data:
        word = [x.strip() for x in line.split(' ')]    
        break        
    
    global linesflagged
    linesflagged = 0
    
    header = {}
    #print(word)

    header['start_Wavelength'] = word[1][:-1]
    header['end_Wavelength'] = word[2][:-1]
    header['lines_selected'] = word[3][:-1]
    header['lines_processed'] = word[4][:-1]
    header['vmicro'] = word[5][:-1]

    #print(header)

    data = open(valdfile, 'r')
    i = -1

    for line in data:
        i = i + 1
        if i == 0 or i == 1:
            continue
        elif i == 2:
            word = [x.strip() for x in line.split(' ' or '  ')]
            col_head = [word[0]+word[1], word[8], word[9], word[10], word[11]+word[12][:-1], word[14], word[17], word[19], 'lande'+word[21], 'central'+ word[23]]
            dataout = pd.DataFrame(columns=col_head)
        else:
            word = [x.strip() for x in line.split(',')][0:10]   
            try:
                dataout.loc[i] = word
            except:
                #print(i+1, word) #lines left out
                None
    try:
        os.remove(os.getcwd() + "/mask.txt") 
    except:
        None
                
    mask = open(os.getcwd() + "/mask.txt","a")
    mask.write(header['lines_selected'] + '\n')

    
    for i in range(3,len(dataout)):
        #if i == 10:
         #   break

        #clear_output(wait=True)
        #print(str("{:.2f}".format((i+1)*100/int(header['lines_selected']))) + ' % completed')

        try:

            mask.write(' ' + str("{:.4f}".format((float(dataout.loc[i][1])/10))) + '  ') #wavelength
            try: #species code
                mask.write(str(ion_state(dataout.loc[i][0][1:-1].split(' ')[0], dataout.loc[i][0][1:-1].split(' ')[1])) + '  ')
            except:
                print('Error: def ion_state has problem!')
            mask.write(str(dataout.loc[i][9]) + '  ')    #central depth
            mask.write(str(dataout.loc[i][2]) + '  ')    #excitation level
            mask.write(str((dataout.loc[i][8])) + '  ' )  #lande factor
            mask.write('0 \n') #flag #str(isLineStrong(dataout.loc[i][1], lines))
        except:
            None       

    mask.close()
#================================================================================
#COMPUTE LSD PROFILES
#dependant on directory: till LSDpy
def inlsdedit(text):
    text = (text.split('/'))[-1]
    count = 0
    with open(os.getcwd() + "/inlsd.dat", 'r') as f:
        lines = f.readlines()    
        f.close()
    try:
        os.remove(os.getcwd() + '/inlsd.dat')
    except:
        print('Old inlsd file not found')
        exit()
    with open('inlsd.dat', 'w') as f:
        for line in lines:
            count = count + 1
            if count == 3:
                f.write(text)
            elif count == 5:
                f.write(mask_path.split('/')[-1] + '\n')                
            else:
                f.write(line)

def lsdcomputer():                

    listdir = spectra_file_paths.split('\n')
    
    filelist = []
    status = []
    progress = []
    if os.path.exists(os.getcwd() + "/dataout"):
        shutil.rmtree(os.getcwd() + "/dataout")
    os.makedirs(os.getcwd() + "/dataout")
    
    #create lsd profile image storage folder and remove old one
    if os.path.exists(os.getcwd() + "/lsd_images"):
        shutil.rmtree(os.getcwd() + "/lsd_images")
    os.makedirs(os.getcwd() + "/lsd_images")
    
    for file in listdir:
        #edit the inlsd file and put the spectra name which you will be creating lsd profile of
        inlsdedit(file.split('/')[-1] + '\n')
        
        #remove the mask file and copy the one given by user to lsdpy folder
        try:
            os.remove(os.getcwd() + '/mask.txt')
        except:
            None
        shutil.copy(mask_path, os.getcwd())
        
        filelist.append(file.split('/')[-1])
        status.append('Processed')
        progress.append('('+ str(listdir.index(file)+1) + '/' + str(len(listdir)) + ')')
        
        #copy the spectra to lsdpy folder for lsd to compute its profile
        shutil.copy(file, os.getcwd())

        #calling terminal to create the terminal
        theproc = subprocess.Popen([sys.executable, "lsdpy.py"])
        theproc.communicate()    
        

        
        #copy lsd profile image to prof folder
        shutil.copy(os.getcwd() + '/lsdproftemp.png', os.getcwd() + '/lsd_images/lsdprof_' + file.split('/')[-1])
        
        #delete lsd profile image from lsdpy folder
        os.remove(os.getcwd() + '/lsdproftemp.png')
        
        #wait time for lsd compution to be over
        time.sleep(3)

        #try to create dataout folder if doesnot exist

        
        #copy the profile file to dataout folder
        shutil.copy(os.getcwd() + "/prof.dat", os.getcwd() + "/dataout/prof_" + file.split('/')[-1])

        #remove the spectra field from lsdpy folder for next spectra to be loaded
        os.remove(os.getcwd() + '/' + file.split('/')[-1])
    
    table2 = pd.DataFrame(columns = ['File', 'Status', 'Progress'])
    table2['File'] = filelist
    table2['Status'] = status
    table2['Progress'] = progress
    result_string = "Computed LSD Profiles:\n\n" + table2.to_string(index=False)
    return result_string + '\n======================================================'

# Function to open a new window and display an image
def open_images_window():
    image_paths = os.listdir(os.getcwd() + '/lsd_images')
    images = []

    for image_path in image_paths:
        img = Image.open(os.getcwd() + '/lsd_images/' + image_path)
        img = img.resize((600, 450), Image.ANTIALIAS)  # Resize the images

        img_tk = ImageTk.PhotoImage(img)
        images.append(img_tk)

    images_window = tk.Toplevel(root)
    images_canvas = tk.Canvas(images_window)
    images_scrollbar = tk.Scrollbar(images_window, orient=tk.VERTICAL, command=images_canvas.yview)
    images_frame = tk.Frame(images_canvas)

    images_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
    images_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
    images_canvas.config(yscrollcommand=images_scrollbar.set)
    images_canvas.create_window((0, 0), window=images_frame, anchor=tk.NW)

    images_frame.bind("<Configure>", lambda event: images_canvas.configure(scrollregion=images_canvas.bbox("all")))

    for img_tk in images:
        image_label = tk.Label(images_frame, image=img_tk)
        image_label.image = img_tk  # Keep a reference to prevent the image from being garbage collected
        image_label.pack(side=tk.TOP, padx=10, pady=5)
        
#================================================================================      
#Bz COMPUTER
#depends on specpolflow location only
def specpolflow_main(direc):
    sys.path.append(specpolflowfolderpath)
    import specpolFlow as pol
    mag = []
    error = []
    spectra = []
    
    direc = str(direc)
    starname = []
    files = [f for f in listdir(direc) if isfile(join(direc, f))]

    for i in range(0,len(files)):
        header = files[i]
        
        lsd = pol.iolsd.read_lsd(direc + '/' + files[i]) #dataout/2612135pu.con_nm.dat
        fig, ax = lsd.plot()
        
        Bz, fig = pol.bz.calcBz(lsd, norm='auto', cog='I', 
                       velrange=[-100,100],bzwidth=30, 
                       geff=1.2, lambda0=5000*u.AA,
                       plot=True, )    
        Bz = pd.DataFrame(data=[Bz])

        # simple display of the pandas dataframe
        Bz.style
        # To extract the entry for the Bz column:
        print('Bz=',Bz.at[0,"V bz (G)"],'+-',Bz.at[0,"V bz sig (G)"],'Gauss')

        # The 0 is the index of the row (there is only one row)
        
        spectra.append(files[i][:-4])
        mag.append(Bz.at[0,"V bz (G)"])
        error.append('+-' + str(Bz.at[0,"V bz sig (G)"]))
    
    table = pd.DataFrame(columns = ['Spectra', 'V bz (G)', 'V bz sig (G)'])
    table['Spectra'] = spectra
    table['V bz (G)'] = mag
    table['V bz sig (G)'] = error
    result_string = table.to_string(index=False)
    return 'Longitudinal Magnetic Field:\n\n' + result_string + '\n======================================================'
#==================================================================================
#==================================================================================


class Tooltip:
    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.tooltip = None
        self.widget.bind("<Enter>", self.show_tooltip)
        self.widget.bind("<Leave>", self.hide_tooltip)
        
    def show_tooltip(self, event):
        x, y, _, _ = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 25

        self.tooltip = tk.Toplevel(self.widget)
        self.tooltip.wm_overrideredirect(True)
        self.tooltip.wm_geometry(f"+{x}+{y}")
        tk.Label(self.tooltip, text=self.text, bg="yellow").pack()

    def hide_tooltip(self, event):
        if self.tooltip:
            self.tooltip.destroy()
            self.tooltip = None

def choose_spectra():
    file_paths = filedialog.askopenfilenames()
    global spectra_file_paths
    if file_paths:
        selected_files.set("\n".join(file_paths))
    # Get the selected file paths
    spectra_file_paths = selected_files.get()   
    if spectra_file_paths:
        return spectra_file_paths + '\n\n' + str(len(file_paths)) + ' Spectra Files Choosen' + '\n======================================================'
def choose_mask_file():
    global mask_path    
    mask_file_path = filedialog.askopenfilename(title='Select MASK File', filetypes=[("Text Files", "*.txt")])
    if mask_file_path:
        mask_file.set(mask_file_path)
    mask_path = mask_file.get()
    return 'Mask File Selected from: ' + mask_path + '\n======================================================' 
    
def choose_vald_file():
    global vald_path
    vald_file_path = filedialog.askopenfilename(title='Select VALD File', filetypes=[("Text Files", "*.txt")])
    if vald_file_path:
        vald_file.set(vald_file_path)
    vald_path = vald_file.get()
        
def install_specpolflow():
    global specpolflowfolderpath
    specpolflowfolderpath = filedialog.askdirectory(title="Select Specpolflow Download Directory")
    if specpolflowfolderpath:
        return specpolflowfolderpath + '\n======================================================'

def install_lsdpypackage():
    # installation logic for LSDpy
    pass

def convert_vald_to_mask():
    try:
        os.remove(os.getcwd() + "/mask.txt")
    except:
        None
    vald_to_mask_main(vald_path)
    # Ask the user where to save the converted MASK file
    save_path = filedialog.asksaveasfilename(
        title="Save Converted MASK File",
        filetypes=[("Text Files", "*.txt")],
        defaultextension=".txt"
    )
    
    if save_path:
        shutil.copy("mask.txt", save_path)
        return "Converted MASK file saved at:\n" + save_path + '\n======================================================'
    else:
        return "Conversion completed, but no file was saved."
    os.remove(os.getcwd() + "/mask.txt")
def count_flagged_lines():
    result = main_flag_counter(mask_path)
    return result


def set_flag_zero_in_mask():
    set_zero_main(mask_path)
    
    # Ask the user where to save the converted MASK file
    save_path = filedialog.asksaveasfilename(
        title="Save MASK File with Flags Set to Zero",
        filetypes=[("Text Files", "*.txt")],
        defaultextension=".txt"
    )
    
    if save_path:
        shutil.copy("mask.txt", save_path)
        return "MASK File with Flags Set to Zero saved at:\n" + save_path + '\n======================================================'
    else:
        return "Conversion completed, but no file was saved."
    os.remove(os.getcwd() + "/mask.txt")
def compute_lsd_profiles():
    return lsdcomputer()

def compute_bz():
    result = specpolflow_main(os.getcwd() + '/dataout')
    return result 
    
def display_bz_variability():
    # logic to display Bz variability
    pass


# Function to display the results 
def display_results(result):
    current_text = result_text.get("1.0", "end-1c")  # Get the current text in the Text widget
    
    if current_text:
        new_text = current_text + "\n\n\n" + result  # Append the new result with a newline character
    else:
        new_text = result
    #print(new_text)    
    result_text.delete("1.0", "end")  # Clear the Text widget
    result_text.insert("1.0", new_text)  # Insert the updated text with the new result
    result_text.see("end")  # Scroll to the end
    
    
# Create the main window
root = tk.Tk()
root.title("MagStar")

# Create a label
label = tk.Label(root, text="Select Spectra Files:")
label.pack(pady=10)

# Create a button to choose spectra files with a tooltip
choose_button = tk.Button(root, text="Choose Spectra", command=lambda: display_results(choose_spectra()))
choose_button.pack()

# Display the selected spectra files
selected_files = tk.StringVar()
#selected_files_label = tk.Label(root, textvariable=selected_files)
#selected_files_label.pack(pady=10)

# Create buttons to choose mask and VALD files
mask_file = tk.StringVar()
vald_file = tk.StringVar()

choose_mask_button = tk.Button(root, text="Choose Mask File", command=lambda: display_results(choose_mask_file()))
choose_mask_button.pack()

choose_vald_button = tk.Button(root, text="Choose VALD File", command=choose_vald_file)
choose_vald_button.pack()

# Create a label
labelmask = tk.Label(root, text="MASK File:")
labelmask.pack(pady=10)

# Display the selected mask files
mask_entry = tk.Entry(root, textvariable=mask_file, bg="white", width=100)
mask_entry.pack(pady=5)

# Create a label
labelvald = tk.Label(root, text="VALD File:")
labelvald.pack(pady=10)

# Display the selected VALD files
vald_entry = tk.Entry(root, textvariable=vald_file, bg="white", width=100)
vald_entry.pack(pady=5)

# Create a label frame to group the additional buttons
button_frame = ttk.LabelFrame(root, text="Functions")
button_frame.pack(pady=20, padx=10)

# Add Specpolflow Installation button
specpolflow_button = tk.Button(button_frame, text="Setup Specpolflow", command=lambda: display_results("Specpolflow set up from: " + install_specpolflow()))
specpolflow_button.pack(pady=5, padx=10, fill=tk.X)

# Add Convert VALD to MASK button
convert_vald_button = tk.Button(button_frame, text="Convert VALD to MASK", 
                                command=lambda: display_results(convert_vald_to_mask()))
convert_vald_button.pack(pady=5, padx=10, fill=tk.X)

# Add Flagged Lines Counter button
flagged_lines_button = tk.Button(button_frame, text="Flagged Lines Counter", command=lambda: display_results(count_flagged_lines()))
flagged_lines_button.pack(pady=5, padx=10, fill=tk.X)

# Add Set Flag Zero in MASK button
set_flag_zero_button = tk.Button(button_frame, text="Set Flag Zero in MASK", command=lambda: display_results(set_flag_zero_in_mask()))
set_flag_zero_button.pack(pady=5, padx=10, fill=tk.X)

# Add Compute LSD Profile(s) button
compute_lsd_button = tk.Button(button_frame, text="Compute LSD Profile(s)", command=lambda: display_results(compute_lsd_profiles()))
compute_lsd_button.pack(pady=5, padx=10, fill=tk.X)

# Add Compute Bz button
compute_bz_button = tk.Button(button_frame, text="Compute Bz using specpolflow", command=lambda: display_results(compute_bz()))
compute_bz_button.pack(pady=5, padx=10, fill=tk.X)

# Add Display Bz Variability button
display_bz_button = tk.Button(button_frame, text="Display Bz Variability", command=display_bz_variability)
display_bz_button.pack(pady=5, padx=10, fill=tk.X)

# Create a Text widget for displaying function results
result_text = tk.Text(root, height=10, width=100, bg="white")
result_text.pack(pady=20)


# Create a button to open the image window
open_button = tk.Button(root, text="Display generated LSD profile plots", command=open_images_window)
open_button.pack(pady=20)

# Create a tooltip for the "Choose Spectra" button
spectra_tooltip = Tooltip(choose_button, "Select Photometric Data from ESPaDOnS")
# Create a tooltip for the "Choose Spectra" button
spectra_tooltip = Tooltip(choose_vald_button, "Select VALD File, as it is downloaded, in short format")
# Create a tooltip for the "Choose Spectra" button
spectra_tooltip = Tooltip(choose_mask_button, "Select MASK File to be given as input to LSD Code")
# Create a tooltip for the "Choose Spectra" button
spectra_tooltip = Tooltip(specpolflow_button, "Install Specpolflow in your working directory")
# Create a tooltip for the "Choose Spectra" button
spectra_tooltip = Tooltip(convert_vald_button, "Convert VALD File format to LSDpy Supported MASK File format, with no lines flagged")
# Create a tooltip for the "Choose Spectra" button
spectra_tooltip = Tooltip(flagged_lines_button, "Count the number of lines flagged in the MASK File")
# Create a tooltip for the "Choose Spectra" button

spectra_tooltip = Tooltip(set_flag_zero_button, "Set all flags in the MASK file as zero")
# Create a tooltip for the "Choose Spectra" button
spectra_tooltip = Tooltip(compute_lsd_button, "Compute LSD Profile(s)")
# Create a tooltip for the "Choose Spectra" button
spectra_tooltip = Tooltip(compute_bz_button, "Compute Longitudinal Magnetic Field from LSD Profiles")
# Create a tooltip for the "Choose Spectra" button
spectra_tooltip = Tooltip(display_bz_button, "Plot Bz versus Phase of Star")


# Start the main event loop
root.mainloop()
