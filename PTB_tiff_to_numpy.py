import os
import numpy as np
import pandas as pd 
from skimage import io 

def write_channel(df, r_path, o_path):
    """
    Writes images for each region of interest within a slide

    Parameters:
    df (DataFrame): DataFrame containing information about the ROIs
    r_path (str): path to the .tiff file of the slide
    o_path (str): path where the .npy files should be saved
    """
    
    xw = 512
    mat = io.imread(r_path)

    img = {}
    # For each row in the DataFrame (each ROI)
    # Extract the part of the slide around the ROI
    for spot in df.iterrows():
        y, x = spot[1][["ROICoordinateX", "ROICoordinateY"]].astype(np.int64)
        roi = spot[1][["ROILabel"]].astype(int).values[0]
        
        img[roi] = mat[(x - xw):(x + xw), (y - xw):(y + xw)]

    for x in img:
        np.save(o_path + str(x) + ".npy", img[x])

    del img
    del mat

# Define base path
path = "/mnt/c/Users/Jonathan/Documents/UCLA/Pelligrini/projects/"

# List of slides and their corresponding colors
slides_and_colors = [
    ("1 PTB-22.1", ["Blue-006", "Green-002", "Red-005", "Yellow-007"], "ptb22.1"),
    ("2 PTB-22.3", ["Blue-004", "Green", "Red-003", "Yellow"], "ptb22.3"),
    ("PTB22.2", ["Blue", "Green", "Red", "Yellow"], "ptb22.2"),
    ("B172914-2", ["Blue", "Green", "Red", "Yellow"], "ptb21.1"),
]

# Process each slide
for slide_name, colors, out_dir in slides_and_colors:
    # Load DataFrame
    if slide_name == "B172914-2":
        df = pd.read_excel(path + "geomx/current_data/All Data WTA 2 batches.xlsx", sheet_name="SegmentProperties")
    elif slide_name == "PTB22.2":
        df = pd.read_excel(path + "geomx/current_data/Initial Dataset.xlsx", sheet_name="SegmentProperties")
    else:
        df = pd.read_excel(path+"geomx/current_data/Initial Dataset S-23-0131.xlsx",sheet_name="SegmentProperties")
    
    # Filter DataFrame for this slide
    df = df.loc[df["SlideName"] == slide_name]
    df = df[["ROILabel", "ROICoordinateX", "ROICoordinateY", "AOISurfaceArea"]]

    # Process each color
    for color in colors:
        write_channel(df = df, 
                      r_path = os.path.join(path, "geomx/geomx_image_anal/img/tiff/" + slide_name + "_" + color + ".tiff"),
                      o_path = os.path.join(path, "geomx/geomx_image_anal/crop/" + out_dir + "/" + color.split("-")[0].upper() + "_"))