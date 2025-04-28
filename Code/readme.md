# Cell Tracking Pipeline for Time-Lapse Microscopy

This codebase implements a cell tracking pipeline to detect, track, and analyze cells over time using their blue channel intensity (e.g., total protein content) and red channel intensity (e.g., lysosomal activity). The output includes annotated images and detailed logs for downstream biological analysis.


## Overview

This pipeline takes a folder of sequential `.tif` images from a single microscopy well (e.g., from day 0, hour 0 to the last day of the experiment) and outputs:

- Detected cell locations and labels  
- Tracked cell lineages across time  
- Cumulative blue and red intensity measurements  
- Detected division events (based on blue intensity drop)   
- Three CSV logs: lineage, intensity, and pairwise distance  


## Dependencies

This pipeline uses the following Python packages:

| Package         | Purpose                            |
|------------------|------------------------------------|
| `skimage`        | Image reading and region labeling  |
| `cv2` (OpenCV)   | Image annotation and drawing       |
| `numpy`          | Numerical operations               |
| `scipy`          | Hungarian algorithm for matching   |
| `os`, `glob`     | File handling and iteration        |

Install dependencies with:

```bash
pip install numpy scipy scikit-image opencv-python
```


## Code Structure

### 1. Utility Functions & Globals

- `euclidean_distance`: Computes Euclidean distance between two cell centroids  
- `_next_new_label`: Internal counter to assign unique cell labels  
- `current_frame`: Keeps track of the frame index during processing  
- `first_detected_blue_intensity`: Stores the blue intensity of the first detected cell (used as reference for division detection)


### 2. Core Functionalities

#### `detect_cells(image_path, ...)`
Identifies cell-like pixels groups based on blue channel intensity (using a pixel mask) and size thresholds, and extracts properties like centroid, area, and color intensity. It also handles cell division detection by searching for new cells (groups of pixels having at least 20% of the first detected cell’s cumulative blue intensity) in a radius around existing detected cells and shape analysis (e.g. elongated cells which implies two cells are stuck together)

#### `update_labels(old_cells, new_dets, ...)`
Matches newly detected cells with previously tracked ones using the Hungarian algorithm (to assign cells labels based on shortest euclidean distances), updating labels accordingly

#### `draw_detected_cells(image_path, tracked_cells, output_path, ...)`
Draws a circle around detected cells on images and labels them for visual tracking

#### `process_images(image_folder, output_folder, ...)`
Main driver function that iterates through images in chronological order. It calls `detect_cells`, `update_labels`, `draw_detected_cells` and saves 3 CSV logs:
  - `*_cell_tracking_log.csv` — cell label, centroid, intensities, and division time  
  - `*_cell_distance_log.csv` — distances between all cells per frame  
  - `*_cell_intensity_log.csv` — cell area and color intensity values per frame


### 3. Script Entry Point
Uses `argparse` to accept command-line arguments for input directories (a folder with all the images for one well from day 0, hour 0 to the last day of the experiment) and configurable parameters like thresholds and size limits. It then calls `process_images()` accordingly.
The script can be run from the command line using:

```bash
python cell_tracking.py \
  --image_folder path/to/images \
  --output_folder path/to/output \
  --blue_thresh 220 \
  --adaptive_thresh 40 \
  --max_blue_area 50 \
  --height_limit 1200 \
  --tracking_thresh 40
```

**Important Note: the images need to be sorted in chronological order in the input folder (i.e. the name of the first image should be lexicography inferior to the name of the second image ect.)**

#### Optional Parameters Description

| Flag               | Default | Description                                                  |
|--------------------|---------|--------------------------------------------------------------|
| `--blue_thresh`     | 180     | Global threshold for blue channel detection                  |
| `--adaptive_thresh` | 40      | Local threshold for adaptive region checking                 |
| `--max_blue_area`   | 30      | Maximum pixel area for a region to count as a cell           |
| `--height_limit`    | 1200    | Process only the top rows of the image                       |
| `--tracking_thresh` | 40      | Max distance for linking detections across frames            |

---

## Output Files

For each well, the following files will be saved to the `output_folder`:

- `marked_time_<t>.png`: Annotated frame with detected cells and labels  
- `<well_name>_cell_tracking_log.csv`: Main tracking log  
- `<well_name>_cell_distance_log.csv`: Pairwise distances between cells  
- `<well_name>_cell_intensity_log.csv`: Cell area, blue/red intensity per frame  
