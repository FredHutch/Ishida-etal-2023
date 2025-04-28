import os
import glob
import numpy as np
from skimage import io, measure
import cv2
from scipy.optimize import linear_sum_assignment

# Global state variables used across frames
_next_new_label = [2]
first_detected_blue_intensity = None
current_frame = 0

def detect_cells(image_path, blue_thresh=220, adaptive_thresh=180, regions=None, height_limit=1200, max_blue_area=50):
    """
    Detect cells and cells divisions. Ensure they are valid cells.

    Args:
      - image_path (str): Path to the image file.
      - blue_thresh (int): Initial blue threshold for first cell detection.
      - adaptive_thresh (int): Local blue threshold used in radius around already detected cells.
      - regions (list of tuples): Optional region centroids for adaptive thresholding.
      - height_limit (int): Limit image processing to top N pixels (to be used to crop the bottom of the image if necessary).
      - max_blue_area (int): Maximum blue size a cell can have.

    Returns a list of dictionaries with keys:
      - centroid: (x, y) - center of a potential cell
      - area: area (in pixels)
      - blue_intensity: cumulative blue channel sum in the cell region
      - red_intensity: cumulative red channel sum in the cell region
    """
    
    print(f"\nProcessing {image_path}...")
    global first_detected_blue_intensity, current_frame
    img = io.imread(image_path)

    if img.ndim < 3:
        raise ValueError("Expected a multi-channel (RGB) image.")

    blue_channel = img[:1100, :, 2]
    red_channel = img[:1100, :, 0]
    frame = 1

    # Search for cells in the area where a cell was previously detected (with a radius of 35)
    if regions:
        frame = 1
        blue_mask = np.zeros_like(blue_channel, dtype=bool)
        for cx, cy in regions:
            x_min = max(0, int(cx) - 35) # 35 is the radious around the existing cell that is searched for a new cell with a lower thresold
            x_max = min(blue_channel.shape[1], int(cx) + 35)
            y_min = max(0, int(cy) - 35)
            y_max = min(blue_channel.shape[0], int(cy) + 35)
            region = blue_channel[y_min:y_max, x_min:x_max]
            blue_mask[y_min:y_max, x_min:x_max] |= (region >= adaptive_thresh)

    # If no cell has been detected on the previous frame, use the "initial blue thresold" on the whole image
    else:
        blue_mask = blue_channel >= blue_thresh

    labeled = measure.label(blue_mask)
    props = measure.regionprops(labeled)

    detections = []

    # Identify the valid cells in the detected centroids list
    for prop in props:
        # Check size of the centroid
        if prop.area > max_blue_area:
            print(f"Rejected cell at {prop.centroid} - Blue area too large: {prop.area} pixels.")
            continue

        # If the cell is elongated, consider it as two cells overlapping
        if prop.minor_axis_length > 0:
            aspect_ratio = prop.major_axis_length / prop.minor_axis_length
        else:
            aspect_ratio = 0
        if aspect_ratio > 2:
            print(f"Elongated cell detected at {prop.centroid} with aspect ratio {aspect_ratio:.2f}. Splitting into two cells.")
            mid_x = (prop.bbox[1] + prop.bbox[3]) / 2
            mid_y = (prop.bbox[0] + prop.bbox[2]) / 2

            offset = prop.major_axis_length / 4
            angle = prop.orientation
            dx = offset * np.cos(angle)
            dy = offset * np.sin(angle)

            centroid1 = (prop.centroid[1] - dx, prop.centroid[0] - dy)
            centroid2 = (prop.centroid[1] + dx, prop.centroid[0] + dy)

            for centroid in [centroid1, centroid2]:
                x, y = int(centroid[0]), int(centroid[1])
                if 0 <= x < blue_channel.shape[1] and 0 <= y < blue_channel.shape[0]:
                    blue_intensity = blue_channel[y, x]
                    red_intensity = red_channel[y, x]
                    detections.append({
                        'centroid': centroid,
                        'area': prop.area / 2,  # Approximate area split
                        'blue_intensity': blue_intensity,
                        'red_intensity': red_intensity
                    })
        
        else:
            cy, cx = prop.centroid
            coords = prop.coords
            print(f"Cell at {prop.centroid}: Area={prop.area}, Pixels Considered={coords.shape[0]}")

            blue_intensity = np.sum(blue_channel[coords[:, 0], coords[:, 1]])
            red_intensity = np.sum(red_channel[coords[:, 0], coords[:, 1]])

            # Save the cummulative blue pixel intensity of the first cell detected of a well (on frame 1)
            if current_frame == 1 and first_detected_blue_intensity is None:
                first_detected_blue_intensity = blue_intensity

            # Ignore if the newly detected centroid has a blue intensity too low compared to its parent (less than 20%)
            if current_frame>= 2 and blue_intensity < 0.2 * first_detected_blue_intensity:
               print(f"Ignoring new detection at ({cx:.2f}, {cy:.2f}) - Blue intensity too low: {blue_intensity} (Threshold: {0.2 * first_detected_blue_intensity})")
               continue 

            detections.append({
                'centroid': (cx, cy),
                'area': prop.area,
                'blue_intensity': blue_intensity,
                'red_intensity': red_intensity
            })      

    # As there is only one cell on the first frame, keep the centroid with the strongest cummulative blue intensity as the cell
    if current_frame == 0 and len(detections) > 1:
        detections.sort(key=lambda x: x['blue_intensity'], reverse=True)
        strongest_cell = detections[0]
        detections = [strongest_cell] 

    # handle teh case when there are multiple potential cells on frame 1
    if current_frame == 1 and len(detections) > 1:
        detections.sort(key=lambda x: x['blue_intensity'], reverse=True)
        strongest_cell = detections[0]
        detections = [strongest_cell]
        first_detected_blue_intensity = strongest_cell['blue_intensity']
        print(f"First detected blue intensity set to: {first_detected_blue_intensity}")
        print(f"Multiple cells detected in Frame 1. Keeping only the strongest: "
              f"Centroid={strongest_cell['centroid']}, Blue Intensity={strongest_cell['blue_intensity']}")

    if detections:
        print(f"Detected {len(detections)} valid cell(s):")
        for idx, cell in enumerate(detections):
            print(f"Cell {idx+1}: Centroid=({cell['centroid'][0]:.2f}, {cell['centroid'][1]:.2f}), "
                  f"Area={cell['area']}, Blue_Intensity={cell['blue_intensity']}, Red_Intensity={cell['red_intensity']}")
    else:
        print("No valid cells detected.")
    
    current_frame += 1
    return detections


def update_labels(old_cells, new_dets, tracking_thresh=40, frame=None):
    """
    Update cell labels across frames based on proximity to model mothers/ daughters relationships.

    This function matches new detections to previously tracked cells and handles:
    - Continuing a cell's label if matched (1-to-1 case)
    - Splitting labels into -1 and -2 suffixes when a cell divides

    Args:
      - old_cells (list of dict): Cells tracked on the previous frame, each with 'centroid' and 'label'.
      - new_dets (list of dict): Newly detected cells to be labeled.
      - tracking_thresh (float): Maximum distance to consider a valid match.
      - frame (int or None): Current frame index, used to record division timing.

    Returns:
      - list of dict: Updated `new_dets`, each with assigned 'label' and optionally a 'division_time'.
    """
    n_old = len(old_cells)
    n_new = len(new_dets)
    
    # Case when there is no division (Hungarian matching)
    if n_new == n_old:
        cost_matrix = np.zeros((n_old, n_new))
        for i, oc in enumerate(old_cells):
            for j, nd in enumerate(new_dets):
                cost_matrix[i, j] = euclidean_distance(oc['centroid'], nd['centroid'])
        row_ind, col_ind = linear_sum_assignment(cost_matrix)
        
        for i, j in zip(row_ind, col_ind):
            new_dets[j]['label'] = old_cells[i]['label']
            new_dets[j]['division_time'] = None
            print(f"Assigned label {new_dets[j]['label']} to cell at {new_dets[j]['centroid']}")
        return new_dets
    
    # Case when there is at least one division
    assignments = {}
    assignment_count = {oc['label']: 0 for oc in old_cells}
    
    # Build a list of potential matching cells for each newly deteted cells
    candidate_lists = {}
    for i, nd in enumerate(new_dets):
        candidates = []
        for oc in old_cells:
            d = euclidean_distance(oc['centroid'], nd['centroid'])
            candidates.append((oc['label'], d))
        candidates.sort(key=lambda x: x[1]) # Sort by shortest distance
        candidate_lists[i] = candidates
    
    # Assign new detections to closest available candidate
    new_indices = list(range(n_new))
    new_indices.sort(key=lambda i: candidate_lists[i][0][1])

    for i in new_indices:
        assigned = False
        for cand_label, d in candidate_lists[i]:
            if assignment_count[cand_label] < 2:
                assignments[i] = cand_label
                assignment_count[cand_label] += 1
                assigned = True
                break
        # If no mother or previous cell can be found
        if not assigned:
            assignments[i] = str(_next_new_label[0])
            _next_new_label[0] += 1

    
    groups = {}
    for i, label in assignments.items():
        groups.setdefault(label, []).append(i)
    
    new_labels = [None] * n_new
    for label, indices in groups.items():
        if label in assignment_count:
            # 1 to 1 direct match
            if len(indices) == 1:
                new_labels[indices[0]] = label
                new_dets[indices[0]]['division_time'] = None
            # Division case
            elif len(indices) == 2:
                new_labels[indices[0]] = label + "-1"
                new_dets[indices[0]]['division_time'] = frame
                new_labels[indices[1]] = label + "-2"
                new_dets[indices[1]]['division_time'] = frame
            # More than two daughters possible - assign two and reassign the others
            else:
                new_labels[indices[0]] = label + "-1"
                new_dets[indices[0]]['division_time'] = frame
                new_labels[indices[1]] = label + "-2"
                new_dets[indices[1]]['division_time'] = frame
                for j in indices[2:]:
                    reassign_label = None
                    for alt_label, d in candidate_lists[j]:
                        if alt_label != label and assignment_count.get(alt_label, 0) < 2:
                            reassign_label = alt_label
                            assignment_count[alt_label] += 1
                            break
                    if reassign_label is None:
                        reassign_label = str(_next_new_label[0])
                        _next_new_label[0] += 1
                    new_labels[j] = reassign_label
                    new_dets[j]['division_time'] = None
        else:
            for i in indices:
                new_labels[i] = assignments[i]
                new_dets[i]['division_time'] = None
                
    for i, nd in enumerate(new_dets):
        nd['label'] = new_labels[i]
        print(f"Assigned label {nd['label']} to cell at {nd['centroid']}")
        
    return new_dets


def process_images(image_folder, output_folder, blue_thresh=220, adaptive_thresh=150,
                   height_limit=1200, max_blue_area=5, tracking_thresh=40):
    """
    Process a sequence of microscopy images to detect and track cells over time.

    For each frame, this function:
      - Detects cells based on blue channel intensity
      - Assigns or updates labels based on proximity to prior cells
      - Logs cell position, intensity, and division events

    Saves outputs to the output folder:
      1. Marked images with drawn cells
      2. Tracking log CSV: frame-by-frame lineage info
      3. Distance log CSV: pairwise distances per frame
      4. Intensity log CSV: area, blue, red intensities per cell

    Args:
      - image_folder (str): Folder containing `.tif` input images for one well. The images need to be sorted in chronological order in the folder. 
      - output_folder (str): Folder to save logs and annotated images.
      - blue_thresh (int): Initial blue threshold for first cell detection.
      - adaptive_thresh (int): Local blue threshold used in radius around already detected cells.
      - height_limit (int): Limit image processing to top N pixels (to be used to crop the bottom of the image if necessary).
      - max_blue_area (int): Maximum blue size a cell can have.
      - tracking_thresh (float): Maximum distance for matching cells between frames.

    """
    
    os.makedirs(output_folder, exist_ok=True)
    image_paths = sorted(glob.glob(os.path.join(image_folder, "*.tif")))
    folder_name = os.path.basename(image_folder)
    print(f"\nProcessing {len(image_paths)} images...")
    tracked_cells = []
    
    tracking_log = []
    distance_log = []
    intensity_log = []
    
    # Process each image
    for t, image_path in enumerate(image_paths):
        print(f"\n--- Frame {t} ---")
        regions = [cell['centroid'] for cell in tracked_cells] if tracked_cells else None
        detections = detect_cells(
            image_path,
            blue_thresh=blue_thresh,
            adaptive_thresh=adaptive_thresh,
            regions=regions,
            height_limit=height_limit,
            max_blue_area=max_blue_area
        )
        # Assign the label of the first cell, on the first frame
        if t == 0:
            if len(detections) == 1:
                detections[0]['label'] = "1"
                detections[0]['division_time'] = None
                print(f"Assigned label {detections[0]['label']} to cell at {detections[0]['centroid']}")
            
            # This case should not happen as there is only one possible cell on the first frame
            else:
                for i, det in enumerate(detections):
                    det['label'] = f"1-{i+1}"
                    det['division_time'] = None
                    print(f"Assigned label {det['label']} to cell at {det['centroid']}")
            tracked_cells = detections
        
        else:
            tracked_cells = update_labels(tracked_cells, detections, tracking_thresh=tracking_thresh, frame=t)
        
        # Collect data for the csv files
        for cell in tracked_cells:
            tracking_log.append([
                t,
                cell['label'],
                cell['centroid'][0],
                cell['centroid'][1],
                cell['blue_intensity'],
                cell['red_intensity'],
                cell.get('division_time')
            ])
            intensity_log.append([
                t,
                cell['label'],
                cell['area'],
                cell['blue_intensity'],
                cell['red_intensity']
            ])
        
        num_cells = len(tracked_cells)
        for i in range(num_cells):
            for j in range(i+1, num_cells):
                cell1 = tracked_cells[i]
                cell2 = tracked_cells[j]
                dist = euclidean_distance(cell1['centroid'], cell2['centroid'])
                distance_log.append([
                    t,
                    cell1['label'],
                    cell1['centroid'][0],
                    cell1['centroid'][1],
                    cell2['label'],
                    cell2['centroid'][0],
                    cell2['centroid'][1],
                    dist
                ])
        
        if tracked_cells:
            marked_path = os.path.join(output_folder, f"marked_time_{t}.png")
            draw_detected_cells(image_path, tracked_cells, marked_path)
        else:
            print(f"No valid cells detected in {os.path.basename(image_path)}. Skipping output image.")

    # Save the csv files
    tracking_csv_path = os.path.join(output_folder, f"{folder_name}_cell_tracking_log.csv")
    with open(tracking_csv_path, "w") as f:
        f.write("Frame,Label,Centroid_X,Centroid_Y,Blue_Intensity,Red_Intensity,Division_Time\n")
        for entry in tracking_log:
            f.write(",".join(map(str, entry)) + "\n")
    print(f"\nTracking log saved as: {tracking_csv_path}")

    distance_csv_path = os.path.join(output_folder, f"{folder_name}_cell_distance_log.csv")
    with open(distance_csv_path, "w") as f:
        f.write("Frame,Cell1_Label,Cell1_Centroid_X,Cell1_Centroid_Y,Cell2_Label,Cell2_Centroid_X,Cell2_Centroid_Y,Distance\n")
        for entry in distance_log:
            f.write(",".join(map(str, entry)) + "\n")
    print(f"Distance log saved as: {distance_csv_path}")

    intensity_csv_path = os.path.join(output_folder, f"{folder_name}_cell_intensity_log.csv")
    with open(intensity_csv_path, "w") as f:
        f.write("Frame,Label,Area,Blue_Intensity,Red_Intensity\n")
        for entry in intensity_log:
            f.write(",".join(map(str, entry)) + "\n")
    print(f"Intensity log saved as: {intensity_csv_path}")


def euclidean_distance(p1, p2):
    """
    Compute the Euclidean distance between two (x,y) points.

    Args: 
      - p1: first point coordinates
      - p2: second point coordinates
    
    Return:
      - Ditstance between point 1 and point 2
    
    """
    return np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)


def draw_detected_cells(image_path, tracked_cells, output_path, circle_radius=10, circle_color=(0,255,0), thickness=2):
    """
    This is purely a visual helper function. It draws circles and labels on the image each tracked cell.
    Args:
      - image_path (str): Path to the image to draw on.
      - tracked_cells (list): List of cell dictionaries with centroid and label.
      - output_path (str): Path to save the output image.
      - circle_radius (int): Radius of the circle to draw.
      - circle_color (tuple): Color of the circle in BGR.
      - thickness (int): Thickness of the circle and label text.
    """
    img = io.imread(image_path)
    img_bgr = cv2.cvtColor(img, cv2.COLOR_RGB2BGR)
    
    # For each detected and valid cell, circle it on an output image.
    for cell in tracked_cells:
        cx, cy = cell['centroid']
        cv2.circle(img_bgr, (int(cx), int(cy)), circle_radius, circle_color, thickness)
        label_text = cell.get('label', '')
        div_time = cell.get('division_time', None)
        if div_time is not None:
            label_text += f" (Div@{div_time})"
        
        cv2.putText(img_bgr, label_text, (int(cx)+circle_radius, int(cy)),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.4, circle_color, 1, cv2.LINE_AA)
    cv2.imwrite(output_path, img_bgr)
    print(f"Marked image saved as: {output_path}")




if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Cell tracking with hierarchical labeling, computed intensities, and division timing.")
    parser.add_argument("--image_folder", required=True, help="Folder containing input images (.tif).")
    parser.add_argument("--output_folder", required=True, help="Folder to save marked images and CSV logs.")
    parser.add_argument("--blue_thresh", type=int, default=180, help="Global blue threshold for initial detection.") #250 for Col3
    parser.add_argument("--adaptive_thresh", type=int, default=40, help="Adaptive threshold for subsequent frames.")
    parser.add_argument("--max_blue_area", type=int, default=50, help="Maximum allowed blue area in pixels.")
    parser.add_argument("--height_limit", type=int, default=1200, help="Only process the top height_limit rows.")
    parser.add_argument("--tracking_thresh", type=int, default=40, help="Distance threshold for matching cells.")
    args = parser.parse_args()

    process_images(
        image_folder=args.image_folder,
        output_folder=args.output_folder,
        blue_thresh=args.blue_thresh,
        adaptive_thresh=args.adaptive_thresh,
        height_limit=args.height_limit,
        max_blue_area=args.max_blue_area,
        tracking_thresh=args.tracking_thresh
    )
