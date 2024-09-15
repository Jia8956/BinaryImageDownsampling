import os
import cv2
from PIL import Image
import os
import sys
import torch
import re
import torchvision.transforms as transforms
from gtda.images import Binarizer
import numpy as np
from gtda.homology import CubicalPersistence
from gtda.diagrams import Scaler
import pandas as pd
import random
import math
from collections import defaultdict
from scipy.spatial import cKDTree
from functools import lru_cache
from numba import jit
from typing import Tuple, List, Optional, Union

CONFIG = {
    'THRESHOLD': 26,
    'WHITE_COLOR': [255, 255, 255],
    'BLACK_COLOR': [0, 0, 0],
    'SCALE' : 4,
    'SIZE' : 128,
    'NEIGHBORHOOD': np.array([(dx, dy) for dx in range(-1, 2) for dy in range(-1, 2)])
}

def should_ignore(feature):

    """
    Determine if a topological feature should be ignored.

    Args:
    feature (array-like): Array of [birth, death] values.

    Returns:
    bool: True if the feature should be ignored, False otherwise.

    Ignores features with zero persistence (birth = death = 0).
    """
    birth = feature[0]
    death = feature[1]
    return birth == 0 and death == 0

def cal_betti(image):

    """
    Calculate Betti numbers (β0 and β1) for a given image using cubical persistence.

    Args:
    image (PIL.Image): Input image.

    Returns:
    tuple: (β0, β1) - Betti numbers for 0-dimensional and 1-dimensional homology.

    """
    # Image preprocessing

    to_tensor = transforms.Compose([transforms.ToTensor()])
    binarizer = Binarizer(threshold=0.9)
    cubical_persistence = CubicalPersistence(reduced_homology=False, n_jobs=-1)

    # Convert and invert image
    image_tensor = to_tensor(image).unsqueeze(0)
    inverted_image = 1 - image_tensor  # (0: black, 1: white)

    # Convert to numpy and reshape
    np_image = inverted_image.squeeze().numpy()[None, :, :]
    
    # Binarize and compute persistence
    binary_image = binarizer.fit_transform(np_image)
    persistence_diagram = cubical_persistence.fit_transform(binary_image)

    # Calculate Betti numbers
    betti_0 = np.sum(np.logical_and(persistence_diagram[:, :, 2] == 0, 
                     np.logical_not(np.apply_along_axis(should_ignore, 2, persistence_diagram))))
    betti_1 = np.sum(np.logical_and(persistence_diagram[:, :, 2] == 1, 
                     np.logical_not(np.apply_along_axis(should_ignore, 2, persistence_diagram))))

    return betti_0, betti_1


@jit(nopython=True)
def bilinear_interpolation(image: np.ndarray, x: float, y: float) -> np.ndarray:
    """執行雙線性插值"""
    height, width = image.shape[:2]
    x1, y1 = int(x), int(y)
    x2, y2 = min(x1 + 1, width - 1), min(y1 + 1, height - 1)
    
    dx, dy = x - x1, y - y1
    
    a = image[y1, x1] * (1 - dx) * (1 - dy)
    b = image[y1, x2] * dx * (1 - dy)
    c = image[y2, x1] * (1 - dx) * dy
    d = image[y2, x2] * dx * dy
    
    return a + b + c + d

def preprocess_image(input_image: Union[str, np.ndarray]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Preprocess image: read (if necessary), convert to grayscale, and binarize.

    Args:
        input_image (Union[str, np.ndarray]): Either a file path to the image or a numpy array of the image.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]: Original image, grayscale image, and thresholded image.

    Raises:
        FileNotFoundError: If the input is a file path and the file is not found or cannot be read.
        ValueError: If the input is neither a valid file path nor a numpy array.
    """
    if isinstance(input_image, str):
        image = cv2.imread(input_image)
        if image is None:
            raise FileNotFoundError(f"Image file not found or unable to read: {input_image}")
    elif isinstance(input_image, np.ndarray):
        image = input_image
    else:
        raise ValueError("Input must be either a file path or a numpy array")

    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    _, thresh = cv2.threshold(gray, CONFIG['THRESHOLD'], 255, cv2.THRESH_BINARY)
    
    return image, thresh


def find_contours(thresh: np.ndarray, 
                  retrieval_mode: int = cv2.RETR_EXTERNAL, 
                  sort: bool = True) -> Tuple[List[np.ndarray], Optional[np.ndarray]]:
    """
    Find contours in the thresholded image.

    Args:
        thresh (np.ndarray): Thresholded input image.
        retrieval_mode (int): Contour retrieval mode. Default is cv2.RETR_EXTERNAL.
        sort (bool): Whether to sort contours by area. Default is True.

    Returns:
        Tuple[List[np.ndarray], Optional[np.ndarray]]: 
            - List of contours
            - Hierarchy (if available, else None)
    """
    contours, hierarchy = cv2.findContours(thresh, retrieval_mode, cv2.CHAIN_APPROX_SIMPLE)
    
    if sort:
        contours = sorted(contours, key=cv2.contourArea)
    
    return contours, hierarchy[0] if hierarchy is not None else None

@lru_cache(maxsize=None)
def calculate_distances(contour_tuple, point):
    """Calculate distances from contour points to a given point"""
    contour_array = np.array(contour_tuple)
    return np.linalg.norm(contour_array - point, axis=1)

def find_nearest_pixel(image, center, color):
    """Find the nearest pixel of a specified color"""
    y, x = np.ogrid[:image.shape[0], :image.shape[1]]
    dist = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    mask = np.all(image == color, axis=-1)
    dist[~mask] = np.inf
    nearest_y, nearest_x = np.unravel_index(np.argmin(dist), dist.shape)
    return nearest_x, nearest_y

def process_initial_contours(image, contours, hierarchy):
    """Process contours, extract center points and parent-child relationships"""
    parent_indices = [idx for idx, h in enumerate(hierarchy) if h[3] == -1]
    connected = []
    parent_child_map = defaultdict(list)

    for idx, contour in enumerate(contours):
        has_child = hierarchy[idx][2] != -1
        if not has_child or idx not in parent_indices:
            x, y, w, h = cv2.boundingRect(contour)
            contour_tuple = tuple(map(tuple, contour.squeeze(axis=1)))
            
            distances = calculate_distances(contour_tuple, (x, y))
            min_distance_top_left = contour_tuple[np.argmin(distances)]
            
            distances = calculate_distances(contour_tuple, (x + w, y + h))
            min_distance_bottom_right = contour_tuple[np.argmin(distances)]

            center_x = (min_distance_top_left[0] + min_distance_bottom_right[0]) // 2
            center_y = (min_distance_top_left[1] + min_distance_bottom_right[1]) // 2
            center_point = (center_x, center_y)

            if idx in parent_indices:
                if not np.all(image[center_y, center_x] == CONFIG['WHITE_COLOR']):
                    center_point = find_nearest_pixel(image[y:y+h, x:x+w], (center_x-x, center_y-y), CONFIG['WHITE_COLOR'])
                    center_point = (center_point[0] + x, center_point[1] + y)
                connected.append(center_point)
            
            else:
                parent_index = hierarchy[idx][3]
                if not np.all(image[center_y, center_x] == CONFIG['BLACK_COLOR']):
                    center_point = find_nearest_pixel(image[y:y+h, x:x+w], (center_x-x, center_y-y), CONFIG['BLACK_COLOR'])
                    center_point = (center_point[0] + x, center_point[1] + y)
                parent_child_map[parent_index].append((round(center_point[0] / CONFIG['SCALE']), round(center_point[1] / CONFIG['SCALE'])))

    return connected, parent_child_map

def generate_squares(parent_child_map, connected, size):
    """Generate squares based on parent-child relationships"""
    squares = []
    hole_squares = []
    downsample_img = np.zeros((size, size, 3), dtype=np.uint8)

    if squares:
        tree = cKDTree(np.array(squares))
    else:
        tree = cKDTree(np.array([[0, 0]]))  # 使用一個佔位符點

    for _, value in sorted(parent_child_map.items(), key=lambda item: len(item[1]), reverse=True):
        draw_multiple = len(value) > 1


        if draw_multiple:
            temp = []
            for p in value:
                x, y = p
                hole_bounding_box = np.array([[x - 1, y - 1], [x + 1, y + 1]])
                while tree.query_ball_point(hole_bounding_box.mean(axis=0), r=2):
                    new_x, new_y = CONFIG['NEIGHBORHOOD'][np.random.choice(len(CONFIG['NEIGHBORHOOD']))]
                    x, y = x + new_x, y + new_y
                    hole_bounding_box = np.array([[x - 1, y - 1], [x + 1, y + 1]])
                hole_squares.append(hole_bounding_box)
                temp.append((x, y))
                #tree = cKDTree(np.vstack([tree.data, hole_bounding_box]))

            top_left_x, top_left_y = np.min(temp, axis=0)
            bottom_right_x, bottom_right_y = np.max(temp, axis=0)
            bounding_box = np.array([[top_left_x - 2, top_left_y - 2], [bottom_right_x + 2, bottom_right_y + 2]])
            while tree.query_ball_point(bounding_box.mean(axis=0), r=4):
                new_x, new_y = CONFIG['NEIGHBORHOOD'][np.random.choice(len(CONFIG['NEIGHBORHOOD']))]
                bounding_box += np.array([new_x, new_y])

            squares.append(bounding_box)
            tree = cKDTree(np.vstack([tree.data, bounding_box]))

            cv2.rectangle(downsample_img, tuple(bounding_box[0]), tuple(bounding_box[1]), CONFIG['WHITE_COLOR'], thickness=cv2.FILLED)
            for point in temp:
                cv2.circle(downsample_img, point, 0, CONFIG['BLACK_COLOR'], thickness=-1)
        else:
            x, y = value[0]
            bounding_box = np.array([[x - 2, y - 2], [x + 2, y + 2]])
            while tree.query_ball_point(bounding_box.mean(axis=0), r=4):
                new_x, new_y = CONFIG['NEIGHBORHOOD'][np.random.choice(len(CONFIG['NEIGHBORHOOD']))]
                bounding_box += np.array([new_x, new_y])
            squares.append(bounding_box)
            tree = cKDTree(np.vstack([tree.data, bounding_box]))
            cv2.rectangle(downsample_img, tuple(bounding_box[0] + [1, 1]), tuple(bounding_box[1] - [1, 1]), CONFIG['WHITE_COLOR'], thickness=cv2.FILLED)
            cv2.circle(downsample_img, tuple(bounding_box.mean(axis=0).astype(int)), 0, CONFIG['BLACK_COLOR'], thickness=-1)

    for point in connected:
        x, y = round(point[0] / CONFIG['SCALE']), round(point[1] / CONFIG['SCALE'])
        bounding_box = np.array([[x - 1, y - 1], [x + 1, y + 1]])
        while tree.query_ball_point(bounding_box.mean(axis=0), r=2):
            new_x, new_y = CONFIG['NEIGHBORHOOD'][np.random.choice(len(CONFIG['NEIGHBORHOOD']))]
            bounding_box += np.array([new_x, new_y])
        squares.append(bounding_box)
        tree = cKDTree(np.vstack([tree.data, bounding_box]))
        cv2.circle(downsample_img, tuple(bounding_box.mean(axis=0).astype(int)), 0, CONFIG['WHITE_COLOR'], thickness=-1)

    return downsample_img

def initial(file_path):
    """Main function to process the image and generate the result"""

    image, thresh = preprocess_image(file_path)
    contours, hierarchy = find_contours(thresh, retrieval_mode=cv2.RETR_TREE, sort = False)
    
    if contours is None or hierarchy is None:
        raise ValueError("No contours found in the image")

    
    connected, parent_child_map = process_initial_contours(image, contours, hierarchy)
    downsample_img = generate_squares(parent_child_map, connected, CONFIG['SIZE'])

    return downsample_img

def process_betti0_contours(contours: np.ndarray) -> np.ndarray:
    """Process contours and create bounding boxes"""
    boxes = []
    for contour in contours:
        x, y, w, h = cv2.boundingRect(contour)
        top_x = int((x - CONFIG['SCALE']) / CONFIG['SCALE'])
        top_y = int((y - CONFIG['SCALE']) / CONFIG['SCALE'])
        bottom_x = int((x + w + CONFIG['SCALE']) / CONFIG['SCALE'])
        bottom_y = int((y + h + CONFIG['SCALE']) / CONFIG['SCALE'])
        boxes.append(((top_x, top_y), (bottom_x, bottom_y)))
    return np.array(boxes)


def dilate_betti0_regions(initial_image: np.ndarray, boxes: np.ndarray, origin_image: np.ndarray, initial_contours: List[np.ndarray], origin_contours: List[np.ndarray]) -> np.ndarray:
    """Perform dilation on white regions"""
    height, width = initial_image.shape[:2]
    check_contour = len(initial_contours)
    initial_betti0_image = initial_image.copy()
    origin_betti0_image = origin_image.copy()

    cv2.drawContours(initial_betti0_image, initial_contours, -1, CONFIG['WHITE_COLOR'], -1)
    cv2.drawContours(origin_betti0_image, origin_contours, -1, CONFIG['WHITE_COLOR'], -1)

    initial_betti0_image_mask = initial_betti0_image.copy()

    for contour in initial_contours:
        point_x, point_y, point_w, point_h = cv2.boundingRect(contour)

        for box in boxes:
            if (box[0][0] <= point_x or box[1][0] >= (point_x + point_w) or
                box[0][1] <= point_y or box[1][1] >= (point_y + point_h)):

                for y in range(box[0][1], box[1][1] + 1):
                    for x in range(box[0][0], box[1][0] + 1):
                        if not (point_x < x < (point_x + point_w) or point_y < y < (point_y + point_h)):
                            in_other_boxes = any(
                                other_box[0][0] + 1 <= x <= other_box[1][0] - 1 and
                                other_box[0][1] + 1 <= y <= other_box[1][1] - 1
                                for other_box in boxes if other_box.all() != box.all()
                            )

                            if not in_other_boxes:
                                x_original = int(x * CONFIG['SCALE'])
                                y_original = int(y * CONFIG['SCALE'])
                                
                                resampled_pixel_value = bilinear_interpolation(origin_image, x_original, y_original)

                                if np.all(resampled_pixel_value == CONFIG['WHITE_COLOR']):
                                    _, betti0_thresh = preprocess_image(initial_betti0_image)
                                    betti0_contours, _ = find_contours(betti0_thresh)

                                    neighbors = [(y-1, x), (y, x-1), (y, x+1), (y+1, x),
                                                 (y-1, x-1), (y-1, x+1), (y+1, x-1), (y+1, x+1)]
                                    
                                    check_neighbor = []
                                    for neighbor_y, neighbor_x in neighbors:
                                        if 0 <= neighbor_y < height and 0 < neighbor_x < width:
                                            if np.all(initial_betti0_image[neighbor_y, neighbor_x] == CONFIG['WHITE_COLOR']):
                                                for idx, cont in enumerate(betti0_contours):
                                                    check = cv2.pointPolygonTest(cont, (neighbor_x, neighbor_y), False)
                                                    if check >= 0:
                                                        check_neighbor.append(idx)
                                    
                                    if len(set(check_neighbor)) == 1:
                                        if 0 <= y < height and 0 < x < width:               
                                            initial_betti0_image[y, x] = CONFIG['WHITE_COLOR']
                                            _, betti0_thresh_check = preprocess_image(initial_betti0_image)
                                            betti0_check_contours, _ = find_contours(betti0_thresh_check)
                                            if len(betti0_check_contours) != check_contour:
                                                initial_betti0_image[y, x] = CONFIG['BLACK_COLOR']

        _, initial_betti0_thresh = preprocess_image(initial_betti0_image)
        initial_betti0_contours, _ = find_contours(initial_betti0_thresh)

        cv2.drawContours(initial_betti0_image, initial_betti0_contours, -1, (255, 255, 255), -1)

        betti0_result = cv2.bitwise_and(initial_betti0_image, cv2.bitwise_not(initial_betti0_image_mask))
        betti0_result = cv2.cvtColor(betti0_result, cv2.COLOR_BGR2GRAY)
        _, betti0_result = cv2.threshold(betti0_result, 1, 255, cv2.THRESH_BINARY)

        white_pixels = cv2.findNonZero(betti0_result)
        if white_pixels is not None:
            for pixel in white_pixels:
                x, y = pixel[0]
                # 将原始图像中相应位置的像素值设为白色
                initial_image[y, x] = [255, 255, 255]

    return initial_image

def white_dilation(origin_image: np.ndarray, initial_image: np.ndarray) -> np.ndarray:
    """
    Perform white region dilation operation.
    
    Args:
        origin_image (np.ndarray): Original input image.
        initial_image (np.ndarray): Initially processed image.
    
    Returns:
        np.ndarray: Dilated image.
    
    Raises:
        ValueError: If input images are None.
    """
    if origin_image is None or initial_image is None:
        raise ValueError("Input images cannot be None")
    
    # Preprocess original image
    _, origin_thresh = preprocess_image(origin_image)
    origin_contours, _= find_contours(origin_thresh)
    
    # Preprocess initial image
    _, initial_thresh = preprocess_image(initial_image)
    initial_contours, _ = find_contours(initial_thresh)
    
    # Process contours
    boxes = process_betti0_contours(origin_contours)
    
    # Perform white region dilation
    initial_betti0_image = dilate_betti0_regions(initial_image, boxes, origin_image, initial_contours, origin_contours)
    
    
    return initial_betti0_image


def process_image(input_path, output_path):
    origin_img = cv2.imread(input_path)
    initial_image = initial(input_path)
    result_img = white_dilation(origin_img, initial_image)
    cv2.imwrite(output_path, result_img)


def process_directory(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for filename in os.listdir(input_dir):
        if filename.lower().endswith(('.png', '.jpg', '.jpeg')):
            input_path = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, f"processed_{filename}")
            process_image(input_path, output_path)


if __name__ == "__main__":

    

    if len(sys.argv) < 2:
        print("Usage: python script.py <input_path> [output_path]")
        sys.exit(1)

    input_path = sys.argv[1]
    
    if os.path.isfile(input_path):
        output_path = sys.argv[2] if len(sys.argv) > 2 else "output_image.jpg"
        process_image(input_path, output_path)
    elif os.path.isdir(input_path):
        output_dir = sys.argv[2] if len(sys.argv) > 2 else "output"
        process_directory(input_path, output_dir)
    else:
        print("Invalid input path")
        sys.exit(1)

    print("Processing complete.")