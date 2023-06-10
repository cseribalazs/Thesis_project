import geojson
import pandas as pd
from shapely.geometry import Polygon, Point, mapping
import json
import geopy
import geopy.distance
import rasterio as rio
from pathlib import Path
import os
from rasterio import mask, plot
from rasterio.warp import calculate_default_transform, reproject, Resampling
import matplotlib.pyplot as plt
import numpy as np

DISTANCE = 2
longitude_coordinates = 0
latitude_coordinates = 0


def find_file_by_point(lon, lat, dictionary):
    for key, val in dictionary.items():
        if Polygon(val).contains(Point(lat, lon)):
            return key
    return None


def polygon_definition(lon, lat, distance):
    polygon_list = []
    start = geopy.Point(lat, lon)
    r = geopy.distance.distance(kilometers=distance / 2)
    d = geopy.distance.distance(kilometers=distance)
    new_point = r.destination(point=start, bearing=0)
    polygon_1 = r.destination(point=new_point, bearing=90)
    coordinates = polygon_1[1], polygon_1[0]
    polygon_list.append(tuple(coordinates))
    for i in range(3):
        polygon_point = d.destination(point=polygon_1, bearing=(90 * (i + 2)))
        polygon_1 = polygon_point
        coordinates = polygon_1[1], polygon_1[0]
        polygon_list.append(tuple(coordinates))
    polygon_list.append(polygon_list[0])
    return polygon_list


with open('recog.json', 'r') as json_file:
    data = json.load(json_file)

file_name = find_file_by_point(longitude_coordinates, latitude_coordinates, data)
geom_list = polygon_definition(longitude_coordinates, latitude_coordinates, DISTANCE)
file_path = os.path.join(Path.cwd().parent, 'raw_DEM\\{}'.format(file_name))

geometry = geojson.Polygon([geom_list])

with rio.open("{}".format(file_path)) as src:
    out_image, out_transform = rio.mask.mask(src, [geometry], crop=True, all_touched=True)
    out_meta = src.meta.copy()
    out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform})

    elevation = src.read(1)
    transform = src.transform

    x = transform[2] + np.arange(src.width) * transform[0]
    y = transform[5] + np.arange(src.height) * transform[4]
    xx, yy = np.meshgrid(x, y)


    lat_min, lat_max = np.percentile(yy.flatten(), 40), np.percentile(yy.flatten(), 60)
    lon_min, lon_max = np.percentile(xx.flatten(), 40), np.percentile(xx.flatten(), 60)
    highlight_mask = np.logical_and(yy >= lat_min, yy <= lat_max) & np.logical_and(xx >= lon_min, xx <= lon_max)
    slope_x = np.gradient(elevation, axis=1)
    slope_y = np.gradient(elevation, axis=0)
    slope_magnitude = np.sqrt(slope_x ** 2 + slope_y ** 2)


    cropped_slope_magnitude = slope_magnitude[highlight_mask]
    slope_std = np.std(cropped_slope_magnitude)

    average_slope = np.mean(cropped_slope_magnitude)
    src.close()



