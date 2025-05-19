import json
import logging
import math
import os
import traceback
import warnings
from datetime import datetime
from math import cos, sin

import geopandas as gpd
import matplotlib.pyplot as plt
import overpy
import pandas as pd
import pyproj
from pyproj import Transformer
from shapely.geometry import LineString, Point, Polygon, box
from shapely.ops import nearest_points, transform

warnings.filterwarnings("ignore")


# this function transforms the polygons artistically. ultimately I want to introduce a way to add multiple transformation parameters to create a map, but for now there is only one optional transfornation and that's a swirl which rotates points around a center. the strength of the rotation depends on the distance to a virtual circle around the centre. note: this function transforms points. for lines to be transformed in a non geometric way, they first need to be "densified", i.e. points need to be added with a certain resolution   for them to then represent non straight line
"""def densifyPolygons(row, maxSegmentLength):  # segmentize the geometry in ech row
    return row["geometry"].segmentize(maxSegmentLength)
"""


def initialize_crs(card):
    card["crs"] = createCRS(card["center"])  # create the card CRS
    card["center_DMS"] = dec2dms(
        card["center"][0], card["center"][1]
    )  # create the card center in degree, minutes, seconds

    # set the wgs to crs and vice versa transformations
    t_crs2WGS84 = Transformer.from_crs(
        card["crs"], "epsg:4326", always_xy=True
    ).transform

    t_WGS842crs = Transformer.from_crs(
        "epsg:4326", card["crs"], always_xy=True
    ).transform
    # set the projected center (with the used crs this will always be 0,0 - but if the card crs changes it wont
    card["center_projected"] = (
        transform(t_WGS842crs, Point(card["center"][1], card["center"][0])).y,
        transform(t_WGS842crs, Point(card["center"][1], card["center"][0])).x,
    )  # transform(t_WGS2CRS,box)
    return card, t_crs2WGS84


# this actually creates a custom crs or more specifically a Mercator CRS around the centre of the map
def createCRS(cardcenter):
    # create a custom crs based on UTM for the card centre
    crs_string = (
        f"+lon_0={cardcenter[1]} +ellps=WGS84 +y_0=0 +no_defs=True +proj=tmerc +x_0=0 "
        f"+units=m  +lat_0={cardcenter[0]} +k=1 +towgs84=0,0,0,0,0,0,0"
    )
    crs = pyproj.CRS.from_string(crs_string)
    return crs


def calculate_card_dimensions(card, ax, fig, cm, cwidthmm):
    """Calculate dimensions and scale for the card."""
    dummyPoly = Polygon(
        [
            [card["width"] / 2, card["height"] / 2],
            [card["width"] / 2, card["height"] / -2],
            [card["width"] / -2, card["height"] / -2],
            [card["width"] / -2, card["height"] / 2],
        ]
    )

    dummyGDF = gpd.GeoDataFrame(index=[0], geometry=[dummyPoly], crs=card["crs"])

    ax[0, 0] = dummyGDF.plot(ax=ax[0, 0], color="None", edgecolor="None")
    axisLimits = ax[0, 0].axes.get_xlim()
    axisincm = abs(axisLimits[1] - axisLimits[0]) * 100
    fig.canvas.draw()  # Force a draw to update transformations
    axisinpixel = (
        ax[0, 0].transData.transform((axisLimits[1], 0))[0]
        - ax[0, 0].transData.transform((axisLimits[0], 0))[0]
    )
    card["scale"] = axisincm / ((axisinpixel) / fig.dpi / cm)

    # =max((card['width'] * 10000)/(ax[0, 0].bbox.width /cm),(card['height'] * 10000)/(ax[0, 0].bbox.height /cm))
    card["scaleText"] = ("Scale 1:{:,}".format(round(card["scale"]))).replace(",", " ")
    cwidth = cwidthmm * card["scale"] / 1000  # connector width is set to mm on paper
    card["layers"][2]["graphbuffer"] = max(
        card["layers"][2]["graphbuffer_mm"] * card["scale"] / 1000, 2
    )  # streeets are min graphbuffer_mm or if wider, 1m buffered on map
    card["layers"][4]["graphbuffer"] = max(
        card["layers"][4]["graphbuffer_mm"] * card["scale"] / 1000, 2
    )  # rail  are min graphbuffer_mm or if wider, 1m buffered on map
    card["frame_width"] = card["frame_width_mm"] * card["scale"] / 1000
    card["hole_radius_for_libro"] = (
        card["hole_radius_for_libro_in_mm"] * card["scale"] / 1000
    )
    card["hole_radius_for_globe"] = (
        card["hole_radius_for_globe_in_mm"] * card["scale"] / 1000
    )
    card["hole_radius_for_frame"] = (
        card["hole_radius_for_frame_in_mm"] * card["scale"] / 1000
    )
    card["visible_dimension_x"] = (
        card["width"] - card["frame_width_mm"] * card["scale"] / 1000
    )
    card["visible_dimension_y"] = (
        card["height"] - card["frame_width_mm"] * card["scale"] / 1000
    )

    return card, cwidth


# this function does not give the correct results with coordinates in NZ - probably a bug with dealing with negative latitude?
def dec2dms(latitude, longitude):
    # Convert latitude
    lat_deg = int(latitude)
    lat_min = int((abs(latitude) * 60) % 60)
    lat_sec = (abs(latitude) * 3600) % 60
    lat_direction = "S" if latitude < 0 else "N"
    lat_dms = f"{abs(lat_deg)}\u00b0{lat_min}'{lat_sec:.0f}\"{lat_direction}"

    # Convert longitude
    long_deg = int(longitude)
    long_min = int((abs(longitude) * 60) % 60)
    long_sec = (abs(longitude) * 3600) % 60
    long_direction = "W" if longitude < 0 else "E"
    long_dms = f"{abs(long_deg)}\u00b0{long_min}'{long_sec:.0f}\"{long_direction}"

    return f"{lat_dms} {long_dms}"


def c_linesold(gdf, rects, layer_name, clinesyes, card, cwidth):
    # Function removed as it is not used anywhere in the backend.
    pass


def c_lines(gdf, layer_name, clinesyes, card, cwidth):
    gdf = gdf.dissolve().explode(index_parts=True)
    all_clines = gpd.GeoDataFrame(columns=["geometry"], crs=card["crs"])
    if cwidth == 0 or clinesyes == 0:  # return an empty cliones gdf
        return gdf.dissolve().explode(index_parts=True), all_clines
    outer_hull = card["frame"].convex_hull
    ring_area = 0
    total_frame_area = outer_hull.area[0]
    buffer_increment = (
        -min(card["width"], card["height"]) / 2 / card["number_cline_buffer"]
    )
    total_buffer = buffer_increment
    while (
        ring_area < total_frame_area * 0.98
    ):  # while the ring area is smaller than the whole area
        ring = outer_hull.difference(
            outer_hull.buffer(total_buffer)
        )  # go from out to in with increasiong ring area
        # ring=outer_hull-outer_hull.buffer(total_buffer) # go from out to in with increasiong ring area
        total_buffer += buffer_increment
        ring_area = ring.area[0]
        clines = gpd.GeoDataFrame(columns=["geometry"], crs=card["crs"])
        subgdf = gdf.clip(ring, keep_geom_type=True)
        subgdf = subgdf[subgdf["geometry"] != "POLYGON EMPTY"]
        subgdf = subgdf.explode(index_parts=True)
        print(
            f"{layer_name} increasing buffer from outside- {ring_area / total_frame_area * 100:.2f}%: connecting {len(subgdf)} polygons with ",
            end=" ",
        )
        while (
            subgdf.shape[0] > 1
        ):  # add connector lines to the dataframe, explode and dissolve until therev is only 1 polygon, i.e. everything is connected
            for index, row in subgdf.iterrows():
                point = row.geometry  # take the next "point" (actually any geometry to check against the combined (unary union) of the rest
                multipoint = subgdf.drop(
                    index, axis=0
                ).geometry.unary_union  # drop the checked geometry from the rest to avoid finding the geometry itself as the nearest neighbor
                queried_geom, nearest_geom = nearest_points(point, multipoint)
                # card['frame'].geometry.contains(nearest_geom).any() and card['frame'].geometry.contains(queried_geom).any()
                if nearest_geom != queried_geom:
                    # print("\nnearest:"+str(card['frame'].geometry.contains(nearest_geom).any())+"-queried:"+str(card['frame'].geometry.contains(queried_geom).any()),end=" ")
                    if (
                        card["frame"].geometry.contains(nearest_geom).any()
                        & card["frame"].geometry.contains(queried_geom).any()
                    ):
                        print("***", end="")

                    queried_geom = point.intersection(
                        queried_geom.buffer(2 * cwidth)
                    ).centroid  # move the point inwards to avoid corners in the cline
                    nearest_geom = multipoint.intersection(
                        nearest_geom.buffer(2 * cwidth)
                    ).centroid
                    cline = LineString([nearest_geom, queried_geom]).buffer(
                        cwidth, cap_style=3
                    )
                    clines = pd.concat(
                        [
                            clines,
                            gpd.GeoDataFrame([{"geometry": cline}], crs=card["crs"]),
                        ],
                        ignore_index=True,
                    )
                    # clines = clines.append(gpd.GeoDataFrame([{'geometry': cline}], crs=card['crs']),ignore_index=True)
                # subgdf = subgdf.append(clines).dissolve().explode(index_parts=True)
            subgdf = (
                pd.concat([subgdf, clines], ignore_index=True)
                .dissolve()
                .explode(index_parts=True)
            )

        # all_clines=all_clines.append(clines)
        all_clines = pd.concat([all_clines, clines], ignore_index=True)
        print(f"{len(clines)} connector lines")
        gdf = (
            pd.concat([gdf, all_clines], ignore_index=True)
            .dissolve()
            .explode(index_parts=True)
        )
        # gdf=gdf.append(all_clines).dissolve().explode(index_parts=True)
        gdf = gdf.dissolve().explode(index_parts=True)
    return gdf, all_clines


def c_linesGrid(gdf, grids, layer_name, clinesyes, card, cwidth):
    all_clines = gpd.GeoDataFrame(columns=["geometry"], crs=card["crs"])
    if cwidth == 0 or clinesyes == 0:  # return an empty cliones gdf
        return gdf.dissolve().explode(index_parts=True), all_clines
    ### duplicate  the algo to go through all the grid polygons to add to the clines and eventually to the gdf
    clines = []
    for grid in grids:
        clines = gpd.GeoDataFrame(columns=["geometry"], crs=card["crs"])
        subgdf = gdf.clip(grid, keep_geom_type=True)
        subgdf = subgdf[subgdf["geometry"] != "POLYGON EMPTY"]
        subgdf = subgdf.explode(index_parts=True)
        print(
            f"{layer_name}-grid {grids.index(grid)}: connecting {len(subgdf)} polygons. Rectangle area {int(grid.area)}m2 with ",
            end=" ",
        )
        while (
            subgdf.shape[0] > 1
        ):  # add connector lines to the dataframe, explode and dissolve until there is only 1 polygon, i.e. everything is connected
            for index, row in subgdf.iterrows():
                point = row.geometry  # take the next "point" (actually any geometry to check against the combined (unary union) of the rest
                multipoint = subgdf.drop(
                    index, axis=0
                ).geometry.unary_union  # drop the checked geometry from the rest to avoid finding the geometry itself as the nearest neighbor
                queried_geom, nearest_geom = nearest_points(point, multipoint)
                if nearest_geom != queried_geom:
                    queried_geom = point.intersection(
                        queried_geom.buffer(2 * cwidth)
                    ).centroid
                    nearest_geom = multipoint.intersection(
                        nearest_geom.buffer(2 * cwidth)
                    ).centroid
                    cline = LineString([nearest_geom, queried_geom]).buffer(
                        cwidth, cap_style=3
                    )
                    clines = pd.concat(
                        [
                            clines,
                            gpd.GeoDataFrame([{"geometry": cline}], crs=card["crs"]),
                        ],
                        ignore_index=True,
                    )
                    # clines = clines.append(gpd.GeoDataFrame([{'geometry': cline}], crs=card['crs']), ignore_index=True)
            subgdf = (
                pd.concat([subgdf, clines], ignore_index=True)
                .dissolve()
                .explode(index_parts=True)
            )
            # subgdf = subgdf.append(clines).dissolve().explode(index_parts=True)
        all_clines = pd.concat([all_clines, clines], ignore_index=True)
        # all_clines = all_clines.append(clines)
        print(f"{len(clines)} connector lines")
        gdf = (
            pd.concat([gdf, all_clines], ignore_index=True)
            .dissolve()
            .explode(index_parts=True)
        )
        # gdf = gdf.append(all_clines).dissolve().explode(index_parts=True)
        gdf = gdf.dissolve().explode(index_parts=True)
    return gdf, all_clines


def get_rectangles(center, card_height, card_width, frame_width, gridx, gridy, card):
    # TODO: set all the variable accessible to card["key"]
    rectangles = []  # these are rectangles from out to in
    grids = []  # this is a grid of rectangles
    # create the outer frame rectangle
    # with the card crs centered around the center of the card, the center here is 0, but not for others crs. i should use the projected card center instead.
    rectangles.append(
        Polygon(
            [
                [center[1] + card_width / 2, center[0] + card_height / 2],
                [center[1] - card_width / 2, center[0] + card_height / 2],
                [center[1] - card_width / 2, center[0] - card_height / 2],
                [center[1] + card_width / 2, center[0] - card_height / 2],
            ]
        )
    )
    # iterate from the outer frame rectangle to the middle with framewidth iterations
    for i in range(
        1, int(((min(card_height, card_width) + frame_width * 2) / 2) // frame_width)
    ):
        rectangles.append(rectangles[0].buffer(-frame_width * i, join_style=2))

    # secondary grid of rectangles for more clines
    x = range(gridx)[1:]
    y = range(gridy)[1:]
    for gridx in x:
        for gridy in y:
            for i in range(gridx):
                for ii in range(gridy):
                    grids.append(
                        Point(
                            center[1] - card_width / 2 + i / gridx * card_width,
                            center[0] - card_height / 2 + ii / gridy * card_height,
                        ).buffer(max(card_width / gridx, card_height / gridy))
                    )  # this makeds it a square, but it should be a rectangle depending on the gridx and grid y and card dimanesion...

    match card["binding"]:
        case "libro":  # building a polygon with a double thick left frame for the binding
            inner_bframe = Polygon(
                [
                    [
                        center[1] + card_width / 2 - frame_width,
                        center[0] + card_height / 2 - frame_width,
                    ],
                    [
                        center[1] - card_width / 2 + 2 * frame_width,
                        center[0] + card_height / 2 - frame_width,
                    ],
                    [
                        center[1] - card_width / 2 + 2 * frame_width,
                        center[0] - card_height / 2 + frame_width,
                    ],
                    [
                        center[1] + card_width / 2 - frame_width,
                        center[0] - card_height / 2 + frame_width,
                    ],
                ]
            )
            bframe = rectangles[0] - inner_bframe

            ### create binding block with holes
            number_of_squares = int(card["height"] // card["frame_width"] - 1)
            square_edge = card["height"] / number_of_squares
            # bblock_outer=Polygon([[center[1] - card_width / 2, center[0] - card_height / 2],
            #                     [center[1] - card_width / 2 + (2*card['frame_width']), center[0] - card_height / 2],
            #                     [center[1] - card_width / 2 + (2*card['frame_width']), center[0] + card_height / 2],
            #                     [center[1] - card_width / 2, center[0] + card_height / 2]])
            # bblock=bblock_outer
            for i in range(number_of_squares):
                if i > 0:
                    bframe = bframe - Point(
                        center[1] - card_width / 2 + square_edge,
                        center[0] - card_height / 2 + i * square_edge,
                    ).buffer(card["hole_radius_for_libro"])
                bframe = bframe - Point(
                    center[1] - card_width / 2 + square_edge / 2,
                    center[0] - card_height / 2 + i * square_edge + square_edge / 2,
                ).buffer(card["hole_radius_for_libro"])
            bframegdf = gpd.GeoDataFrame(index=[0], geometry=[bframe], crs=card["crs"])
        case "carta":
            inner_bframe = Polygon(
                [
                    [
                        center[1] + card_width / 2 - frame_width,
                        center[0] + card_height / 2 - frame_width,
                    ],
                    [
                        center[1] - card_width / 2 + frame_width,
                        center[0] + card_height / 2 - frame_width,
                    ],
                    [
                        center[1] - card_width / 2 + frame_width,
                        center[0] - card_height / 2 + frame_width,
                    ],
                    [
                        center[1] + card_width / 2 - frame_width,
                        center[0] - card_height / 2 + frame_width,
                    ],
                ]
            )
            bframe = rectangles[0] - inner_bframe

            x1 = -card_width / 2 + frame_width
            x2 = -frame_width / 2
            x3 = frame_width / 2
            x4 = card_width / 2 - frame_width
            y1 = -card_height / 2 + frame_width * 2
            y2 = card_height / 2 - frame_width

            bframe = rectangles[0] - box(x1, y1, x2, y2) - box(x3, y1, x4, y2)
            bframegdf = gpd.GeoDataFrame(index=[0], geometry=[bframe], crs=card["crs"])
        case "frame":  # building a polygon with for the binding
            inner_bframe = Polygon(
                [
                    [
                        center[1] + card_width / 2 - frame_width,
                        center[0] + card_height / 2 - frame_width,
                    ],
                    [
                        center[1] - card_width / 2 + frame_width,
                        center[0] + card_height / 2 - frame_width,
                    ],
                    [
                        center[1] - card_width / 2 + frame_width,
                        center[0] - card_height / 2 + frame_width,
                    ],
                    [
                        center[1] + card_width / 2 - frame_width,
                        center[0] - card_height / 2 + frame_width,
                    ],
                ]
            )
            bframe = rectangles[0] - inner_bframe

            ### create holes
            num_of_holes = card["number_of_holes_for_frame"]
            # sizeOfHolesmm=card['hole_radius_for_frame'] #in card mm radius
            # sizeOfHoles=sizeOfHolesmm*card['scale']/1000
            sizeOfHoles = card["hole_radius_for_frame"]
            for i in range(1, num_of_holes):  # number of holes
                bframe = bframe - Point(
                    center[1] - card_width / 2 + frame_width / 4,
                    center[0]
                    - card_height / 2
                    + frame_width / 2
                    + i * (card_height - frame_width) / num_of_holes,
                ).buffer(sizeOfHoles)
                bframe = bframe - Point(
                    center[1] + card_width / 2 - frame_width / 4,
                    center[0]
                    - card_height / 2
                    + frame_width / 2
                    + i * (card_height - frame_width) / num_of_holes,
                ).buffer(sizeOfHoles)
                bframe = bframe - Point(
                    center[1]
                    - card_width / 2
                    + frame_width / 2
                    + i * (card_width - frame_width) / num_of_holes,
                    center[0] - card_height / 2 + frame_width / 4,
                ).buffer(sizeOfHoles)
                bframe = bframe - Point(
                    center[1]
                    - card_width / 2
                    + frame_width / 2
                    + i * (card_width - frame_width) / num_of_holes,
                    center[0] + card_height / 2 - frame_width / 4,
                ).buffer(sizeOfHoles)
                bframe = bframe - Point(
                    center[1] - card_width / 2 + frame_width * (3 / 4),
                    center[0]
                    - card_height / 2
                    + frame_width / 2
                    + i * (card_height - frame_width) / num_of_holes,
                ).buffer(sizeOfHoles)
                bframe = bframe - Point(
                    center[1] + card_width / 2 - frame_width * (3 / 4),
                    center[0]
                    - card_height / 2
                    + frame_width / 2
                    + i * (card_height - frame_width) / num_of_holes,
                ).buffer(sizeOfHoles)
                bframe = bframe - Point(
                    center[1]
                    - card_width / 2
                    + frame_width / 2
                    + i * (card_width - frame_width) / num_of_holes,
                    center[0] - card_height / 2 + frame_width * (3 / 4),
                ).buffer(sizeOfHoles)
                bframe = bframe - Point(
                    center[1]
                    - card_width / 2
                    + frame_width / 2
                    + i * (card_width - frame_width) / num_of_holes,
                    center[0] + card_height / 2 - frame_width * (3 / 4),
                ).buffer(sizeOfHoles)
            # cornerholes
            bframe = bframe - Point(
                center[1] - card_width / 2 + frame_width / 4 * 1,
                center[0] - card_height / 2 + frame_width / 4 * 3,
            ).buffer(sizeOfHoles)
            bframe = bframe - Point(
                center[1] - card_width / 2 + frame_width / 4 * 3,
                center[0] - card_height / 2 + frame_width / 4 * 1,
            ).buffer(sizeOfHoles)
            bframe = bframe - Point(
                center[1] - card_width / 2 + frame_width / 4 * 3,
                center[0] - card_height / 2 + frame_width / 4 * 3,
            ).buffer(sizeOfHoles)

            bframe = bframe - Point(
                center[1] - card_width / 2 + frame_width / 4 * 1,
                center[0] + card_height / 2 - frame_width / 4 * 3,
            ).buffer(sizeOfHoles)
            bframe = bframe - Point(
                center[1] - card_width / 2 + frame_width / 4 * 3,
                center[0] + card_height / 2 - frame_width / 4 * 1,
            ).buffer(sizeOfHoles)
            bframe = bframe - Point(
                center[1] - card_width / 2 + frame_width / 4 * 3,
                center[0] + card_height / 2 - frame_width / 4 * 3,
            ).buffer(sizeOfHoles)

            bframe = bframe - Point(
                center[1] + card_width / 2 - frame_width / 4 * 1,
                center[0] + card_height / 2 - frame_width / 4 * 3,
            ).buffer(sizeOfHoles)
            bframe = bframe - Point(
                center[1] + card_width / 2 - frame_width / 4 * 3,
                center[0] + card_height / 2 - frame_width / 4 * 1,
            ).buffer(sizeOfHoles)
            bframe = bframe - Point(
                center[1] + card_width / 2 - frame_width / 4 * 3,
                center[0] + card_height / 2 - frame_width / 4 * 3,
            ).buffer(sizeOfHoles)

            bframe = bframe - Point(
                center[1] + card_width / 2 - frame_width / 4 * 1,
                center[0] - card_height / 2 + frame_width / 4 * 3,
            ).buffer(sizeOfHoles)
            bframe = bframe - Point(
                center[1] + card_width / 2 - frame_width / 4 * 3,
                center[0] - card_height / 2 + frame_width / 4 * 3,
            ).buffer(sizeOfHoles)
            bframe = bframe - Point(
                center[1] + card_width / 2 - frame_width / 4 * 3,
                center[0] - card_height / 2 + frame_width / 4 * 1,
            ).buffer(sizeOfHoles)

            bframegdf = gpd.GeoDataFrame(index=[0], geometry=[bframe], crs=card["crs"])
        ### experimental
        case _:  # removed 'globe' unused variable
            sizeOfHoles = card["hole_radius_for_globe"]
            circle = Point(0, 0).buffer(card_width / 2) - Point(0, 0).buffer(
                card_width / 2 - frame_width
            )
            inner_bframe = Point(0, 0).buffer(card_width / 2 - frame_width * 1.1)
            num_of_holes = card["number_of_holes_for_globe"]
            for i in range(0, num_of_holes):
                circle = circle - Point(
                    center[1]
                    + sin(i * 6.28318543 / num_of_holes)
                    * (card_width / 2 - frame_width / 2),
                    center[0]
                    + cos(i * 6.28318543 / num_of_holes)
                    * (card_height / 2 - frame_width / 2),
                ).buffer(sizeOfHoles)
            bframegdf = gpd.GeoDataFrame(index=[0], geometry=[circle], crs=card["crs"])

    return inner_bframe, rectangles, bframegdf, grids


def simplify(layer, gdfname, card_area):
    # global card_area
    gdf = layer[gdfname]
    if layer["type"] == "graph":
        gdf["geometry"] = gdf["geometry"].buffer(
            layer["graphbuffer"]
        )  # converting line streets into polygons. should work with a streetwidth dictionary
    else:
        gdf["geometry"] = gdf["geometry"].buffer(layer["emph_buffer"], join_style=1)
        gdf["geometry"] = gdf["geometry"].buffer(layer["smooth"], join_style=1)

        gdf = gdf.dissolve()
        gdf["geometry"] = gdf["geometry"].unary_union
        gdf = gdf.explode(
            index_parts=True
        )  # not exactly sure what this does but i think it makes multiploygons one geometry and possibly even makes joins all polygons together

        gdf["geometry"] = gdf["geometry"].simplify(
            layer["simplify_tolerance"]
        )  # straighten the edges
        gdf["geometry"] = gdf["geometry"].buffer(
            -layer["smooth"], join_style=1
        )  # reverse the applied buffer

        # reduce buffer and epand again to remove small thin pieces
        gdf["geometry"] = gdf["geometry"].buffer(
            -layer["smooth"], join_style=1
        )  # reverse the applied buffer
        gdf["geometry"] = gdf["geometry"].buffer(layer["smooth"], join_style=1)

        gdf["area"] = gdf["geometry"].area  # add a column with area
        gdf = gdf[
            gdf["area"] > (card_area * layer["pat"])
        ]  # filter out polygons that dont meet the treshold criteria of pat% or total area
        gdf = gdf.drop("area", axis=1)
    return gdf


def generate_map(lat, lon, width, height, formats, style, location, layers):
    local_path = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(local_path, "../card.json")) as f:
        card_template = json.load(f)
    with open(os.path.join(local_path, "../lightburncolor.json")) as f:
        color_map = json.load(f)
    output_dir = os.path.abspath(os.path.join(local_path, "../output"))
    os.makedirs(output_dir, exist_ok=True)
    dt_string = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    card = card_template.copy()
    card["center"] = [lat, lon]
    card["width"] = width
    card["height"] = height
    card["location"] = location
    card["layers"] = [
        layer_obj
        for layer_obj in card_template["layers"]
        if layer_obj["name"] in layers
    ]
    cm = 1 / 2.54
    from pyproj import Transformer
    from shapely.geometry import Point, box
    from shapely.ops import transform

    def createCRS(cardcenter):
        crs_string = (
            f"+lon_0={cardcenter[1]} +ellps=WGS84 +y_0=0 +no_defs=True +proj=tmerc +x_0=0 "
            f"+units=m  +lat_0={cardcenter[0]} +k=1 +towgs84=0,0,0,0,0,0,0"
        )
        import pyproj

        return pyproj.CRS.from_string(crs_string)

    def initialize_crs(card):
        card["crs"] = createCRS(card["center"])
        t_crs2WGS84 = Transformer.from_crs(
            card["crs"], "epsg:4326", always_xy=True
        ).transform
        t_WGS842crs = Transformer.from_crs(
            "epsg:4326", card["crs"], always_xy=True
        ).transform
        card["center_projected"] = (
            transform(t_WGS842crs, Point(card["center"][1], card["center"][0])).y,
            transform(t_WGS842crs, Point(card["center"][1], card["center"][0])).x,
        )
        return card, t_crs2WGS84

    card, t_crs2WGS84 = initialize_crs(card)
    fig, ax = plt.subplots(
        figsize=(card["paper_width"] * cm, card["paper_height"] * cm)
    )
    ax.axis("off")
    fig.set_facecolor("None")
    fig.set_edgecolor("None")
    fig.set_frameon(True)

    # Calculate bbox as shapely box
    # 1 deg latitude ≈ 111.32 km; 1 deg longitude ≈ 111.32 * cos(lat) km

    lat_km = 111.32
    lon_km = 111.32 * math.cos(math.radians(card["center"][0]))
    half_height_deg = card["height"] / 2 / lat_km
    half_width_deg = card["width"] / 2 / lon_km
    bbox = box(
        card["center"][1] - half_width_deg,
        card["center"][0] - half_height_deg,
        card["center"][1] + half_width_deg,
        card["center"][0] + half_height_deg,
    )
    minx, miny, maxx, maxy = bbox.bounds
    bbox_tuple = (maxy, miny, maxx, minx)  # (north, south, east, west)
    max_layers = 6
    plotted_any = False
    all_clipped_gdfs = []
    for layer in card["layers"][:max_layers]:
        logging.info(f"Processing layer: {layer['name']}")
        logging.info(f"BBOX: {bbox_tuple}")
        logging.info(f"OSM tags: {layer['osm_tags']}")
        try:
            gdf = fetch_osm_layer_overpy(layer, bbox_tuple)
            if gdf.empty:
                logging.warning(
                    f"Layer {layer['name']} is empty for bbox {bbox_tuple} and tags {layer['osm_tags']}"
                )
                continue
            gdf = gdf.to_crs(card["crs"])
            # Clip layer to bounding box (ensure bbox is in same CRS)
            bbox_geom = box(minx, miny, maxx, maxy)
            import pyproj
            from shapely.ops import transform

            print(f"Layer {layer['name']} features before clip: {len(gdf)}")
            # Transform bbox_geom to gdf CRS if needed
            if gdf.crs != "EPSG:4326":
                project = pyproj.Transformer.from_crs(
                    "EPSG:4326", gdf.crs, always_xy=True
                ).transform
                bbox_geom_proj = transform(project, bbox_geom)
            else:
                bbox_geom_proj = bbox_geom
            gdf_clipped = gdf.clip(bbox_geom_proj)
            # Ensure CRS is set
            if gdf_clipped.crs is None or str(gdf_clipped.crs) != str(card["crs"]):
                gdf_clipped = gdf_clipped.set_crs(card["crs"], allow_override=True)
            print(f"Layer {layer['name']} features after clip: {len(gdf_clipped)}")
            all_clipped_gdfs.append(gdf_clipped)
            # Eksportiraj svaki sloj kao GeoJSON
            # geojson_path = os.path.join(
            #     output_dir, f"{dt_string}_{location}_{layer['name']}.geojson"
            # )
            # gdf_clipped.to_file(geojson_path, driver="GeoJSON")
            logging.info(f"Layer {layer['name']} features: {len(gdf)}")
            color_hex = (
                layer["lcolor"][1]
                if len(layer["lcolor"]) > 1
                else color_map.get(layer["lcolor"][0], "#000000")
            )
            edgecolor_hex = color_map.get(layer["lcolor"][0], "#000000")
            if not gdf_clipped.empty:
                gdf_clipped.plot(
                    ax=ax, edgecolor=edgecolor_hex, facecolor="none", linewidth=0.01
                )
            else:
                logging.warning(
                    f"Layer {layer['name']} is empty after clipping, skipping plot."
                )
                continue
            plotted_any = True
        except Exception as e:
            logging.error(f"Error processing layer {layer['name']}: {e}")
            logging.error(traceback.format_exc())
    if not plotted_any:
        ax.text(
            0.5,
            0.5,
            "No map data found for this area/layers.",
            ha="center",
            va="center",
            fontsize=18,
            color="gray",
            transform=ax.transAxes,
        )
    # --- Fit/crop to envelope of all objects with margin ---
    if all_clipped_gdfs:
        all_gdf = gpd.GeoDataFrame(
            pd.concat(all_clipped_gdfs, ignore_index=True), crs=card["crs"]
        )
        if not all_gdf.empty:
            minx_env, miny_env, maxx_env, maxy_env = all_gdf.total_bounds
            # Add margin (5%)
            margin_x = (maxx_env - minx_env) * 0.05
            margin_y = (maxy_env - miny_env) * 0.05
            # Check for finite and positive bounds
            if (
                all(
                    map(
                        lambda v: pd.notnull(v) and abs(v) != float("inf"),
                        [minx_env, miny_env, maxx_env, maxy_env],
                    )
                )
                and (maxx_env - minx_env) > 0
                and (maxy_env - miny_env) > 0
            ):
                ax.set_xlim(minx_env - margin_x, maxx_env + margin_x)
                ax.set_ylim(miny_env - margin_y, maxy_env + margin_y)
    # Set aspect to auto to avoid aspect errors
    ax.set_aspect("auto")
    # --- Set PNG output to 3000x3000 px ---
    try:
        svg_laser_path = os.path.join(output_dir, f"{dt_string}_{location}_laser.svg")
        svg_preview_path = os.path.join(
            output_dir, f"{dt_string}_{location}_preview.svg"
        )
        png_path = os.path.join(output_dir, f"{dt_string}_{location}.png")
        # SVG (Laser Ready)
        if "svg_laser" in formats:
            fig.savefig(svg_laser_path, format="svg", bbox_inches="tight", pad_inches=0)
        # SVG (Preview)
        if "svg_preview" in formats:
            fig_svg, ax_svg = plt.subplots(
                figsize=(card["paper_width"] * cm, card["paper_height"] * cm)
            )
            ax_svg.axis("off")
            fig_svg.patch.set_alpha(0.0)
            ax_svg.set_facecolor((0, 0, 0, 0))
            for idx, layer in enumerate(card["layers"][:max_layers]):
                color_hex = (
                    layer["lcolor"][1]
                    if len(layer["lcolor"]) > 1
                    else color_map.get(layer["lcolor"][0], "#000000")
                )
                edgecolor_hex = color_hex
                if idx < len(all_clipped_gdfs):
                    gdf_clipped = all_clipped_gdfs[idx]
                    if not gdf_clipped.empty:
                        gdf_clipped.plot(
                            ax=ax_svg,
                            edgecolor=edgecolor_hex,
                            facecolor=color_hex,
                            linewidth=1,
                        )
                    else:
                        logging.warning(
                            f"[SVG Preview] Layer {layer['name']} is empty after clipping, skipping plot."
                        )
                        continue
            if all_clipped_gdfs:
                all_gdf = gpd.GeoDataFrame(
                    pd.concat(all_clipped_gdfs, ignore_index=True), crs=card["crs"]
                )
                if not all_gdf.empty:
                    minx_env, miny_env, maxx_env, maxy_env = all_gdf.total_bounds
                    margin_x = (maxx_env - minx_env) * 0.05
                    margin_y = (maxy_env - miny_env) * 0.05
                    ax_svg.set_xlim(minx_env - margin_x, maxx_env + margin_x)
                    ax_svg.set_ylim(miny_env - margin_y, maxy_env + margin_y)
            fig_svg.savefig(
                svg_preview_path,
                format="svg",
                bbox_inches="tight",
                pad_inches=0,
                transparent=True,
            )
            plt.close(fig_svg)
        # PNG (Preview)
        if "png" in formats:
            fig_png, ax_png = plt.subplots(
                figsize=(card["paper_width"] * cm, card["paper_height"] * cm)
            )
            ax_png.axis("off")
            fig_png.patch.set_alpha(0.0)
            ax_png.set_facecolor((0, 0, 0, 0))
            plotted_any_png = False
            valid_gdfs_png = []
            for idx, layer in enumerate(card["layers"][:max_layers]):
                color_hex = (
                    layer["lcolor"][1]
                    if len(layer["lcolor"]) > 1
                    else color_map.get(layer["lcolor"][0], "#000000")
                )
                edgecolor_hex = color_hex
                if idx < len(all_clipped_gdfs):
                    gdf_clipped = all_clipped_gdfs[idx]
                    gdf_clipped = gdf_clipped[
                        gdf_clipped.is_valid & ~gdf_clipped.is_empty
                    ]
                    if gdf_clipped.crs is None:
                        gdf_clipped.set_crs(card["crs"], inplace=True)
                    if not gdf_clipped.empty:
                        gdf_clipped.plot(
                            ax=ax_png,
                            edgecolor=edgecolor_hex,
                            facecolor=color_hex,
                            linewidth=1,
                        )
                        plotted_any_png = True
                        valid_gdfs_png.append(gdf_clipped)
                    else:
                        logging.warning(
                            f"[PNG Preview] Layer {layer['name']} is empty or invalid after filtering, skipping plot."
                        )
            if plotted_any_png:
                all_gdf = gpd.GeoDataFrame(
                    pd.concat(valid_gdfs_png, ignore_index=True), crs=card["crs"]
                )
                if not all_gdf.empty:
                    minx_env, miny_env, maxx_env, maxy_env = all_gdf.total_bounds
                    margin_x = (maxx_env - minx_env) * 0.05
                    margin_y = (maxy_env - miny_env) * 0.05
                    ax_png.set_xlim(minx_env - margin_x, maxx_env + margin_x)
                    ax_png.set_ylim(miny_env - margin_y, maxy_env + margin_y)
            else:
                ax_png.text(
                    0.5,
                    0.5,
                    "No map data found for this area/layers.",
                    ha="center",
                    va="center",
                    fontsize=18,
                    color="gray",
                    transform=ax_png.transAxes,
                )
            fig_png.savefig(
                png_path,
                format="png",
                dpi=600,
                bbox_inches="tight",
                pad_inches=0,
                transparent=True,
            )
            plt.close(fig_png)
        # SVG (Preview)
        if "svg_preview" in formats:
            fig_svg, ax_svg = plt.subplots(
                figsize=(card["paper_width"] * cm, card["paper_height"] * cm)
            )
            ax_svg.axis("off")
            fig_svg.patch.set_alpha(0.0)
            ax_svg.set_facecolor((0, 0, 0, 0))
            plotted_any_svg = False
            valid_gdfs_svg = []
            for idx, layer in enumerate(card["layers"][:max_layers]):
                color_hex = (
                    layer["lcolor"][1]
                    if len(layer["lcolor"]) > 1
                    else color_map.get(layer["lcolor"][0], "#000000")
                )
                edgecolor_hex = color_hex
                if idx < len(all_clipped_gdfs):
                    gdf_clipped = all_clipped_gdfs[idx]
                    gdf_clipped = gdf_clipped[
                        gdf_clipped.is_valid & ~gdf_clipped.is_empty
                    ]
                    if gdf_clipped.crs is None:
                        gdf_clipped.set_crs(card["crs"], inplace=True)
                    if not gdf_clipped.empty:
                        gdf_clipped.plot(
                            ax=ax_svg,
                            edgecolor=edgecolor_hex,
                            facecolor=color_hex,
                            linewidth=1,
                        )
                        plotted_any_svg = True
                        valid_gdfs_svg.append(gdf_clipped)
                    else:
                        logging.warning(
                            f"[SVG Preview] Layer {layer['name']} is empty or invalid after filtering, skipping plot."
                        )
            if plotted_any_svg:
                all_gdf = gpd.GeoDataFrame(
                    pd.concat(valid_gdfs_svg, ignore_index=True), crs=card["crs"]
                )
                if not all_gdf.empty:
                    minx_env, miny_env, maxx_env, maxy_env = all_gdf.total_bounds
                    margin_x = (maxx_env - minx_env) * 0.05
                    margin_y = (maxy_env - miny_env) * 0.05
                    ax_svg.set_xlim(minx_env - margin_x, maxx_env + margin_x)
                    ax_svg.set_ylim(miny_env - margin_y, maxy_env + margin_y)
            else:
                ax_svg.text(
                    0.5,
                    0.5,
                    "No map data found for this area/layers.",
                    ha="center",
                    va="center",
                    fontsize=18,
                    color="gray",
                    transform=ax_svg.transAxes,
                )
            fig_svg.savefig(
                svg_preview_path,
                format="svg",
                bbox_inches="tight",
                pad_inches=0,
                transparent=True,
            )
            plt.close(fig_svg)
        # SVG (Laser Ready)
        if "svg_laser" in formats:
            fig_laser, ax_laser = plt.subplots(
                figsize=(card["paper_width"] * cm, card["paper_height"] * cm)
            )
            ax_laser.axis("off")
            fig_laser.patch.set_alpha(0.0)
            ax_laser.set_facecolor((0, 0, 0, 0))
            plotted_any_laser = False
            valid_gdfs_laser = []
            for idx, layer in enumerate(card["layers"][:max_layers]):
                edgecolor_hex = color_map.get(layer["lcolor"][0], "#000000")
                if idx < len(all_clipped_gdfs):
                    gdf_clipped = all_clipped_gdfs[idx]
                    gdf_clipped = gdf_clipped[
                        gdf_clipped.is_valid & ~gdf_clipped.is_empty
                    ]
                    if gdf_clipped.crs is None:
                        gdf_clipped.set_crs(card["crs"], inplace=True)
                    if not gdf_clipped.empty:
                        gdf_clipped.plot(
                            ax=ax_laser,
                            edgecolor=edgecolor_hex,
                            facecolor="none",
                            linewidth=0.01,
                        )
                        plotted_any_laser = True
                        valid_gdfs_laser.append(gdf_clipped)
                    else:
                        logging.warning(
                            f"[SVG Laser] Layer {layer['name']} is empty or invalid after filtering, skipping plot."
                        )
            if plotted_any_laser:
                all_gdf = gpd.GeoDataFrame(
                    pd.concat(valid_gdfs_laser, ignore_index=True), crs=card["crs"]
                )
                if not all_gdf.empty:
                    minx_env, miny_env, maxx_env, maxy_env = all_gdf.total_bounds
                    margin_x = (maxx_env - minx_env) * 0.05
                    margin_y = (maxy_env - miny_env) * 0.05
                    ax_laser.set_xlim(minx_env - margin_x, maxx_env + margin_x)
                    ax_laser.set_ylim(miny_env - margin_y, maxy_env + margin_y)
            else:
                ax_laser.text(
                    0.5,
                    0.5,
                    "No map data found for this area/layers.",
                    ha="center",
                    va="center",
                    fontsize=18,
                    color="gray",
                    transform=ax_laser.transAxes,
                )
            fig_laser.savefig(
                svg_laser_path,
                format="svg",
                bbox_inches="tight",
                pad_inches=0,
            )
            plt.close(fig_laser)
            # --- REMOVE CLIP-PATH FROM LASER SVG ---
            with open(svg_laser_path, "r", encoding="utf-8") as f:
                svg_data = f.read()
            import re

            svg_data = re.sub(r' clip-path="[^"]*"', "", svg_data)
            with open(svg_laser_path, "w", encoding="utf-8") as f:
                f.write(svg_data)
        if "crs" in card:
            del card["crs"]
        plt.close(fig)
    except Exception as e:
        logging.error(f"Error saving output: {e}")
        logging.error(traceback.format_exc())
        raise
    result = {}
    if "svg_laser" in formats:
        result["svgLaserUrl"] = f"/api/output/{os.path.basename(svg_laser_path)}"
    if "svg_preview" in formats:
        result["svgPreviewUrl"] = f"/api/output/{os.path.basename(svg_preview_path)}"
    if "png" in formats:
        result["pngUrl"] = f"/api/output/{os.path.basename(png_path)}"
    return result


def fetch_osm_layer_overpy(layer, bbox):
    """
    Dohvati OSM podatke za zadani layer i bbox koristeći overpy i vrati GeoDataFrame.
    bbox: (north, south, east, west)
    layer: dict s osm_tags i name
    """
    api = overpy.Overpass()
    # Turbo-style QL upiti po sloju
    if layer["name"] == "buildings":
        ql = f"""
        (
          way["building"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});
          relation["building"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});
        );
        (._;>;);
        out body;
        """
    elif layer["name"] == "greens":
        ql = f"""
        (
          way["leisure"~"park|garden|golf_course|recreation_ground"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});
          relation["leisure"~"park|garden|golf_course|recreation_ground"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});
          way["landuse"~"forest|grass"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});
          relation["landuse"~"forest|grass"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});
        );
        (._;>;);
        out body;
        """
    elif layer["name"] == "streets":
        ql = f"""
        (
          way["highway"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});
          relation["highway"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});
        );
        (._;>;);
        out body;
        """
    elif layer["name"] == "blues":
        ql = f"""
        (
          way["waterway"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});
          relation["waterway"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});
          way["natural"="water"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});
          relation["natural"="water"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});
          way["landuse"="reservoir"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});
          relation["landuse"="reservoir"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});
          way["water"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});
          relation["water"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});
        );
        (._;>;);
        out body;
        """
    elif layer["name"] == "rail":
        ql = f"""
        (
          way["railway"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});
          relation["railway"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});
        );
        (._;>;);
        out body;
        """
    else:
        # fallback: stari način
        ql_parts = []
        for key, value in layer["osm_tags"].items():
            if value is True:
                ql_parts.append(
                    f'way["{key}"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]}); relation["{key}"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});'
                )
            elif isinstance(value, list):
                for v in value:
                    ql_parts.append(
                        f'way["{key}"="{v}"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]}); relation["{key}"="{v}"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});'
                    )
            else:
                ql_parts.append(
                    f'way["{key}"="{value}"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]}); relation["{key}"="{value}"]({bbox[1]},{bbox[3]},{bbox[0]},{bbox[2]});'
                )
        ql = f"""
        (
        {chr(10).join(ql_parts)}
        );
        (._;>;);
        out body;
        """
    print(f"[Overpy QL for layer {layer['name']}]:\n{ql}")
    logging.info(f"[Overpy QL for layer {layer['name']}]:\n{ql}")
    try:
        result = api.query(ql)
    except Exception as e:
        print(f"Overpy error for layer {layer['name']}: {e}")
        return gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")
    # Log broj ways i tipove geometrije
    print(
        f"[Overpy] {layer['name']}: {len(result.ways)} ways, {len(result.relations)} relations"
    )
    geom_types = {}
    records = []
    # Ways: Polygon/LineString
    for way in result.ways:
        coords = [(float(node.lon), float(node.lat)) for node in way.nodes]
        if len(coords) < 2:
            continue
        if coords[0] == coords[-1] and len(coords) > 3:
            geom = Polygon(coords)
            geom_types["Polygon"] = geom_types.get("Polygon", 0) + 1
        elif len(coords) >= 2:
            geom = LineString(coords)
            geom_types["LineString"] = geom_types.get("LineString", 0) + 1
        else:
            continue
        # Spremi OSM tagove u properties
        rec = {"geometry": geom}
        rec.update(way.tags)
        records.append(rec)
    # Relations: MultiPolygon
    from shapely.geometry import MultiPolygon

    for rel in result.relations:
        outers = []
        for member in rel.members:
            if member.role == "outer" and hasattr(member, "resolve"):
                way = member.resolve()
                coords = [(float(node.lon), float(node.lat)) for node in way.nodes]
                if len(coords) > 3 and coords[0] == coords[-1]:
                    outers.append(Polygon(coords))
        if len(outers) == 1:
            rec = {"geometry": outers[0]}
            rec.update(rel.tags)
            records.append(rec)
            geom_types["Polygon"] = geom_types.get("Polygon", 0) + 1
        elif len(outers) > 1:
            rec = {"geometry": MultiPolygon(outers)}
            rec.update(rel.tags)
            records.append(rec)
            geom_types["MultiPolygon"] = geom_types.get("MultiPolygon", 0) + 1
    print(f"[Overpy] {layer['name']}: geom types {geom_types}")
    logging.info(f"[Overpy] {layer['name']}: geom types {geom_types}")
    gdf = gpd.GeoDataFrame(records, crs="EPSG:4326")
    return gdf
