import logging
import os
import traceback
import warnings
from math import cos, sin

import cairosvg
import geopandas as gpd
import overpy
import pandas as pd
import pyproj
import svgwrite
from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from pyproj import Transformer
from shapely.geometry import (
    LineString,
    MultiLineString,
    MultiPolygon,
    Point,
    Polygon,
    box,
)
from shapely.ops import nearest_points, transform

warnings.filterwarnings("ignore")

# Ensure output directory exists
os.makedirs("output", exist_ok=True)


"""def densifyPolygons(row, maxSegmentLength):  # segmentize the geometry in ech row
    return row["geometry"].segmentize(maxSegmentLength)
"""


def initialize_crs(card):
    """Initialize the coordinate reference system for the card.

    Args:
        card (dict): Card configuration.

    Returns:
        tuple: Updated card dict and transformation function.
    """
    card["crs"] = createCRS(card["center"])
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
    """Create a custom Transverse Mercator CRS centered at the given coordinates."""
    # Create a custom CRS based on UTM for the card center.
    crs_string = (
        f"+lon_0={cardcenter[1]} +ellps=WGS84 +y_0=0 +no_defs=True +proj=tmerc +x_0=0 "
        f"+units=m  +lat_0={cardcenter[0]} +k=1 +towgs84=0,0,0,0,0,0,0"
    )
    crs = pyproj.CRS.from_string(crs_string)
    return crs


def calculate_card_dimensions(card, ax, fig, cm, cwidthmm):
    """Calculate card scale and dimensions for plotting.

    Args:
        card (dict): Card configuration.
        ax: Matplotlib axis.
        fig: Matplotlib figure.
        cm (float): Centimeter to inch conversion.
        cwidthmm (float): Connector width in mm.

    Returns:
        tuple: Updated card dict and connector width in map units.
    """
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
    cwidth = (
        cwidthmm * card["scale"] / 1000
    )  # Connector width in map units (mm on paper).
    card["layers"][2]["graphbuffer"] = max(
        card["layers"][2]["graphbuffer_mm"] * card["scale"] / 1000, 2
    )  # Streets: min graphbuffer_mm or, if wider, 1m buffered on map.
    card["layers"][4]["graphbuffer"] = max(
        card["layers"][4]["graphbuffer_mm"] * card["scale"] / 1000, 2
    )  # Rail: min graphbuffer_mm or, if wider, 1m buffered on map.
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
    """Convert decimal latitude and longitude to DMS string.
    Note: This function may not handle negative latitudes correctly for New Zealand coordinates. Needs review for southern hemisphere edge cases.
    """
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
    """Deprecated. Not used."""
    # Function removed as it is not used anywhere in the backend.
    pass


def c_lines(gdf, layer_name, clinesyes, card, cwidth):
    """Generate connector lines between polygons in a layer.

    Args:
        gdf (GeoDataFrame): Input geometries.
        layer_name (str): Name of the layer.
        clinesyes (bool): Whether to generate connector lines.
        card (dict): Card configuration.
        cwidth (float): Connector width.

    Returns:
        tuple: Updated GeoDataFrame and connector lines GeoDataFrame.
    """
    gdf = gdf.dissolve().explode(index_parts=True)
    all_clines = gpd.GeoDataFrame(columns=["geometry"], crs=card["crs"])
    if (
        cwidth == 0 or clinesyes == 0
    ):  # Return an empty connector lines GeoDataFrame if not needed.
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
    ):  # Continue until the ring area covers most of the frame area.
        ring = outer_hull.difference(
            outer_hull.buffer(total_buffer)
        )  # Shrink the ring area from outside to inside.
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
        ):  # Add connector lines and dissolve until only one polygon remains (all connected).
            for index, row in subgdf.iterrows():
                point = row.geometry  # Geometry to check against the rest.
                multipoint = subgdf.drop(
                    index, axis=0
                ).geometry.unary_union  # Union of all other geometries.
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
                    ).centroid  # Move the point inwards to avoid sharp corners in the connector line.
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
    """Generate connector lines for polygons within grid cells."""
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
    # TODO: Make all variables accessible via card["key"] if needed for future features.
    rectangles = []  # List of rectangles from outer to inner frame.
    grids = []  # List of grid rectangles for connector lines.
    # Create the outer frame rectangle.
    # Note: With the card CRS centered around the card center, the center here is 0, but not for other CRSs. Consider using the projected card center if needed.
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
    # Iterate from the outer frame rectangle to the center with frame width steps.
    for i in range(
        1, int(((min(card_height, card_width) + frame_width * 2) / 2) // frame_width)
    ):
        rectangles.append(rectangles[0].buffer(-frame_width * i, join_style=2))

    # Secondary grid of rectangles for additional connector lines.
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
                    )  # TODO: Refactor to use rectangles instead of squares for grid cells.

    match card["binding"]:
        case (
            "libro"
        ):  # Create a polygon with a double-thick left frame for the binding.
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
        case "frame":  # Create a polygon for the binding frame.
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
        case _:  # Default: handle 'globe' or other types.
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


def simplify(layer, gdf, card_area):
    """Simplify geometries in a layer based on configuration.

    Args:
        layer (dict): Layer configuration
        gdf (GeoDataFrame): Input geometries
        card_area (float): Total card area for size filtering

    Returns:
        GeoDataFrame: Simplified geometries
    """
    if layer["type"] == "graph":
        gdf["geometry"] = gdf["geometry"].buffer(
            layer["graphbuffer"]
        )  # Convert line streets into polygons
    else:
        gdf["geometry"] = gdf["geometry"].buffer(layer["emph_buffer"], join_style=1)
        gdf["geometry"] = gdf["geometry"].buffer(layer["smooth"], join_style=1)

        gdf = gdf.dissolve()
        gdf["geometry"] = gdf["geometry"].unary_union
        gdf = gdf.explode(
            index_parts=True
        )  # Explode MultiPolygon geometries into individual Polygon geometries

        gdf["geometry"] = gdf["geometry"].simplify(
            layer["simplify_tolerance"]
        )  # Simplify geometry edges
        gdf["geometry"] = gdf["geometry"].buffer(
            -layer["smooth"], join_style=1
        )  # Reverse the applied buffer

        # Reduce buffer and expand again to remove small thin pieces
        gdf["geometry"] = gdf["geometry"].buffer(
            -layer["smooth"], join_style=1
        )  # Reverse the applied buffer
        gdf["geometry"] = gdf["geometry"].buffer(layer["smooth"], join_style=1)

        gdf["area"] = gdf["geometry"].area  # Add a column with area
        gdf = gdf[
            gdf["area"] > (card_area * layer["pat"])
        ]  # Filter out polygons that do not meet the threshold criteria
        gdf = gdf.drop("area", axis=1)
    return gdf


def write_svg_laser_output(gdf, output_path, width=1000, height=1000):
    """Write geometries to SVG file in laser-ready format.

    Args:
        gdf: GeoDataFrame containing geometries
        output_path: Path to save SVG file
        width: SVG width in points (default 1000)
        height: SVG height in points (default 1000)
    """
    # Get bounds of all geometries
    bounds = gdf.total_bounds

    # Create SVG document with full profile
    dwg = svgwrite.Drawing(output_path, size=(width, height), profile="full")

    # Add XML declaration and DOCTYPE
    dwg.attribs["xmlns"] = "http://www.w3.org/2000/svg"
    dwg.attribs["xmlns:xlink"] = "http://www.w3.org/1999/xlink"
    dwg.attribs["version"] = "1.1"
    dwg.attribs["viewBox"] = f"0 0 {width} {height}"

    # Draw each geometry as a path with round linecap and linejoin
    for geom in gdf.geometry:
        path_data = geometry_to_svg_path(geom, bounds, width, height)
        if path_data:
            dwg.add(
                dwg.path(
                    d=path_data,
                    fill="none",
                    stroke=layer["color"],
                    stroke_width=layer.get("stroke_width", 1),
                    stroke_linecap="round",
                    stroke_linejoin="round",
                )
            )

    # Write SVG file
    with open(output_path, "w") as f:
        f.write('<?xml version="1.0" encoding="utf-8" standalone="no"?>\n')
        f.write(
            '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n'
        )
        f.write(dwg.tostring())


def geometry_to_svg_path(geometry, width=1000, height=1000, bounds=None):
    """Convert a Shapely geometry to SVG path data.

    Args:
        geometry: Shapely geometry object
        width: SVG width in points (default 1000)
        height: SVG height in points (default 1000)
        bounds: Bounding box of all geometries (minx, miny, maxx, maxy)

    Returns:
        str: SVG path data string
    """
    if bounds is None:
        bounds = geometry.bounds

    minx, miny, maxx, maxy = bounds
    scale_x = width / (maxx - minx)
    scale_y = height / (maxy - miny)

    def transform_coords(x, y):
        # Transform coordinates to SVG space
        svg_x = (x - minx) * scale_x
        svg_y = height - (y - miny) * scale_y  # Flip Y axis
        return f"{svg_x:.2f} {svg_y:.2f}"

    if isinstance(geometry, Polygon):
        paths = []
        # Only one exterior ring per geometry
        coords = list(geometry.exterior.coords)
        path_data = f"M {transform_coords(coords[0][0], coords[0][1])}"
        for x, y in coords[1:]:
            path_data += f"\nL {transform_coords(x, y)}"
        path_data += "\nz"
        paths.append(path_data)
        # Each interior ring (hole) as its own path segment
        for interior in geometry.interiors:
            coords = list(interior.coords)
            path_data = f"M {transform_coords(coords[0][0], coords[0][1])}"
            for x, y in coords[1:]:
                path_data += f"\nL {transform_coords(x, y)}"
            path_data += "\nz"
            paths.append(path_data)
        return paths
    elif isinstance(geometry, MultiPolygon):
        # Only one path per exterior ring of each polygon
        paths = []
        for poly in geometry.geoms:
            coords = list(poly.exterior.coords)
            path_data = f"M {transform_coords(coords[0][0], coords[0][1])}"
            for x, y in coords[1:]:
                path_data += f"\nL {transform_coords(x, y)}"
            path_data += "\nz"
            paths.append(path_data)
            # Each interior ring (hole) as its own path segment
            for interior in poly.interiors:
                coords = list(interior.coords)
                path_data = f"M {transform_coords(coords[0][0], coords[0][1])}"
                for x, y in coords[1:]:
                    path_data += f"\nL {transform_coords(x, y)}"
                path_data += "\nz"
                paths.append(path_data)
        return paths
    elif isinstance(geometry, LineString):
        coords = list(geometry.coords)
        path_data = f"M {transform_coords(coords[0][0], coords[0][1])}"
        for x, y in coords[1:]:
            path_data += f"\nL {transform_coords(x, y)}"
        return [path_data]
    elif isinstance(geometry, MultiLineString):
        paths = []
        for line in geometry.geoms:
            coords = list(line.coords)
            path_data = f"M {transform_coords(coords[0][0], coords[0][1])}"
            for x, y in coords[1:]:
                path_data += f"\nL {transform_coords(x, y)}"
            paths.append(path_data)
        return paths
    return []


def generate_map(lat, lon, width, height, formats, style, location, layers):
    """Generate a map at the specified location and size.

    Args:
        lat (float): Latitude of map center
        lon (float): Longitude of map center
        width (float): Width of map in km (will be converted to meters)
        height (float): Height of map in km (will be converted to meters)
        formats (list): List of output formats ('png', 'svg_laser')
        style (str): Map style
        location (str): Location name
        layers (list): List of layer names (strings)
    """
    try:
        CANONICAL_LAYERS = [
            {
                "name": "buildings",
                "active": True,
                "type": "polygon",
                "lcolor": ["02red", "#FF0000"],
                "osm_tags": {"building": True},
            },
            {
                "name": "greens",
                "active": True,
                "type": "polygon",
                "lcolor": ["03green", "#00FF00"],
                "osm_tags": {
                    "leisure": ["park", "garden", "golf_course", "recreation_ground"],
                    "landuse": ["forest", "grass"],
                },
            },
            {
                "name": "streets",
                "active": True,
                "type": "graph",
                "lcolor": ["08lightgrey", "#B4B4B4"],
                "osm_tags": {"highway": True},
            },
            {
                "name": "blues",
                "active": True,
                "type": "polygon",
                "lcolor": ["01blue", "#0000FF"],
                "osm_tags": {
                    "waterway": True,
                    "natural": "water",
                    "landuse": "reservoir",
                    "water": True,
                },
            },
            {
                "name": "rail",
                "active": True,
                "type": "graph",
                "lcolor": ["09darkgrey", "#808080"],
                "osm_tags": {"railway": True},
            },
            {
                "name": "front_cover",
                "active": True,
                "type": "polygon",
                "lcolor": ["20darkred", "#D33F6A"],
                "osm_tags": {},
            },
            {
                "name": "back_cover",
                "active": True,
                "type": "polygon",
                "lcolor": ["20darkred", "#D33F6A"],
                "osm_tags": {},
            },
        ]
        LAYER_ORDER = [
            "back_cover",
            "blues",
            "greens",
            "rail",
            "streets",
            "buildings",
            "front_cover",
        ]
        # Convert width and height from kilometers to meters
        width = width * 1000
        height = height * 1000
        # Build the canonical layer list, only including requested layers
        requested_layer_names = set(layers)
        layers = [
            layer.copy()
            for layer in CANONICAL_LAYERS
            if layer["name"] in requested_layer_names
        ]
        card = {
            "center": [lat, lon],
            "width": width,
            "height": height,
            "layers": layers,
            "style": style,
            "location": location,
            "formats": formats,
            "frame_width_mm": 5,
            "hole_radius_for_libro_in_mm": 1.5,
            "hole_radius_for_globe_in_mm": 1.5,
            "hole_radius_for_frame_in_mm": 1.5,
            "number_cline_buffer": 10,
        }
        card["layers"] = sorted(
            card["layers"],
            key=lambda layer: LAYER_ORDER.index(layer["name"])
            if layer["name"] in LAYER_ORDER
            else 999,
        )
        # Initialize CRS and transformations
        card, t_crs2WGS84 = initialize_crs(card)
        # Create frame
        card["frame"] = gpd.GeoDataFrame(
            geometry=[box(-width / 2, -height / 2, width / 2, height / 2)],
            crs=card["crs"],
        )
        # Process each layer
        all_layers_gdf = {}
        for layer in card["layers"]:
            if layer["active"]:
                bbox = card["frame"].to_crs("epsg:4326").total_bounds
                gdf = fetch_osm_layer_overpy(layer, bbox)
                if not gdf.empty:
                    gdf = gdf.to_crs(card["crs"])
                    gdf = gdf.clip(card["frame"])
                    all_layers_gdf[layer["name"]] = gdf
        output_urls = {}
        # Get card frame bounds for consistent SVG scaling
        card_frame_bounds = card["frame"].total_bounds
        # Generate SVG laser ready (outlines only)
        if "svg_laser" in formats:
            svg_laser_path = f"output/{location}_laser.svg"
            dwg = svgwrite.Drawing(svg_laser_path, size=(1000, 1000), profile="full")
            dwg.attribs["xmlns"] = "http://www.w3.org/2000/svg"
            dwg.attribs["xmlns:xlink"] = "http://www.w3.org/1999/xlink"
            dwg.attribs["version"] = "1.1"
            dwg.attribs["viewBox"] = "0 0 1000 1000"
            for layer in card["layers"]:
                layer_name = layer["name"]
                if layer_name in all_layers_gdf:
                    gdf = all_layers_gdf[layer_name]
                    layer_group = dwg.g(id=f"layer_{layer_name}")
                    stroke_width = 0.1  # Fixed thin stroke for laser cutting
                    for i, (_, row) in enumerate(gdf.iterrows()):
                        path_datas = geometry_to_svg_path(
                            row.geometry, 1000, 1000, card_frame_bounds
                        )
                        if path_datas:
                            if not isinstance(path_datas, list):
                                path_datas = [path_datas]
                            for j, path_data in enumerate(path_datas):
                                if path_data:
                                    path = dwg.path(
                                        d=path_data,
                                        stroke=layer["lcolor"][1],
                                        fill="none",
                                        stroke_width=stroke_width,
                                        id=f"{layer_name}_path{i}_{j}",
                                    )
                                    layer_group.add(path)
                    dwg.add(layer_group)
            with open(svg_laser_path, "w") as f:
                f.write('<?xml version="1.0" encoding="utf-8" standalone="no"?>\n')
                f.write(
                    '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n'
                )
                f.write(dwg.tostring())
            output_urls["svgLaserUrl"] = f"/download/{location}_laser.svg"
        # Generate PNG preview from SVG laser if requested
        if "png" in formats:
            # Calculate output size in mm (assuming width/height are in meters)
            width_mm = width
            height_mm = height
            # Convert mm to inches
            width_in = width_mm / 25.4
            height_in = height_mm / 25.4
            # 150 DPI for preview
            dpi = 150
            # Calculate pixel dimensions at 150 DPI
            px_width = width_in * dpi
            px_height = height_in * dpi
            # Adaptive max dimension based on map size
            map_km = max(width, height) / 1000
            if map_km <= 3:
                max_dim = 3000
            elif map_km <= 10:
                max_dim = 4000
            elif map_km <= 20:
                max_dim = 5000
            else:
                max_dim = 6000
            if px_width >= px_height:
                scale = min(1.0, max_dim / px_width)
            else:
                scale = min(1.0, max_dim / px_height)
            png_width = int(px_width * scale)
            png_height = int(px_height * scale)
            print(
                f"[PNG Export] width_mm={width_mm}, height_mm={height_mm}, px=({png_width},{png_height}), dpi={dpi}, max_dim={max_dim}"
            )
            png_path = f"output/{location}_preview.png"
            cairosvg.svg2png(
                url=svg_laser_path,
                write_to=png_path,
                output_width=png_width,
                output_height=png_height,
                dpi=dpi,
            )
            output_urls["pngUrl"] = f"/download/{location}_preview.png"
        # Debug: Print blues layer bounds and geometry count
        if "blues" in all_layers_gdf:
            blues_gdf = all_layers_gdf["blues"]
            print(f"[DEBUG] Blues layer bounds: {blues_gdf.total_bounds}")
            print(f"[DEBUG] Blues layer geometry count: {len(blues_gdf)}")
        return output_urls

    except Exception as e:
        logging.error(f"Error generating map: {str(e)}")
        logging.error(traceback.format_exc())
        raise


def fetch_osm_layer_overpy(layer, bbox):
    """
    Fetch OSM data for the given layer and bbox using overpy and return a GeoDataFrame.
    bbox: (minx, miny, maxx, maxy)
    layer: dict with osm_tags and name
    """
    api = overpy.Overpass()
    # Correct bbox order: minx, miny, maxx, maxy = west, south, east, north
    west, south, east, north = bbox

    if layer["name"] == "buildings":
        ql = f"""
        (
          way["building"]({south},{west},{north},{east});
          relation["building"]({south},{west},{north},{east});
        );
        (._;>;);
        out body;
        """
    elif layer["name"] == "greens":
        ql = f"""
        (
          way["leisure"~"park|garden|golf_course|recreation_ground"]({south},{west},{north},{east});
          relation["leisure"~"park|garden|golf_course|recreation_ground"]({south},{west},{north},{east});
          way["landuse"~"forest|grass"]({south},{west},{north},{east});
          relation["landuse"~"forest|grass"]({south},{west},{north},{east});
        );
        (._;>;);
        out body;
        """
    elif layer["name"] == "streets":
        ql = f"""
        (
          way["highway"]({south},{west},{north},{east});
          relation["highway"]({south},{west},{north},{east});
        );
        (._;>;);
        out body;
        """
    elif layer["name"] == "blues":
        ql = f"""
        (
          way["waterway"]({south},{west},{north},{east});
          relation["waterway"]({south},{west},{north},{east});
          way["natural"="water"]({south},{west},{north},{east});
          relation["natural"="water"]({south},{west},{north},{east});
          way["landuse"="reservoir"]({south},{west},{north},{east});
          relation["landuse"="reservoir"]({south},{west},{north},{east});
          way["water"]({south},{west},{north},{east});
          relation["water"]({south},{west},{north},{east});
        );
        (._;>;);
        out body;
        """
    elif layer["name"] == "rail":
        ql = f"""
        (
          way["railway"]({south},{west},{north},{east});
          relation["railway"]({south},{west},{north},{east});
        );
        (._;>;);
        out body;
        """
    else:
        ql_parts = []
        for key, value in layer["osm_tags"].items():
            if value is True:
                ql_parts.append(
                    f'way["{key}"]({south},{west},{north},{east}); relation["{key}"]({south},{west},{north},{east});'
                )
            elif isinstance(value, list):
                for v in value:
                    ql_parts.append(
                        f'way["{key}"="{v}"]({south},{west},{north},{east}); relation["{key}"="{v}"]({south},{west},{north},{east});'
                    )
            else:
                ql_parts.append(
                    f'way["{key}"="{value}"]({south},{west},{north},{east}); relation["{key}"="{value}"]({south},{west},{north},{east});'
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
    print(
        f"[Overpy] {layer['name']}: {len(result.ways)} ways, {len(result.relations)} relations"
    )
    geom_types = {}
    records = []
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
        rec = {"geometry": geom}
        rec.update(way.tags)
        records.append(rec)
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
    if not records:
        return gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")
    gdf = gpd.GeoDataFrame(records, crs="EPSG:4326")
    return gdf


app = FastAPI()


@app.get("/download/{filename}")
async def download_file(filename: str):
    file_path = f"output/{filename}"
    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="File not found")
    return FileResponse(
        file_path, media_type="application/octet-stream", filename=filename
    )


@app.post("/api/generate-map")
async def generate_map_endpoint(request: dict):
    try:
        result = generate_map(
            lat=request["lat"],
            lon=request["lon"],
            width=request["width"],
            height=request["height"],
            formats=request["formats"],
            style=request["style"],
            location=request["location"],
            layers=request["layers"],
        )
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
