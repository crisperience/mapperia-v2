from math import atan, cos, pi, sin, sqrt

from shapely import transform


def artisticTransformation(layer, aTs):
    """Apply artistic transformations to a GeoDataFrame layer.

    This function applies a series of geometric transformations to the geometries in a GeoDataFrame.
    Supported transformations include swirl, shift, and lens effects.

    Mathematical Background:
    - Swirl: Rotates points around a center with decreasing angle based on distance
    - Shift: Moves points along a line with decreasing effect based on distance
    - Lens: Creates a lens-like distortion from a center point

    Args:
        layer (GeoDataFrame): Input layer containing geometries to transform.
        aTs (list): List of transformation dictionaries. Each dict should contain:
            - type: str - One of "swirl", "shift", or "lens"
            - on: int - 1 to enable, 0 to disable the transformation
            For swirl transformations:
                - swirlCenter: tuple - Center point of the swirl (x,y)
                - swirlRadius: float - Radius of the swirl circle
                - maxRadiusSwirled: float - Maximum distance from swirl circle where points are affected
                - maxSwirlAngle: float - Maximum rotation angle in radians
            For shift transformations:
                - shiftPoint1: tuple - First point of the shift line (x,y)
                - shiftPoint2: tuple - Second point of the shift line (x,y)
                - shiftStrength: float - Distance to shift points along the line
                - maxDistanceShifted: float - Maximum distance from shift line where points are affected
            For lens transformations:
                - lensCentre: tuple - Center point of the lens effect (x,y)
                - lensRadius: float - Radius of the lens effect
                - lensStrength: float - Strength of the lens effect (0-100)

    Returns:
        GeoDataFrame: Layer with transformed geometries.

    Example:
        transformations = [
            {
                "type": "swirl",
                "on": 1,
                "swirlCenter": (0, 0),
                "swirlRadius": 1110,
                "maxRadiusSwirled": 1100,
                "maxSwirlAngle": 0.6
            }
        ]
        transformed_layer = artisticTransformation(layer, transformations)

    Notes:
        - The function densifies geometries before transformation for smoother results
        - Transformations are applied in sequence
        - Each transformation can be enabled/disabled independently
    """

    def gdf_transform(aT, x):
        """Apply a single transformation to a set of coordinates.

        This is an internal helper function that applies one transformation to a set of coordinates.
        It handles the coordinate transformation logic for all supported transformation types.

        Args:
            aT (dict): Transformation parameters dictionary.
            x (list): List of [x,y] coordinate pairs to transform.

        Returns:
            list: Transformed coordinate pairs.
        """

        def swirlTransformation(narray, aT):
            """Apply a swirl effect to coordinates.

            The swirl effect rotates points around a center point, with the rotation angle
            decreasing as points get further from the swirl circle.

            Mathematical formula:
            - Distance from swirl ring = |sqrt((x-x0)² + (y-y0)²) - r|
            - Angle = maxAngle * (sin((1 - d/maxRadius) * π - π/2) + 1) / 2
            - New coordinates use rotation matrix with calculated angle

            Args:
                narray (list): List of [x,y] coordinate pairs.
                aT (dict): Swirl parameters including center, radius, and angle.

            Returns:
                list: Swirled coordinate pairs.
            """
            i = 0
            for n in narray:
                aT["distanceFromSwirlRing"] = abs(
                    sqrt(
                        pow(n[0] - aT["swirlCenter"][0], 2)
                        + pow(n[1] - aT["swirlCenter"][1], 2)
                    )
                    - aT["swirlRadius"]
                )
                if aT["maxRadiusSwirled"] < aT["distanceFromSwirlRing"]:
                    aT["swirlAngle"] = 0
                else:
                    aT["swirlAngle"] = (
                        aT["maxSwirlAngle"]
                        * (
                            sin(
                                (
                                    (
                                        (
                                            aT["maxRadiusSwirled"]
                                            - aT["distanceFromSwirlRing"]
                                        )
                                        / aT["maxRadiusSwirled"]
                                    )
                                    * pi
                                    - pi / 2
                                )
                            )
                            + 1
                        )
                        / 2
                    )
                narray[i] = [
                    aT["swirlCenter"][0]
                    + (n[0] - aT["swirlCenter"][0]) * cos(aT["swirlAngle"])
                    - (n[1] - aT["swirlCenter"][1]) * sin(aT["swirlAngle"]),
                    aT["swirlCenter"][1]
                    + (n[0] - aT["swirlCenter"][0]) * sin(aT["swirlAngle"])
                    + (n[1] - aT["swirlCenter"][1]) * cos(aT["swirlAngle"]),
                ]
                i += 1
            return narray

        def shiftTransformation(narray, aT):
            """Apply a shift effect along a line.

            Points are shifted along a line defined by two points, with the shift amount
            decreasing as points get further from the line.

            Mathematical formula:
            - Line equation: ax + by + c = 0
            - Distance from line = |ax0 + by0 + c| / sqrt(a² + b²)
            - Shift amount = strength * (sin((1 - d/maxDist) * π - π/2) + 1) / 2

            Args:
                narray (list): List of [x,y] coordinate pairs.
                aT (dict): Shift parameters including line points and strength.

            Returns:
                list: Shifted coordinate pairs.
            """
            lineparam_a = (aT["shiftPoint2"][1] - aT["shiftPoint1"][1]) / (
                aT["shiftPoint2"][0] - aT["shiftPoint1"][1]
            )
            lineparam_b = -1
            lineparam_c = (
                (aT["shiftPoint2"][0] - aT["shiftPoint1"][0]) * aT["shiftPoint1"][1]
            ) - ((aT["shiftPoint2"][1] - aT["shiftPoint1"][1]) * aT["shiftPoint1"][0])
            steigung = (aT["shiftPoint2"][1] - aT["shiftPoint1"][1]) / (
                aT["shiftPoint2"][0] - aT["shiftPoint1"][0]
            )
            lineAngle = atan(steigung)
            i = 0
            for n in narray:
                aT["distance"] = abs(
                    lineparam_a * n[0] + lineparam_b * n[1] + lineparam_c
                ) / sqrt(lineparam_a * lineparam_a + lineparam_b * lineparam_b)
                if aT["distance"] > aT["maxDistanceShifted"]:
                    aT["shift"] = 0
                else:
                    aT["shift"] = (
                        aT["shiftStrength"]
                        * (
                            sin(
                                (
                                    (
                                        (aT["maxDistanceShifted"] - aT["distance"])
                                        / aT["maxDistanceShifted"]
                                    )
                                    * pi
                                    - pi / 2
                                )
                            )
                            + 1
                        )
                        / 2
                    )
                xinc = cos(lineAngle) * aT["shift"]
                yinc = sin(lineAngle) * aT["shift"]
                narray[i] = [n[0] + xinc, n[1] + yinc]
                i += 1
            return narray

        def lensTransformation(narray, aT):
            """Apply a lens distortion effect.

            Creates a lens-like distortion where points are pulled toward or pushed away
            from a center point, with the effect decreasing toward the edge of the lens.

            Mathematical formula:
            - Distance from center = sqrt((x-x0)² + (y-y0)²)
            - Shift amount = (strength/100) * radius * (sin((1 - d/radius) * π - π/2) + 1) / 2
            - Direction vector = (x-x0, y-y0) / distance

            Args:
                narray (list): List of [x,y] coordinate pairs.
                aT (dict): Lens parameters including center, radius, and strength.

            Returns:
                list: Distorted coordinate pairs.
            """
            i = 0
            for n in narray:
                aT["distanceFromLensCentre"] = sqrt(
                    pow(n[0] - aT["lensCentre"][0], 2)
                    + pow(n[1] - aT["lensCentre"][1], 2)
                )
                if aT["distanceFromLensCentre"] > aT["lensRadius"]:
                    aT["shift"] = 0
                else:
                    aT["shift"] = (
                        aT["lensStrength"]
                        / 100
                        * aT["lensRadius"]
                        * (
                            sin(
                                (
                                    (
                                        (
                                            aT["lensRadius"]
                                            - aT["distanceFromLensCentre"]
                                        )
                                        / aT["lensRadius"]
                                    )
                                    * pi
                                    - pi / 2
                                )
                            )
                            + 1
                        )
                        / 2
                    )
                a = n[1] - aT["lensCentre"][1]
                b = n[0] - aT["lensCentre"][0]
                c = sqrt(
                    pow(n[1] - aT["lensCentre"][1], 2)
                    + pow(n[0] - aT["lensCentre"][0], 2)
                )
                d = aT["shift"]
                xinc = b * d / c
                yinc = a * d / c
                narray[i] = [n[0] + xinc, n[1] + yinc]
                i += 1
            return narray

        if aT["on"] == 0:
            return x
        match aT["type"]:
            case "swirl":
                x = swirlTransformation(x, aT)
            case "shift":
                x = shiftTransformation(x, aT)
            case "lens":
                x = lensTransformation(x, aT)
            case _:
                print("transformation not known")
        return x

    for aT in aTs:
        layer["geometry"] = layer["geometry"].segmentize(1)
        layer["geometry"] = transform(layer["geometry"], lambda x: gdf_transform(aT, x))
    return layer
