from math import atan, cos, pi, sin, sqrt

from shapely import transform


def artisticTransformation(layer, aTs):
    def gdf_transform(aT, x):
        def swirlTransformation(narray, aT):
            # print(len(narray))
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

                """
                if aT["shift"]!=0:
                    aT["newDistanceFromLensCentre"] = sqrt(pow(narray[i][0] - aT["lensCentre"][0], 2) + pow(narray[i][1] - aT["lensCentre"][1], 2))
                    print(i, "lenscentre: ", aT["lensCentre"], " lensRadius: ", aT["lensRadius"], " shift: ", aT["shift"], "- olddistance: ", aT["distanceFromLensCentre"]," new distance",aT["newDistanceFromLensCentre"])
                """
                i += 1
                ##"type": "lens", "lensCentre": [100,100], "lensRadius": 50, "lensStrength": 5
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
        # layer=layer.explode(ignore_index=True)
        layer["geometry"] = layer["geometry"].segmentize(1)
        layer["geometry"] = transform(layer["geometry"], lambda x: gdf_transform(aT, x))
        # layer["geometry"]=layer["geometry"].simplify(.5)

    return layer


"""
I am expecting 
- the row of the geopandas dataframe and
- a dictionary with list of transformations like so
[
    {"type":"swirl","row"=, "swirlCenter":(0,0), "swirlRadius": 1110, "maxRadiusSwirled": 1100, "maxSwirlAngle":0.6},
    {"type":"swirl","row"=, "swirlCenter":(5,5), "swirlRadius": 1300, "maxRadiusSwirled": 1000, "maxSwirlAngle":0.3},
    {"type":"shift, "shiftPoint1":(2,2), shiftPoint2:(4.4), shiftDirection:"12", "shiftStrength":20, "maxDistanceShifted":30}
]

with
Type: Swirl
swirlCenter =  center of the swirlcircle. The Swirlcircle is the circle on which the points on the polygons are rotated with the maxSwirlAngle
swirlRadius =  radius of the swirlcircle
maxSwirlAngle = max angle of rotation on the swirlcircle
maxRadiusSwirled = maximal distance from swirl circle where points on the Polygons are still swirled

Type:Shift
shiftPoint1 = The first point of a shiftline along which the shift takes place
shiftPoint2 = The 2nd point of the shiftline along which the shift takes place
shiftDirection = a value of either 12 or two one meaning from Point 1 to Point 2 or the other way around
shiftStrength = the distance the points on the line a shifter in the shiftDirection
maxDistanceShifted = Points on the shiftLine will be shidted with by shiftStrength, all other points with in the maxShiftDistance from the ShiftLine will be shifted less. Points further away will not be shifted


Questions:

I am passing the whole row of the densified dataframe, transforming the geometries and returning the new gemetry, not the whole row...does that make sense?
I am using the densify polygon function out of the main script on the whole dataframe before any transformation. so the densify fucntion should be moved out of the transform_polygon script, right


Y=(y2-y1)/(x2-x1)*(x-x2)+y2

Y=(y2-y1)/(x2-x1)*X  +   Y1-x1(y2-y1)/(x2-x1)

Distance
a line equation  ax+by+c=0 for a line through two points shiftpoint1 ans shiftpoint 2 is
(shiftPoint2[1]-shiftPoint1[1])/(shiftPoint2[0]-shiftPoint1[1])*x + (-1)*y +  shiftPoint1[1]-shiftPoint1[0](shiftPoint2[1]-shiftPoint1[2])/(shiftPoint2[0]-shiftPoint1[0])

lineparam_a=(shiftPoint2[1]-shiftPoint1[1])/(shiftPoint2[0]-shiftPoint1[1])
lineparam_b=-1
lineparam_c=shiftPoint1[1]-shiftPoint1[0](shiftPoint2[1]-shiftPoint1[2])/(shiftPoint2[0]-shiftPoint1[0])

and
for a line ax+by+c=0 the distance of a point x0,y0=
distance=|Ax0+by0+c|/sqrt(a*a+b*b)

y=(-ax-c)/b= 
y=-a/bx - c/b




def transformPolygons(row,aTs):
    def transformCoordinates(aT,n):  # transform an individual coordinate
        def swirlTupleTransformation(n, artisticTransformation):
            aT["distanceFromSwirlRing"] = abs(sqrt(pow(n[0] - aT["swirlCenter"][0], 2) + pow(n[1] - aT["swirlCenter"][1], 2))- aT["swirlRadius"])
            if aT["maxRadiusSwirled"] < aT["distanceFromSwirlRing"]:
                aT["swirlAngle"] = 0
            else:
                aT["swirlAngle"] = (
                    aT["maxSwirlAngle"]
                    * (
                        sin(
                            (
                                (
                                    (aT["maxRadiusSwirled"] - aT["distanceFromSwirlRing"])
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
            newTuple = (
                aT["swirlCenter"][0]
                + (n[0] - aT["swirlCenter"][0]) * cos(aT["swirlAngle"])
                - (n[1] - aT["swirlCenter"][1]) * sin(aT["swirlAngle"]),
                aT["swirlCenter"][1]
                + (n[0] - aT["swirlCenter"][0]) * sin(aT["swirlAngle"])
                + (n[1] - aT["swirlCenter"][1]) * cos(aT["swirlAngle"]),
            )
            return newTuple
        def shiftTupleTransformation(n, aT):
            lineparam_a = (aT.shiftPoint2[1] - aT.shiftPoint1[1]) / (aT.shiftPoint2[0] - aT.shiftPoint1[1])
            lineparam_b = -1
            lineparam_c = aT.shiftPoint1[1] - aT.shiftPoint1[0](aT.shiftPoint2[1] - aT.shiftPoint1[2]) / (aT.shiftPoint2[0] - aT.shiftPoint1[0])
            aT["distance"] = abs(lineparam_a*n[0] +lineparam_b*n[1] + lineparam_c) / sqrt(lineparam_a * lineparam_a + lineparam_b* lineparam_b)
            if aT["distance"] > aT["maxDistanceShifted"]:
                aT["shift"]=0
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
            #### still needed: shift the point along the line and create a new tuple
            newTuple=n
            return newTuple
        def xyzTupleTransformation(n):
            ##new transformations are defined here for each tuple
            return newTuple

        match aT["type"]:
            case "Swirl":
                newTuple = swirlTupleTransformation(n, artisticTransformation)
            case "shift":
                newTuple = shiftTupleTransformation(n, artisticTransformation)
            case _:
                print("transformation not known")
                newTuple = n
        return newTuple

## applying the transformations one after the other
    for aT in aTs:
        exteriorP = Polygon(list(map(transformCoordinates, row["geometry"].exterior.coords, aT=aT)))  # go through exterior coodinates
        interior_rings = []
        for interior in row["geometry"].interiors:
            interior_rings.append(LinearRing(list(map(transformCoordinates, interior.coords[:]), aT=aT)))
        newgeom = Polygon(exteriorP.exterior, interior_rings)
    return newgeom




## old function renamed (it  worked)###
def transformPolygonsold(row, swirlCenter, swirlRadius, maxRadiusSwirled, maxSwirlAngle):
    def transformCoordinates(n):  # transform an individual coordinate
        distanceFromSwirlRing = abs(
            sqrt(pow(n[0] - swirlCenter[0], 2) + pow(n[1] - swirlCenter[1], 2))
            - swirlRadius
        )
        if maxRadiusSwirled < distanceFromSwirlRing:
            swirlAngle = 0
        else:
            swirlAngle = (
                maxSwirlAngle
                * (
                    sin(
                        (
                            (
                                (maxRadiusSwirled - distanceFromSwirlRing)
                                / maxRadiusSwirled
                            )
                            * pi
                            - pi / 2
                        )
                    )
                    + 1
                )
                / 2
            )
        newTuple = (
            swirlCenter[0]
            + (n[0] - swirlCenter[0]) * cos(swirlAngle)
            - (n[1] - swirlCenter[1]) * sin(swirlAngle),
            swirlCenter[1]
            + (n[0] - swirlCenter[0]) * sin(swirlAngle)
            + (n[1] - swirlCenter[1]) * cos(swirlAngle),
        )
        return newTuple

    exteriorP = Polygon(
        list(map(transformCoordinates, row["geometry"].exterior.coords))
    )  # go through exterior coodinates
    interior_rings = []
    for interior in row["geometry"].interiors:
        interior_rings.append(
            LinearRing(list(map(transformCoordinates, interior.coords[:])))
        )
        t1 = 1
    newgeom = Polygon(exteriorP.exterior, interior_rings)
    return newgeom
"""
