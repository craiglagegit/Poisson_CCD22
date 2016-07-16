#!/usr/bin/env python

#Author: Craig Lage, NYU;
#Date: 16-Nov-15

#This program plots the pixel area plots from the Poisson CCD solver
import matplotlib
matplotlib.use("PDF")
from pylab import *
import sys, time, h5py

#****************SUBROUTINES*****************


class Point:
    def __init__(self,x,y):
        self.x = x
        self.y = y

    def Rotate(self,theta):
        xold = self.x
        yold = self.y
        self.x = xold * cos(theta) - yold * sin(theta)
        self.y = xold * sin(theta) + yold * cos(theta)

    def Shift(self,deltax,deltay):
        self.x = self.x + deltax
        self.y = self.y + deltay

class LineSegment:
    def __init__(self,point1,point2):
        self.p1 = point1
        self.p2 = point2
    def Rotate(self,theta):
        self.p1.Rotate(theta)
        self.p2.Rotate(theta)
    def Shift(self,deltax,deltay):
        self.p1.Shift(deltax,deltay)
        self.p2.Shift(deltax,deltay)

class Polygon:
    def __init__(self,pointlist):
        self.npoints = len(pointlist)
        self.pointlist = pointlist
        self.sorted = False

    def AddPoint(self,point):
        self.pointlist.append(point)
        self.npoints = self.npoints + 1
        self.sorted = False

    def Rotate(self,theta):
        for point in self.pointlist:
            point.Rotate(theta)

    def Shift(self,deltax,deltay):
        for point in self.pointlist:
            point.Shift(deltax,deltay)

    def Sort(self):
        # calculate centroid of the polygon
        if self.npoints < 3:
            self.sorted = True
            return
        else:
            cx = float(sum(point.x for point in self.pointlist)) / self.npoints
            cy = float(sum(point.y for point in self.pointlist)) / self.npoints
            # create a new list of poly which includes angles
            sortedpoints = []
            for point in self.pointlist:
                an = (arctan2(point.y - cy, point.x - cx) + 2.0 * pi) % (2.0 * pi)
                sortedpoints.append((point, an))
            # sort it using the angles

            sortedpoints.sort(key = lambda tup: tup[1])
            # replace the pointlist with the sorted list w/ angles removed
	    for i, sortedpoint in enumerate(sortedpoints):
	    	self.pointlist[i] = sortedpoint[0]
            self.sorted = True
            return

    def Area(self):
        # Calculates area of a polygon using the shoelace algorithm
        if not self.sorted:
            self.Sort() # Polygon points must be in CCW order
        self.area = 0.0
        for i in range(self.npoints):
            j = (i + 1) % self.npoints
            self.area += self.pointlist[i].x * self.pointlist[j].y
            self.area -= self.pointlist[j].x * self.pointlist[i].y
        self.area = abs(self.area) / 2.0
        return self.area

    def PointInside(self,point):
    # Determines if a given point is inside the polygon
        if not self.sorted:
            self.Sort() # Polygon points must be in CCW order
        inside =False
        x1 = self.pointlist[0].x
        y1 = self.pointlist[0].y
        for i in range(self.npoints+1):
            x2 = self.pointlist[i % self.npoints].x
            y2 = self.pointlist[i % self.npoints].y
            if point.y > min(y1,y2):
                if point.y <= max(y1,y2):
                    if point.x <= max(x1,x2):
                        if y1 != y2:
                            xinters = (point.y-y1)*(x2-x1)/(y2-y1)+x1
                        if x1 == x2 or point.x <= xinters:
                            inside = not inside
            x1,y1 = x2,y2
        return inside

def CCW(p1,p2,p3):
    # Determines whether three points are in CCW order
    return (p3.y-p1.y)*(p2.x-p1.x) > (p2.y-p1.y)*(p3.x-p1.x)

def DoesIntersect(line1,line2):
    # Determines whether two line segments intersect
    return CCW(line1.p1,line2.p1,line2.p2) != CCW(line1.p2,line2.p1,line2.p2) and \
        CCW(line1.p1,line1.p2,line2.p1) != CCW(line1.p1,line1.p2,line2.p2)

def LineSegmentIntersection(l1, l2):
    # Finds the intersection point of two line segments
    if not DoesIntersect(l1,l2):
        return None
    else:
        denominator = (l1.p1.x - l1.p2.x) * (l2.p1.y - l2.p2.y) - (l1.p1.y - l1.p2.y) * (l2.p1.x - l2.p2.x)
        x = ((l1.p1.x * l1.p2.y - l1.p1.y * l1.p2.x) * (l2.p1.x - l2.p2.x) - (l1.p1.x - l1.p2.x) *\
             (l2.p1.x * l2.p2.y - l2.p1.y * l2.p2.x) ) / (denominator)
        y = ((l1.p1.x * l1.p2.y - l1.p1.y * l1.p2.x) * (l2.p1.y - l2.p2.y) - (l1.p1.y - l1.p2.y) *\
             (l2.p1.x * l2.p2.y - l2.p1.y * l2.p2.x) ) / (denominator)
        return Point(x, y)

def PolygonIntersectionArea(poly1, poly2):
    # Calculates area of intersection of two polygons
    # Polygon of intersection is the union of:
    #   (A) Vertices of poly1 inside poly2
    #   (B) Vertices of poly2 inside poly1
    #   (C) Intersection points of the polygon segments
    poly_int = Polygon([]) # this is the polygon of intersection
    for point in poly1.pointlist:
        if poly2.PointInside(point):
            poly_int.AddPoint(point)
    for point in poly2.pointlist:
        if poly1.PointInside(point):
            poly_int.AddPoint(point)
    for i, point1 in enumerate(poly1.pointlist):
        line1 = LineSegment(point1, poly1.pointlist[(i+1) % poly1.npoints])
        for j, point2 in enumerate(poly2.pointlist):
            line2 = LineSegment(point2, poly2.pointlist[(j+1) % poly2.npoints])
            if DoesIntersect(line1,line2):
                int_point = LineSegmentIntersection(line1, line2)
                poly_int.AddPoint(int_point)
    poly_int.Sort()
    if poly_int.npoints < 3:
        return 0.0
    else:
        return poly_int.Area()


def ReadConfigFile(filename):
    # This reads the config file for the necessary settings
    # and returns a dictionary with the values
    file = open(filename,'r')
    lines=file.readlines()
    file.close()
    ConfigData = {}
    try:
        for line in lines:
            ThisLine=line.strip().split()
            ThisLineLength=len(ThisLine)
            if ThisLineLength < 3:
                continue
            if list(ThisLine[0])[0]=='#' or ThisLine[0]=='\n':
                continue
            try:
                ParamName = ThisLine[0]
                ThisLine.remove(ThisLine[0])
                for counter,item in enumerate(ThisLine):
                    if list(item)[0] == '#':
                        del ThisLine[counter:] # Strip the rest of the line as a comment
                        continue
                    if item == '=':
                        ThisLine.remove(item)
                        continue
                if len(ThisLine) == 0:
                    continue
                elif len(ThisLine) == 1:
                    ThisParam = ThisLine[0]
                    try: ConfigData[ParamName] = int(ThisParam)
                    except ValueError:
                        try: ConfigData[ParamName] = float(ThisParam)
                        except ValueError:
                            try:
                                ConfigData[ParamName] = ThisParam
                            except ValueError:
                                print "Error reading .cfg file"
                else:
                    ThisParam = []
                    for item in ThisLine:
                        try: ThisParam.append(int(item))
                        except ValueError:
                            try: ThisParam.append(float(item))
                            except ValueError:
                                ThisParam.append(item)
                    ConfigData[ParamName] = ThisParam
            except:
                continue
    except:
        print "Error reading .cfg file"

    return ConfigData

def ReadAreaFile(filename, nx, ny):
    area = zeros([nx, ny])
    file = open(filename, 'r')
    lines = file.readlines()
    file.close()
    for line in lines:
        items = line.split()
        if items[0] == "Nx":
            continue
        i = int(items[0])
        j = int(items[1])
        area[i,j] = float(items[2])
    return area

def ReadVertexFile(filename, configdata, poly):
    file = open(filename, 'r')
    lines = file.readlines()
    file.close()
    lasti = 3457
    lastj = 3457
    for line in lines:
        items = line.split()
        if items[0] == "X0":
            continue
        i = (int(float(items[0])) - int(1.5*configdata["PixelSize"])) / int(configdata["PixelSize"])
        j = (int(float(items[1])) - int(1.5*configdata["PixelSize"])) / int(configdata["PixelSize"])
        point = Point(float(items[3]), float(items[4]))
        poly[i][j].AddPoint(point)

    for i in range(Nx):
        for j in range(Ny):
            poly[i][j].Sort()
    return

def CreatePolyList(nx, ny):
    polylist = [] # This is the undistorted pixels
    for i in range(nx):
        xlist = []
        for j in range(ny):
            pointlist = []
            xlist.append(Polygon(pointlist))
        polylist.append(xlist)
    return polylist

def FillUndistortedPolyList(configdata, poly):
    nx = configdata["PixelBoundaryNx"]
    ny = configdata["PixelBoundaryNy"]

    pixelsize = configdata["PixelSize"]
    numvertices = configdata["NumVertices"]
    dtheta = pi / (2.0 * (numvertices + 1.0))
    theta0 = - pi / 4.0
    for i in range(nx):
        for j in range(ny):
            xcenter = configdata["PixelBoundaryLowerLeft"][0] + (i + 0.5) * pixelsize
            ycenter = configdata["PixelBoundaryLowerLeft"][1] + (j + 0.5) * pixelsize
            # First the corners
            for xsign in [-1.0,1.0]:
                for ysign in [-1.0,1.0]:
                    point = Point(xcenter + xsign * pixelsize / 2.0, ycenter + ysign * pixelsize / 2.0)
                    poly[i][j].AddPoint(point)
            # Then the edges
            for xsign in [-1.0,1.0]:
                for n in range(NumVertices):
                    theta = theta0 + (n + 1.0) * dtheta
                    point = Point(xcenter + xsign * pixelsize / 2.0, ycenter + ysign * tan(theta) * pixelsize / 2.0)
                    poly[i][j].AddPoint(point)
            for ysign in [-1.0,1.0]:
                for n in range(NumVertices):
                    theta = theta0 + (n + 1.0) * dtheta
                    point = Point(xcenter + xsign * tan(theta) * pixelsize / 2.0, ycenter + ysign * pixelsize / 2.0)
                    poly[i][j].AddPoint(point)
            poly[i][j].Sort()
    return

def FillPolyShiftsList(configdata, poly0, polyq, polyd):
    # Calculates the displacement of each vertex per electron.
    # These aren't really polygons (don't sort!), but is a convenient
    # way to store the data
    nx = configdata["PixelBoundaryNx"]
    ny = configdata["PixelBoundaryNy"]
    numelec = configdata["CollectedCharge_0_0"]
    for i in range(nx):
        for j in range(ny):
            for k in range(poly0[i][j].npoints):
                dx = (polyq[i][j].pointlist[k].x - poly0[i][j].pointlist[k].x) / float(numelec)
                dy = (polyq[i][j].pointlist[k].y - poly0[i][j].pointlist[k].y) / float(numelec)
                point = Point(dx, dy)
                polyd[i][j].AddPoint(point)
    return

def FillModelPolyList(configdata, charges, poly0, polyd, polym):
    # Calculates the modeled poly shapes, given the charges,
    # The undistorted polys, and the poly shifts
    nx = configdata["PixelBoundaryNx"]
    ny = configdata["PixelBoundaryNy"]
    nxcenter = (nx - 1) / 2
    nycenter = (ny - 1) / 2
    for i in range(nx):
        for j in range(ny):
            for k in range(poly0[i][j].npoints):
                x = poly0[i][j].pointlist[k].x
                y = poly0[i][j].pointlist[k].y
                point = Point(x, y)
                polym[i][j].AddPoint(point)

            for chargei in range(nx):
                for chargej in range(ny):
                    if charges[chargei, chargej] == 0:
                        continue
                    ii = i - (chargei - nxcenter)
                    jj = j - (chargej - nycenter)
                    if ii < 0 or ii > (nx - 1) or jj < 0 or jj > (ny - 1):
                        continue
                    for k in range(poly0[i][j].npoints):
                        dx = polyd[ii][jj].pointlist[k].x * charges[chargei, chargej]
                        dy = polyd[ii][jj].pointlist[k].y * charges[chargei, chargej]
                        polym[i][j].pointlist[k].x += dx
                        polym[i][j].pointlist[k].y += dy
    return

def PlotPixels(ax, title, poly):
    ax.set_title(title)
    for i in range(NxCenter-PlotDelta, NxCenter+PlotDelta+1):
        ax.plot([10.0 + 10.0 * i, 10.0 + 10.0 * i], [40.0, 70.0], color = 'black', ls = '--')
    for j in range(NyCenter-PlotDelta, NyCenter+PlotDelta+1):
        ax.plot([40.0, 70.0], [10.0 + 10.0 * j, 10.0 + 10.0 * j], color = 'black', ls = '--')

    for i in range(NxCenter-PlotDelta, NxCenter+PlotDelta+1):
        for j in range(NyCenter-PlotDelta, NyCenter+PlotDelta+1):
            xs = []
            ys = []
            for k in range(NumAngles):
                xs.append(poly[i][j].pointlist[k].x)
                ys.append(poly[i][j].pointlist[k].y)
            xs.append(poly[i][j].pointlist[0].x)
            ys.append(poly[i][j].pointlist[0].y)
            ax.plot(xs, ys, lw = 0.5)
            area = poly[i][j].Area()

            textcolor = 'black'
            if charges[i,j] > 0:
                numelec = charges[i,j]
                textcolor = 'red'
                ax.text(12.0 + 10.0 * i, 12.0 + 10.0 * j, "%d"%numelec, color = textcolor, fontsize = 6, fontweight = 'bold')

            ax.text(12.0 + 10.0 * i, 14.0 + 10.0 * j, "%.4f"%area, color = textcolor, fontsize = 6, fontweight = 'bold')
    return

def CalculatePixelErrors(configdata, poly1, poly2):
    # Calculates the errors between two pixel arrays
    nx = configdata["PixelBoundaryNx"]
    ny = configdata["PixelBoundaryNy"]
    AveAreaError = 0.0
    AveVertexError = 0.0
    WCAreaError = 0.0
    WCVertexError = 0.0
    for i in range(nx):
        for j in range(ny):
            area1 = poly1[i][j].Area()
            area2 = poly2[i][j].Area()
            area_error = abs(area1 - area2) / Area_0
            AveAreaError += area_error
            if area_error > WCAreaError:
                WCAreaError = area_error
            for k in range(poly1[i][j].npoints):
                x1 = poly1[i][j].pointlist[k].x
                y1 = poly1[i][j].pointlist[k].y
                x2 = poly2[i][j].pointlist[k].x
                y2 = poly2[i][j].pointlist[k].y
                vertex_error = (x1 - x2)**2 + (y1 - y2)**2
                AveVertexError += vertex_error
                if vertex_error > WCVertexError:
                    WCVertexError = vertex_error
    AveAreaError /= float(nx * ny)
    AveVertexError /= float(nx * ny * poly1[0][0].npoints)

    return [AveAreaError, WCAreaError, AveVertexError, WCVertexError]


#****************MAIN PROGRAM*****************

# First, read the .cfg files
run1 = sys.argv[1] # Path to data with charges in central pixel
configfile1 = "data/run"+run1+"/bf.cfg"
ConfigData1 = ReadConfigFile(configfile1)
filebase1 = ConfigData1["outputfilebase"]
filedir1 = ConfigData1["outputfiledir"]
Nx = ConfigData1["PixelBoundaryNx"]
Ny = ConfigData1["PixelBoundaryNy"]
PixelSize = ConfigData1["PixelSize"]
NxCenter = (Nx - 1) / 2
NyCenter = (Ny - 1) / 2
Area_0 = 100.0
NumVertices = ConfigData1["NumVertices"]
NumAngles = 4 * NumVertices + 4
NumElec = ConfigData1["CollectedCharge_0_0"]

run2 = sys.argv[2] # Path to data with charges to be modeled
configfile2 = "data/run"+run2+"/bf.cfg"
ConfigData2 = ReadConfigFile(configfile2)
filebase2 = ConfigData2["outputfilebase"]
filedir2 = ConfigData2["outputfiledir"]

# Create the polygon lists
poly0 = CreatePolyList(Nx, Ny) # Undistorted polygons
polyq = CreatePolyList(Nx, Ny) # Polygons with NumElec in central pixel
polyd = CreatePolyList(Nx, Ny) # Polygon deltas per electron
polym = CreatePolyList(Nx, Ny) # Modeled polygon shapes
polys = CreatePolyList(Nx, Ny) # Simulated polygon shapes

# Set up the undistorted polygons
FillUndistortedPolyList(ConfigData1, poly0)

# Read in the distorted polygons
filename1 = filedir1 + '/' + filebase1 +'_0_Vertices' + '.dat'
ReadVertexFile(filename1, ConfigData1, polyq)
filename2 = filedir2 + '/' + filebase2 +'_0_Vertices' + '.dat'
ReadVertexFile(filename2, ConfigData2, polys)

# Calculate the pixel shifts
FillPolyShiftsList(ConfigData1, poly0, polyq, polyd)

# Now read in the charges to be modeled
charges = zeros([Nx,Ny],dtype = int)

for i in range(ConfigData2["NumberofFilledWells_0"]):
    [x,y] = ConfigData2["FilledPixelCoords_0_%d"%i]
    charge = ConfigData2["CollectedCharge_0_%d"%i]
    chargei = int((x - ConfigData2["PixelBoundaryLowerLeft"][0]) / PixelSize)
    chargej = int((y - ConfigData2["PixelBoundaryLowerLeft"][1]) / PixelSize)
    charges[chargei, chargej] = charge

# Now build the modeled polygons
FillModelPolyList(ConfigData1, charges, poly0, polyd, polym)

# Now calculate the polygon errors
[AveAreaError, WCAreaError, AveVertexError, WCVertexError] = CalculatePixelErrors(ConfigData2, polym, polys)

# Now plot the results
PlotDelta = 2 # Number of pixels on either side of the center to be plotted
fig = figure()
fig.text(0.3, 0.30, "Average Area Error = %.4f percent"%(AveAreaError*100.0), fontsize = 10)
fig.text(0.3, 0.27, "Average Vertex Error = %.4f microns"%(AveVertexError), fontsize = 10)
fig.text(0.3, 0.24, "Worst Case Area Error = %.4f percent"%(WCAreaError*100.0), fontsize = 10)
fig.text(0.3, 0.21, "Worst Case Vertex Error = %.4f microns"%(WCVertexError), fontsize = 10)

ax1=axes([0.01,0.40,0.48,0.48],aspect=1)
PlotPixels(ax1, "Simulated Pixels", polys)
ax2=axes([0.50,0.40,0.48,0.48],aspect=1)
PlotPixels(ax2, "Modeled Pixels", polym)
savefig(filedir2+"/plots/PixelModel_SupTest.pdf")
