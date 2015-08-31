#!/usr/bin/env python
#******************************************************************************
#  $Id: gdal2cesium.py 2014-10-01 12:01:23Z $
# 
# Project:  Cesium terrain generator for GDAL raster formats - S.I.T. Comune di Prato (Italy)
# Support:  Gis3w s.a.s. (http://gis3w.it)
# Purpose:  Convert a raster into a heightmap terrain for Cesium 3D Javascript library
#           - generate a global geodetic TMS data structure
#           - tiles are generated according to the Cesium heightmap binary format v1.0 (http://cesiumjs.org/data-and-assets/terrain/formats/heightmap-1.0.html)
#           - the max zoom level is calculated on the base of the raster horizontal resolution
#           - zoom levels up to the 0 level are always created to complete the parent-child relationships required by the Cesium format
# Author:   Giovanni Allegri (http://giovanniallegri.it, http://gis3w.it)
#
###############################################################################
# Copyright (c) 2014, S.I.T. Comune di Prato (Italy)

#  This program is free software; you can redistribute it and/or modify  
#  it under the terms of the GNU General Public License as published by  
#  the Free Software Foundation; either version 2 of the License, or     
#  (at your option) any later version.                                   
# 
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
#******************************************************************************

import os,sys,math,glob,struct,shutil
import subprocess

try:
    from osgeo import gdal
    from osgeo import osr
    from osgeo import gdal_array
except:
    import gdal
    print('You are using "old gen" bindings. gdal2cesium needs "new gen" bindings.')
    sys.exit(1)

from shapely.geometry import mapping, Polygon, LineString
from osgeo import ogr


try:
    import numpy
    import osgeo.gdal_array as gdalarray
except:
    print('gdal2cesium needs Numpy.')
    sys.exit(1)

MAXZOOMLEVEL = 32
resampling_list = ('average','near','bilinear','cubic','cubicspline','lanczos')

def makepoly(ulx,uly,lrx,lry):
    return Polygon([(ulx, uly), (lrx, uly), (lrx, lry), (ulx, lry), (ulx, uly)])
    
def makeline(ulx,uly,lrx,lry):
    return LineString([(ulx, uly), (lrx, uly), (lrx, lry), (ulx, lry), (ulx, uly)])
	
def splitpath(path):
    parts=[]
    (path, tail)=os.path.split( path)
    while path and tail:
         parts.append( tail)
         (path,tail)=os.path.split(path)
    parts.append( os.path.join(path,tail) )
    return map( os.path.normpath, parts)[::-1]

class GlobalGeodetic(object):
    def __init__(self, tileSize = 256):
        self.tileSize = tileSize

    def LatLonToPixels(self, lat, lon, zoom):
        "Converts lat/lon to pixel coordinates in given zoom of the EPSG:4326 pyramid"

        res = 180.0 / self.tileSize / 2**zoom
        px = (180 + lat) / res
        py = (90 + lon) / res
        return px, py

    def PixelsToTile(self, px, py):
        "Returns coordinates of the tile covering region in pixel coordinates"

        tx = int( math.ceil( px / float(self.tileSize) ) - 1 )
        ty = int( math.ceil( py / float(self.tileSize) ) - 1 )
        return tx, ty

    def LatLonToTile(self, lat, lon, zoom):
        "Returns the tile for zoom which covers given lat/lon coordinates"

        px, py = self.LatLonToPixels( lat, lon, zoom)
        return self.PixelsToTile(px,py)

    def Resolution(self, zoom ):
        "Resolution (arc/pixel) for given zoom level (measured at Equator)"

        return 180.0 / self.tileSize / 2**zoom
    def ZoomForPixelSize(self, pixelSize ):
        "Maximal scaledown zoom of the pyramid closest to the pixelSize."

        for i in range(MAXZOOMLEVEL):
            if pixelSize > self.Resolution(i):
                if i!=0:
                    return i-1
                else:
                    return 0 # We don't want to scale up

    def TileBounds(self, tx, ty, zoom):
        "Returns bounds of the given tile"
        res = 180.0 / self.tileSize / 2**zoom
        return (
            tx*self.tileSize*res - 180,
            ty*self.tileSize*res - 90,
            (tx+1)*self.tileSize*res - 180,
            (ty+1)*self.tileSize*res - 90
        )
        
    def TileBoundsForTileSize(self, tx, ty, zoom, extrapixels):
        res = 180.0 / self.tileSize / 2**zoom
        # we have to calculate a wider bound to consider the overlapping pixel according to Cesium format
        extrafactor = res*extrapixels
        return (
            tx*self.tileSize*res - 180,
            (ty*self.tileSize*res) - extrafactor - 90,
            ((tx+1)*self.tileSize*res - 180) + extrafactor,
            (ty+1)*self.tileSize*res - 90
        )

    def TileLatLonBounds(self, tx, ty, zoom):
        "Returns bounds of the given tile in the SWNE form"
        b = self.TileBounds(tx, ty, zoom)
        return (b[1],b[0],b[3],b[2])

import fnmatch
class PostProcessor(object):
    def __init__(self,outPath = "./",tmpOutPath = "./tmp"):
        self.rootPath = tmpOutPath
        self.pattern = '*.terrain'
        self.processedPath = outPath
        self.rtype = numpy.dtype('int16')
    
    def walk_tiles(self,folder = '.'):
        for root, _, files in os.walk(self.rootPath):
            for filename in fnmatch.filter(files, self.pattern):
                yield( os.path.join(root, filename))
        
    def get_tiles(self):
        terrains = []
        for terrain in self.walk_tiles():
            terrains.append(terrain)
        return terrains
    
    def extract_data(self,fin,rb):
        data = numpy.fromfile(fin,dtype=self.rtype,count=4096)
        data_mat = data.reshape((64,64))
        if rb == 'r':
            data_slice = data_mat[:,0] # left col
        else:
            data_slice = data_mat[0] # top row
        return data_slice
        
    def augment_tile(self,fin,data_slice_r,data_slice_b,pixel_br):
        data_complete = numpy.fromfile(fin,dtype=self.rtype,count=4097)
        data = data_complete[:4096]
        maskbytes = data_complete[4096]
        data_mat = data.reshape((64,64))
    
        data_slice_b_br = numpy.c_[[data_slice_b],[pixel_br]] # add bottom right pixel to bottom row as new col
        data_mat_r = numpy.c_[data_mat,data_slice_r] # add right col to data
        data_mat_brpr = numpy.r_[data_mat_r,data_slice_b_br] # add bottom row to data
        
        return data_mat_brpr,maskbytes
    
    def write_tile(self,tilename,new_tile,maskbytes):
        tilepath = os.path.join(self.processedPath,tilename)
        if not os.path.exists(os.path.dirname(tilepath)):
            os.makedirs(os.path.dirname(tilepath))
         
         
        tilearrayint = new_tile.astype(numpy.int16)
        data = tilearrayint.flatten()
        data_with_mask = numpy.append(data,maskbytes)
        data_with_mask.tofile(tilepath)
        
    def run(self):
        for terrain in self.get_tiles():
            pathparts = splitpath(terrain)
            idx = len(pathparts)
            root = os.path.join(*pathparts[:idx-3])
            y = int(pathparts[idx-1].split(".")[0])
            x = int(pathparts[idx-2])
            z = int(pathparts[idx-3])
            right_tile = os.path.join(root,str(z),str(x+1),"%s.terrain" % y)
            bottom_tile = os.path.join(root,str(z),str(x),"%s.terrain" % str(y-1))
            bottom_right_tile = os.path.join(root,str(z),str(x+1),"%s.terrain" % str(y-1))
            
            if os.path.exists(right_tile):
                with open(right_tile, 'rb') as right_tile_f:
                    data_slice_r = self.extract_data(right_tile_f,'r')
            else:
                data_slice_r = numpy.empty(64)
                data_slice_r.fill(5000)
            
            if os.path.exists(bottom_tile):
                with open(bottom_tile, 'rb') as bottom_tile_f:
                    data_slice_b = self.extract_data(bottom_tile_f,'t')
            else:
                data_slice_b = numpy.empty(64)
                data_slice_b.fill(5000)
            
            if os.path.exists(bottom_right_tile):
                with open(bottom_right_tile, 'rb') as bottom_right_tile_f:
                    data = numpy.fromfile(bottom_right_tile_f,dtype=self.rtype,count=1)
                    pixel_br = data[0]
            else:
                pixel_br = 5000
        
            with open(terrain, 'rb') as terrain_f:
                new_tile,maskbytes = self.augment_tile(terrain_f,data_slice_r,data_slice_b,pixel_br)
                tilename = os.path.join(*pathparts[idx-3:idx])
                self.write_tile(tilename,new_tile,maskbytes)

class GDAL2Cesium(object):
    # -------------------------------------------------------------------------    
    def process(self):
        for inumpyut_file in self.inumpyuts:
            self.inumpyut = inumpyut_file
            self.pre_process_inumpyut(inumpyut_file)
        
        self.merge_inumpyuts_data()
        
        self.make_tiles()
        
        print """Running post processing"""
        pp = PostProcessor(self.output,self.tmpoutput)
        pp.run()
        print """Post processing terminated"""
        shutil.rmtree(self.tmpoutput)
        
        
    def merge_inumpyuts_data(self):
        # Merge tminmax. We will use the extent containing all the layers for the lower zooms and only the higher resolution layer for the highest zooms
        global_tminmax = []
        for _inumpyut,inumpyut_data in self.inumpyuts_data.iteritems():
            #print "Inumpyut: %s" % _inumpyut
            minz = inumpyut_data[0]
            maxz = inumpyut_data[1]
            tminmax = inumpyut_data[2]
            for tz,tminmax_values in enumerate(tminmax):
                if (self.user_tminz is not None and tz < self.user_tminz) or (self.user_tmaxz is not None and tz > self.user_tmaxz):
                    continue
                if tz <= maxz:
                    #print "  tz: %s, tminmax: %s" % (tz,tminmax_values)
                    if len(global_tminmax)<=tz:
                        global_tminmax.append(list(tminmax_values))
                    else:
                        tminx = tminmax_values[0]
                        tminy = tminmax_values[1]
                        tmaxx = tminmax_values[2]
                        tmaxy = tminmax_values[3]
                        if tminx < global_tminmax[tz][0]:
                            global_tminmax[tz][0] = tminx
                        if tminy < global_tminmax[tz][1]:
                            global_tminmax[tz][1] = tminy
                        if tmaxx > global_tminmax[tz][2]:
                            global_tminmax[tz][2] = tmaxx
                        if tmaxy > global_tminmax[tz][3]:
                            global_tminmax[tz][3] = tmaxy
                            
        self.tminmax = global_tminmax
        
        # Split zooms in zoom ranges based on resolutions (to build the related vrt files)
        for _inumpyut,inumpyut_data in self.inumpyuts_data.iteritems():
            minz = inumpyut_data[0]
            maxz = inumpyut_data[1]
            if self.tminz is None or minz < self.tminz:
                self.tminz = minz
            if self.tmaxz is None or maxz > self.tmaxz:
                self.tmaxz = maxz
            for zoom in range(minz,maxz+1):
                if (self.user_tminz is not None and tz < self.user_tminz) or (self.user_tmaxz is not None and tz > self.user_tmaxz):
                    continue
                if self.zoom_resolutions.get(zoom) is None:
                    self.zoom_resolutions[zoom] = (inumpyut_data[3],inumpyut_data[4])
                else:
                    # the worst resolution is assigned to the common zoom levels (we check only resx, because resy will be consequently correlated)
                    if self.zoom_resolutions[zoom][0] < inumpyut_data[3]:
                        self.zoom_resolutions[zoom] = (inumpyut_data[3],inumpyut_data[4])
        '''print "MERGED"
        for tz,tminmax_values in enumerate(self.global_tminmax): 
            print "  tz: %s, tminmax: %s" % (tz,tminmax_values)
        '''
        
    # -------------------------------------------------------------------------
    def error(self, msg, details = "" ):
        """Print an error message and stop the processing"""
        if details:
            self.parser.error(msg + "\n\n" + details)
        else:
            self.parser.error(msg)
        exit(1)

    # -------------------------------------------------------------------------
    def progressbar(self, complete = 0.0):
        """Print progressbar for float value 0..1"""

        gdal.TermProgress_nocb(complete)

    # -------------------------------------------------------------------------
    def stop(self):
        """Stop the rendering immediately"""
        self.stopped = True

    # -------------------------------------------------------------------------
    def __init__(self, arguments ):
        """Constructor function - initialization"""
        try:
            subprocess.call(["gdalbuildvrt","--help"])
        except:
            print "gdalbuildvrt is required to run gdal2cesium in multi inumpyuts mode"
            exit(1)
        
        self.stopped = False
        self.multi_suffix = ''
        self.inumpyut = None
        self.default_base_output = 'tiles'
        self.min_tile_tz = None
        self.inumpyuts_data = {}
        self.inumpyuts_files_or_vrt = []
        self.vrts = {}
        self.tminmax = None
        self.zoom_resolutions = {}
        self.tminz = None
        self.tmaxz = None
        
        gdal.AllRegister()
        self.mem_drv = gdal.GetDriverByName( 'MEM' )
        self.geodetic = GlobalGeodetic()

        # Tile format
        self.tilesize = 64
        self.tileext = 'terrain'
        
        self.epsg4326 = "EPSG:4326"
        
        self.tilelayer = None

        self.scaledquery = True
        # How big should be query window be for scaling down
        # Later on reset according the chosen resampling algorightm
        self.querysize = 4 * self.tilesize
        
        # pixel overlap between tiles according to Ceiusm heightmap format
        self.extrapixels = 0 

        # RUN THE ARGUMENT PARSER:
        self.optparse_init()
        self.options, self.args = self.parser.parse_args(args=arguments)
        self.options.srcnodata = None
        if not self.args:
            self.error("No inumpyut file specified")

        # POSTPROCESSING OF PARSED ARGUMENTS:
        # Workaround for old versions of GDAL
        try:
            if (self.options.verbose and self.options.resampling == 'near') or gdal.TermProgress_nocb:
                pass
        except:
            self.error("This version of GDAL is not supported. Please upgrade to 1.6+.")
            #,"You can try run crippled version of gdal2tiles with parameters: -v -r 'near'")
        
        self.inumpyuts = [i for i in self.args]

        # Default values for not given options
        if self.options.output:
            self.output = self.options.output
        else:
            if len(self.inumpyuts)>0:
                self.multi_suffix = '_multi'
            self.output = os.path.join(self.default_base_output,os.path.basename( self.inumpyuts[0] ).split('.')[0]+self.multi_suffix)
            self.options.title = os.path.basename( self.inumpyuts[0]+self.multi_suffix )
        self.tmpoutput = os.path.join(self.output,'tmp')

        # Supported options
        self.resampling = None

        if self.options.resampling == 'average':
            try:
                if gdal.RegenerateOverview:
                    pass
            except:
                self.error("'average' resampling algorithm is not available.", "Please use -r 'near' argument or upgrade to newer version of GDAL.")
        elif self.options.resampling == 'near':
            self.resampling = gdal.GRA_NearestNeighbour
            self.querysize = self.tilesize
        elif self.options.resampling == 'bilinear':
            self.resampling = gdal.GRA_Bilinear
            self.querysize = self.tilesize * 2
        elif self.options.resampling == 'cubic':
            self.resampling = gdal.GRA_Cubic
        elif self.options.resampling == 'cubicspline':
            self.resampling = gdal.GRA_CubicSpline
        elif self.options.resampling == 'lanczos':
            self.resampling = gdal.GRA_Lanczos

        # User specified zoom levels
        self.user_tminz = None
        self.user_tmaxz = None
        if self.options.zoom:
            minmax = self.options.zoom.split('-',1)
            minmax.extend([''])
            min, max = minmax[:2]
            self.user_tminz = int(min)
            if max:
                self.user_tmaxz = int(max)
            else:
                self.user_tmaxz = int(min) 

        # Output the results
        if self.options.verbose:
            print("Options:", self.options)
            print("Inumpyut:", self.inumpyuts[0]+self.multi_suffix)
            print("Output:", self.output)
            print("Cache: %s MB" % (gdal.GetCacheMax() / 1024 / 1024))
            print('')

    # -------------------------------------------------------------------------
    def optparse_init(self):
        """Prepare the option parser for inumpyut (argv)"""
        from optparse import OptionParser, OptionGroup
        usage = "Usage: %prog [options] inumpyut_file(s)"
        p = OptionParser(usage, version="%prog ")
        
        p.add_option("-s", "--s_srs", dest="s_srs",
                          help="Define inumpyut raster CRS (eg EPSG:3003)")
        p.add_option('-z', '--zoom', dest="zoom",
                          help="Zoom levels to render (format:'2-5' or '10').")
        p.add_option("-r", "--resampling", dest="resampling", type='choice', choices=resampling_list,
                        help="Resampling method (%s) - default 'average'" % ",".join(resampling_list))
        p.add_option('-e', '--resume', dest="resume", action="store_true",
                          help="Resume mode. Generate only missing files.")        
        p.add_option("-v", "--verbose",
                          action="store_true", dest="verbose",
                          help="Print status messages to stdout")
        p.add_option("-o", "--o_dir",dest="output",
                          help="Root output directory")
        p.add_option("-i", "--index",dest="createtileindexshp",action="store_true",default=False,
                          help="Create the shapefile of tiles index (True or False)")
        p.add_option("-k", "--keep",dest="keepfiles",action="store_true",default=False,
                          help="Keep temporary files reated by gdal2cesium")
              
        p.set_defaults(resume=False,verbose=False,resampling='average')
        self.parser = p
                
    # -------------------------------------------------------------------------
    def pre_process_inumpyut(self,_inumpyut):
        """Initialization of the inumpyut raster, reprojection if necessary"""
        print "Processing: %s" % _inumpyut
        
        inumpyut_or_vrt = _inumpyut
        
        if not self.mem_drv:
            raise Exception("The 'MEM' driver was not found, is it available in this GDAL build?")

        # Open the inumpyut file
        if self.inumpyut:
            in_ds = gdal.Open(_inumpyut, gdal.GA_ReadOnly)
        else:
            raise Exception("No inumpyut file was specified")

        if self.options.verbose:
            print("Inumpyut file:", "( %sP x %sL - %s bands)" % (self.in_ds.RasterXSize, self.in_ds.RasterYSize, self.in_ds.RasterCount))

        if not in_ds:
            # Note: GDAL prints the ERROR message too
            self.error("It is not possible to open the inumpyut file '%s'." % _inumpyut )

        # Read metadata from the inumpyut file
        if in_ds.RasterCount == 0:
            self.error( "Inumpyut file '%s' has no raster band" % _inumpyut )

        if in_ds.GetRasterBand(1).GetRasterColorTable():
            # TODO: Process directly paletted dataset by generating VRT in memory
            self.error( "Please convert this file to RGB/RGBA and run gdal2tiles on the result.",
            """From paletted file you can create RGBA file (temp.vrt) by:
gdal_translate -of vrt -expand rgba %s temp.vrt
then run:
gdal2tiles temp.vrt""" % _inumpyut )

        # Get NODATA value
        in_nodata = []
        for i in range(1, in_ds.RasterCount+1):
            if in_ds.GetRasterBand(i).GetNoDataValue() != None:
                ndata = in_ds.GetRasterBand(i).GetNoDataValue()
                if math.isnan(ndata):
                    ndata = 'none'
                in_nodata.append( ndata )
        if self.options.srcnodata:
            nds = list(map( float, self.options.srcnodata.split(',')))
            if len(nds) < in_ds.RasterCount:
                in_nodata = (nds * in_ds.RasterCount)[:in_ds.RasterCount]
            else:
                in_nodata = nds

        if self.options.verbose:
            print("NODATA: %s" % in_nodata)

        #
        # Here we should have RGBA inumpyut dataset opened in in_ds
        #
        if self.options.verbose:
            print("Preprocessed file:", "( %sP x %sL - %s bands)" % (in_ds.RasterXSize, in_ds.RasterYSize, in_ds.RasterCount))

        # Spatial Reference System of the inumpyut raster
        self.in_srs = None

        if self.options.s_srs:
            self.in_srs = osr.SpatialReference()
            self.in_srs.SetFromUserInumpyut(self.options.s_srs)
            self.in_srs_wkt = self.in_srs.ExportToWkt()
        else:
            self.in_srs_wkt = in_ds.GetProjection()
            if not self.in_srs_wkt and in_ds.GetGCPCount() != 0:
                self.in_srs_wkt = in_ds.GetGCPProjection()
            if self.in_srs_wkt:
                self.in_srs = osr.SpatialReference()
                self.in_srs.ImportFromWkt(self.in_srs_wkt)

        # Spatial Reference System of tiles
        self.out_srs = osr.SpatialReference()
        self.out_srs.ImportFromEPSG(4326)

        # Are the reference systems the same? Reproject if necessary.
        out_ds = None
                              
        if (in_ds.GetGeoTransform() == (0.0, 1.0, 0.0, 0.0, 0.0, 1.0)) and (in_ds.GetGCPCount() == 0):
            self.error("There is no georeference - neither affine transformation (worldfile) nor GCPs. You can generate only 'raster' profile tiles.",
            "Either gdal2tiles with parameter -p 'raster' or use another GIS software for georeference e.g. gdal_transform -gcp / -a_ullr / -a_srs")
        
        in_srs_code = self.in_srs.GetAttrValue("AUTHORITY", 0)
        in_ds_srs = osr.SpatialReference()
        res = in_ds_srs.ImportFromWkt(in_ds.GetProjection())
        
        if res != 0 and in_srs_code is None:
                print "ERROR! The inumpyut file %s has no SRS associated and no SRS has been defined in inumpyut (-s parameter)" % _inumpyut
                exit(1)
 
        if self.in_srs:
            if in_ds_srs.ExportToProj4() != self.out_srs.ExportToProj4():
                if (self.in_srs.ExportToProj4() != self.out_srs.ExportToProj4()) or (in_ds.GetGCPCount() != 0):
                    print "WARNING! Inumpyut file %s has a SR different from EPSG:4326 (WGS84). This can make the processing significantly slow." % _inumpyut
                    # Generation of VRT dataset in tile projection, default 'nearest neighbour' warping
                    out_ds = gdal.AutoCreateWarpedVRT( in_ds, self.in_srs_wkt, self.out_srs.ExportToWkt() )

                    # TODO: HIGH PRIORITY: Correction of AutoCreateWarpedVRT according the max zoomlevel for correct direct warping!!!

                    if self.options.verbose:
                        print("Warping of the raster by AutoCreateWarpedVRT (result saved into 'tiles.vrt')")
                    out_ds.GetDriver().CreateCopy("%s.vrt" % _inumpyut, out_ds)
                    inumpyut_or_vrt = "%s.vrt" % _inumpyut

                    # Note: self.in_srs and self.in_srs_wkt contain still the non-warped reference system!!!

        else:
            self.error("Inumpyut file has unknown SRS.", "Use --s_srs ESPG:xyz (or similar) to provide source reference system." )

        if out_ds and self.options.verbose:
            print("Projected file:", "tiles.vrt", "( %sP x %sL - %s bands)" % (out_ds.RasterXSize, out_ds.RasterYSize, out_ds.RasterCount))
        
        if not out_ds:
            out_ds = in_ds

        #
        # Here we should have a raster (out_ds) in the correct Spatial Reference system
        #

        # Get alpha band (either directly or from NODATA value)
        alphaband = out_ds.GetRasterBand(1).GetMaskBand()
        if (alphaband.GetMaskFlags() & gdal.GMF_ALPHA) or out_ds.RasterCount==4 or out_ds.RasterCount==2:
            # TODO: Better test for alpha band in the dataset
            dataBandsCount = out_ds.RasterCount - 1
        else:
            dataBandsCount = out_ds.RasterCount
        
        # Read the georeference 
        out_gt = out_ds.GetGeoTransform()

        # Report error in case rotation/skew is in geotransform (possible only in 'raster' profile)
        if (out_gt[2], out_gt[4]) != (0,0):
            self.error("Georeference of the raster contains rotation or skew. Such raster is not supported. Please use gdalwarp first.")
            # TODO: Do the warping in this case automaticaly

        #
        # Here we expect: pixel is square, no rotation on the raster
        #

        # Output Bounds - coordinates in the output SRS
        ominx = out_gt[0]
        omaxx = out_gt[0]+out_ds.RasterXSize*out_gt[1]
        omaxy = out_gt[3]
        ominy = out_gt[3]-out_ds.RasterYSize*out_gt[1]
        # Note: maybe round(x, 14) to avoid the gdal_translate behaviour, when 0 becomes -1e-15

        if self.options.verbose:
            print("Bounds (output srs):", round(ominx, 13), ominy, omaxx, omaxy)

        #
        # Calculating ranges for tiles in different zoom levels
        #

        geodetic = GlobalGeodetic() # from globalmaptiles.py

        # Generate table with min max tile coordinates for all zoomlevels
        tminmax = list(range(0,32))
        for tz in range(0, 32):
            tminx, tminy = geodetic.LatLonToTile( ominx, ominy, tz )
            tmaxx, tmaxy = geodetic.LatLonToTile( omaxx, omaxy, tz )
            # crop tiles extending world limits (+-180,+-90)
            tminx, tminy = max(0, tminx), max(0, tminy)
            tmaxx, tmaxy = min(2**(tz+1)-1, tmaxx), min(2**tz-1, tmaxy)
            tminmax[tz] = (tminx, tminy, tmaxx, tmaxy)

        # Get the maximal zoom level (closest possible zoom level up on the resolution of raster)
        tminz = geodetic.ZoomForPixelSize( out_gt[1] * max( out_ds.RasterXSize, out_ds.RasterYSize) / float(self.tilesize) )
        if self.options.verbose:
            print ('Min Zoom: %s' % tminz)

        # Get the maximal zoom level (closest possible zoom level up on the resolution of raster)
        tmaxz = geodetic.ZoomForPixelSize( out_gt[1] )
        if self.options.verbose:
            print ('Max Zoom: %s' % tmaxz)
        
        self.inumpyuts_data[_inumpyut] = [tminz,tmaxz,tminmax,out_gt[1],out_gt[5]]
        
        self.inumpyuts_files_or_vrt.append(inumpyut_or_vrt)

        if self.options.verbose:
            print("Bounds (latlong):", ominx, ominy, omaxx, omaxy)

    def make_vrt(self,resx,resy,i):
        inumpyuts = " ".join(self.inumpyuts_files_or_vrt)
        if self.options.verbose:
            print "Building VRT file cesium_%s.vrt" % s
        try:
            res = subprocess.check_output("gdalbuildvrt -srcnodata 0 -resolution user -tr %s %s cesium_%s.vrt %s" % (abs(resx),abs(resy),i,inumpyuts), shell=True)
        except:
            exit(1)
    
    def make_tiles(self):
        # Generate the vrt files for zoom ranges
        i = 0
        tmp_res = -1
        vrt_file = None
        for tz in range(self.tminz,self.tmaxz+1):
            res = self.zoom_resolutions[tz][0]  # I check only with resx, because resy will be positively correlated
            if res != tmp_res:
		if i>0:
			self.vrts[vrt_file][1] = tz-1
                tmp_res = res
                resx = self.zoom_resolutions[tz][0]
                resy = self.zoom_resolutions[tz][1]
                self.make_vrt(resx,resy,i)
                vrt_file = "cesium_%s.vrt" % i
                self.vrts[vrt_file] = [tz,None]
                i += 1
            if tz == self.tmaxz:
                self.vrts[vrt_file][1] = tz
        
        self.ti_cum = 0
        if self.options.createtileindexshp and self.tilelayer is None:
            driver = ogr.GetDriverByName('Esri Shapefile')
            shptileindexfile = os.path.join(self.output,'tilesindex.shp')
            if os.path.exists(shptileindexfile):
                for f in glob.glob(self.output+'/tilesindex.*'):
                    os.remove(f)
            shptileindex = driver.CreateDataSource(shptileindexfile)
            self.tilelayer = shptileindex.CreateLayer('tiles', None, ogr.wkbLineString)
            self.tilelayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
            self.tilelayer.CreateField(ogr.FieldDefn('zoom', ogr.OFTInteger))
            self.tilelayer.CreateField(ogr.FieldDefn('tile', ogr.OFTString))
            self.tilelayer.CreateField(ogr.FieldDefn('children', ogr.OFTInteger))
            
        # Generate parent tiles
        self.generate_fake_parent_tiles()
        
        # For each vrt (i.e. zoom range) generate the tiles
        self.steps = len(self.vrts)
        self.step = 1
        for vrt in self.vrts.keys():
            self.process_vrt(vrt)
            if not self.options.keepfiles:
                try:
                    os.remove(vrt)
                except:
                    pass
            self.step += 1
        
        self.create_layerjsonfile()
        
        if self.options.createtileindexshp and self.tilelayer is not None:
            shptileindex.Destroy()
            shptileindex = self.tilelayer = feat = geom = None
            
        print """Processing finished. Tiles written to "%s".""" % self.output
    
    def process_vrt(self,vrt):
        self.open_inumpyut(vrt)
        self.generate_tiles(vrt)

        
    def open_inumpyut(self,vrt):
        if vrt:
            self.in_ds = gdal.Open(vrt, gdal.GA_ReadOnly)
        else:
            raise Exception("No vrt file was specified")
            
        if self.options.verbose:
            print("Inumpyut file:", "( %sP x %sL - %s bands)" % (self.in_ds.RasterXSize, self.in_ds.RasterYSize, self.in_ds.RasterCount))

        if not self.in_ds:
            # Note: GDAL prints the ERROR message too
            self.error("It is not possible to open the inumpyut file '%s'." % vrt )
            
        if self.in_ds.RasterCount == 0:
            self.error( "Inumpyut file '%s' has no raster band" % vrt )
        
        self.out_ds = self.in_ds
        
        # Get alpha band (either directly or from NODATA value)
        self.alphaband = self.out_ds.GetRasterBand(1).GetMaskBand()
        if (self.alphaband.GetMaskFlags() & gdal.GMF_ALPHA) or self.out_ds.RasterCount==4 or self.out_ds.RasterCount==2:
            self.dataBandsCount = self.out_ds.RasterCount - 1
        else:
            self.dataBandsCount = self.out_ds.RasterCount

    # -------------------------------------------------------------------------
    def make_child_flags(self,N,S,E,W):
        # Cesium format neighbor tiles flags
        HAS_SW = 0x01
        HAS_SE = 0x02
        HAS_NW = 0x04
        HAS_NE = 0x08
        
        NB_FLAGS = 0x00
        
        if N & W:
            NB_FLAGS = NB_FLAGS | HAS_NW
        if N & E:
            NB_FLAGS = NB_FLAGS | HAS_NE
        if S & W:
            NB_FLAGS = NB_FLAGS | HAS_SW
        if S & E:
            NB_FLAGS = NB_FLAGS | HAS_SE
        
        return NB_FLAGS
        
    def generate_fake_parent_tiles(self):
        tx = None
        for tz in range(self.tminz-1,-1,-1):
            tminx, tminy, tmaxx, tmaxy = self.tminmax[tz]
            tminx_c, tminy_c, tmaxx_c, tmaxy_c = self.tminmax[tz+1]
            for ty in range(tmaxy, tminy-1, -1):
                for tx in range(tminx, tmaxx+1):
                    tminx_cpot = tx * 2
                    tmaxx_cpot = tminx_cpot + 1
                    tminy_cpot = ty * 2
                    tmaxy_cpot = tminy_cpot + 1
                    
                    N = S = E = W = False
                    if tminx_cpot >= tminx_c:
                        W = True
                    if tmaxx_cpot <= tmaxx_c:
                        E = True
                    if tminy_cpot >= tminy_c:
                        S = True              
                    if tmaxy_cpot <= tmaxy_c:
                        N = True
                    
                    NB_FLAGS = self.make_child_flags(N,S,E,W)
                    if self.options.verbose:
                        print "Fake tile %s,%s,%s" % (tz,tx,ty)
                    self.write_fake_tile(tz,tx,ty,NB_FLAGS)
   
        # Write missing zero level tile with no children, tx 0 in case the zero level parent tileX is 1, 1 otherwise
        if tx:
            tx = 1-tx
        if tx is None:
            tx = 0

        self.write_fake_tile(0,tx,0,0x00)
            
    def write_fake_tile(self,tz,tx,ty,NB_FLAGS):
        tilefilename = os.path.join(self.tmpoutput, str(tz), str(tx), "%s.%s" % (ty, self.tileext))
        # Create directories for the tile
        if not os.path.exists(os.path.dirname(tilefilename)):
            os.makedirs(os.path.dirname(tilefilename))
            
        if self.options.createtileindexshp and self.tilelayer is not None:
            self.ti_cum += 1
            tilelayerdefn = self.tilelayer.GetLayerDefn()
            feat = ogr.Feature(tilelayerdefn)
            feat.SetField('id', self.ti_cum)
            feat.SetField('zoom', tz)
            feat.SetField('tile', "%s_%s_%s" % (tz, tx, ty))
            feat.SetField('children', NB_FLAGS)
            b = self.geodetic.TileBounds(tx, ty, tz)
            geom = ogr.CreateGeometryFromWkb(makeline(b[0], b[3], b[2], b[1]).wkb)
            feat.SetGeometry(geom)
            self.tilelayer.CreateFeature(feat)
            feat = geom = None
        
        # convert to integer representation of heightmap accordind to Cesium format and append children flags byte
        tilearrayint = (numpy.zeros(4096,numpy.dtype('int16')) + 1000) * 5
        tilearrayint.tofile(tilefilename)
        child_water_bytes = struct.pack('<BB',NB_FLAGS,0x00)
        with open(tilefilename,'ab') as outfile:
            outfile.write(child_water_bytes)
    
    def generate_tiles(self,vrt):
        """Generation of the Csium tiles from the inumpyut raster"""
        print("Generating Tiles (round %s of %s):" % (self.step,self.steps))

        # Cesium format neighbor tiles flags
        HAS_SW = 0x01
        HAS_SE = 0x02
        HAS_NW = 0x04
        HAS_NE = 0x08
        
        tminz = self.vrts[vrt][0]
        tmaxz = self.vrts[vrt][1]
        
        tcount = 0
        for tz in range(tmaxz, tminz-1, -1):
            tminx, tminy, tmaxx, tmaxy = self.tminmax[tz]
            tcount += (1+abs(tmaxx-tminx)) * (1+abs(tmaxy-tminy))
         
        ti = 0
        for tz in range(tmaxz, tminz-1, -1):
            # do not overwrite any real tile with successive inputs' fake tiles
            tminx, tminy, tmaxx, tmaxy = self.tminmax[tz]
            if tz < self.tmaxz:
                tminx_c, tminy_c, tmaxx_c, tmaxy_c = self.tminmax[tz+1]
            
            if self.options.verbose:
                tcount_zoom = (1+abs(tmaxx-tminx)) * (1+abs(tmaxy-tminy))
                print ("Tminx - Tmax: %s-%s" % (tminx,tmaxx))
                print ("Tminy - Tmaxy: %s-%s" % (tminy,tmaxy))
                print("Tile count for zoom %s: %s" % (tz,tcount_zoom))
            
            for ty in range(tmaxy, tminy-1, -1):
                for tx in range(tminx, tmaxx+1):
                    if self.options.resume and os.path.exists(os.path.join(self.tmpoutput, str(tz), str(tx), "%s.%s" % (ty, self.tileext))):
                        continue
                    # By Default the children flags are set to 0, which means no childern tiles (Cesium format)
                    NB_FLAGS = 0x00
                    # Child flags are calculated for all the zooms except the higher one (which has not children tiles)
                    if tz < self.tmaxz:
                        tminx_cpot = tx * 2
                        tmaxx_cpot = tminx_cpot + 1
                        tminy_cpot = ty * 2
                        tmaxy_cpot = tminy_cpot + 1
                        
                        N = S = E = W = False
                        if tminx_cpot >= tminx_c and tminx_cpot <= tmaxx_c:
                            W = True
                        if tmaxx_cpot >= tminx_c and tmaxx_cpot <= tmaxx_c:
                            E = True
                        if tminy_cpot >= tminy_c and tminy_cpot <= tmaxy_c:
                            S = True              
                        if tmaxy_cpot >= tminy_c and tmaxy_cpot <= tmaxy_c:
                            N = True
                            
                        NB_FLAGS = self.make_child_flags(N,S,E,W)
                    if self.stopped:
                        break
                    ti += 1
                    self.ti_cum += 1
                    tilearray  = self.process_tile(tz,tx,ty,ti,NB_FLAGS)
                    self.write_tile(tilearray,tz,tx,ty,NB_FLAGS)
                    if not self.options.verbose:
                        self.progressbar( ti / float(tcount) )


    def process_tile(self,tz,tx,ty,ti,NB_FLAGS):
        ds = self.out_ds
        tilebands = self.dataBandsCount
        querysize = self.querysize

        b = self.geodetic.TileBounds(tx, ty, tz)
        tilesize_aug = self.tilesize + self.extrapixels
        b_aug = self.geodetic.TileBoundsForTileSize(tx, ty, tz, self.extrapixels)
        if self.options.verbose:
            print "Tile bounds %s,%s,%s,%s" % (b[0], b[1], b[2], b[3])
            print "Tile bounds augomented %s,%s,%s,%s" % (b_aug[0], b_aug[1], b_aug[2], b_aug[3])
            #print "Tile poly: %s" % makepoly(b_aug[0], b_aug[1], b_aug[2], b_aug[3]).wkt

        if self.options.createtileindexshp and self.tilelayer is not None:
            '''
            shptileindex.write({
                'geometry': mapping(makepoly(b[0], b[3], b[2], b[1])),
                'properties': {'id': 123},
            })
            '''
            tilelayerdefn = self.tilelayer.GetLayerDefn()
            feat = ogr.Feature(tilelayerdefn)
            feat.SetField('id', self.ti_cum)
            feat.SetField('zoom', tz)
            feat.SetField('tile', "%s_%s_%s" % (tz, tx, ty))
            feat.SetField('children', NB_FLAGS)
            geom = ogr.CreateGeometryFromWkb(makeline(b[0], b[3], b[2], b[1]).wkb)
            feat.SetGeometry(geom)
            self.tilelayer.CreateFeature(feat)
            feat = geom = None
        rb, wb = self.geo_query( ds, b_aug[0], b_aug[3], b_aug[2], b_aug[1])
        nativesize = wb[0]+wb[2] # Pixel size in the raster covering query geo extent
        if self.options.verbose:
            print("\tNative Extent (querysize",nativesize,"): ", rb, wb)

        # Tile bounds in raster coordinates for ReadRaster query with extrapixels for Cesium tiles overlap
        querysize = self.querysize + ((self.querysize/self.tilesize) * self.extrapixels)
        rb, wb = self.geo_query( ds, b_aug[0], b_aug[3], b_aug[2], b_aug[1], querysize=querysize)

        rx, ry, rxsize, rysize = rb
        wx, wy, wxsize, wysize = wb
        if wxsize == 0:
            wxsize = 1
        if wysize == 0:
            wysize = 1

        if self.options.verbose:
            print("\tReadRaster Extent: ", (rx, ry, rxsize, rysize), (wx, wy, wxsize, wysize))

        # Query is in 'nearest neighbour' but can be bigger in then the tilesize
        # We scale down the query to the tilesize by supplied algorithm.
        # Tile dataset in memory
        dstile = self.mem_drv.Create('', tilesize_aug, tilesize_aug, tilebands, gdal.GDT_Float32)
        data = ds.ReadRaster(rx, ry, rxsize, rysize, wxsize, wysize, band_list=list(range(1,self.dataBandsCount+1)))
        datatype = gdal_array.GDALTypeCodeToNumericTypeCode(ds.GetRasterBand(1).DataType)
        if datatype != numpy.float32:
            data = numpy.frombuffer(data, dtype=datatype).astype(numpy.float32).tostring()
        
        if tilesize_aug == querysize:
            # Use the ReadRaster result directly in tiles ('nearest neighbour' query)
            dstile.WriteRaster(wx, wy, wxsize, wysize, data, band_list=list(range(1,self.dataBandsCount+1)))
        else:
            # Big ReadRaster query in memory scaled to the tilesize - all but 'near' algo
            dsquery = self.mem_drv.Create('', querysize, querysize, tilebands, gdal.GDT_Float32)
            dsquery.WriteRaster(wx, wy, wxsize, wysize, data, band_list=list(range(1,self.dataBandsCount+1)))
            self.scale_query_to_tile(dsquery, dstile)
            del dsquery
        del data
        
        tilearray = numpy.array(dstile.ReadAsArray())
        del dstile
        return tilearray
        #return None
    
    def write_tile(self,tilearray,tz,tx,ty,NB_FLAGS,WATER_MASK=0):
        tilefilename = os.path.join(self.tmpoutput, str(tz), str(tx), "%s.%s" % (ty, self.tileext))
        # Create directories for the tile
        if not os.path.exists(os.path.dirname(tilefilename)):
            os.makedirs(os.path.dirname(tilefilename))
        
        # convert to integer representation of heightmap accordind to Cesium format and append children flags byte       
        tilearray = (tilearray+1000) * 5
        tilearrayint = tilearray.astype(numpy.int16)
        tilearrayint.tofile(tilefilename)
        child_water_bytes = struct.pack('<BB',NB_FLAGS,WATER_MASK)
        with open(tilefilename,'ab') as outfile:
            outfile.write(child_water_bytes)
            
    def create_layerjsonfile(self):
        with open(os.path.join(self.output,'layer.json'),'w') as lj:
            lj.write("""{
              "tilejson": "2.1.0",
              "format": "heightmap-1.0",
              "version": "1.0.0",
              "scheme": "tms",
              "tiles": ["{z}/{x}/{y}.terrain"]
            }""")

    # -----------------------------------------------------------------------
    def geo_query(self, ds, ulx, uly, lrx, lry, querysize = 0):
        """For given dataset and query in cartographic coordinates
        returns parameters for ReadRaster() in raster coordinates and
        x/y shifts (for border tiles). If the querysize is not given, the
        extent is returned in the native resolution of dataset ds."""

        geotran = ds.GetGeoTransform()
        rx= int((ulx - geotran[0]) / geotran[1] + 0.001)
        ry= int((uly - geotran[3]) / geotran[5] + 0.001)
        rxsize= int((lrx - ulx) / geotran[1] + 0.5)
        rysize= int((lry - uly) / geotran[5] + 0.5)

        if not querysize:
            wxsize, wysize = rxsize, rysize
        else:
            wxsize, wysize = querysize, querysize

        # Coordinates should not go out of the bounds of the raster
        wx = 0
        if rx < 0:
            rxshift = abs(rx)
            wx = int( wxsize * (float(rxshift) / rxsize) )
            wxsize = wxsize - wx
            rxsize = rxsize - int( rxsize * (float(rxshift) / rxsize) )
            rx = 0
        if rx+rxsize > ds.RasterXSize:
            wxsize = int( wxsize * (float(ds.RasterXSize - rx) / rxsize) )
            rxsize = ds.RasterXSize - rx

        wy = 0
        if ry < 0:
            ryshift = abs(ry)
            wy = int( wysize * (float(ryshift) / rysize) )
            wysize = wysize - wy
            rysize = rysize - int( rysize * (float(ryshift) / rysize) )
            ry = 0
        if ry+rysize > ds.RasterYSize:
            wysize = int( wysize * (float(ds.RasterYSize - ry) / rysize) )
            rysize = ds.RasterYSize - ry

        return (rx, ry, rxsize, rysize), (wx, wy, wxsize, wysize)

    # -------------------------------------------------------------------------
    def scale_query_to_tile(self, dsquery, dstile):
        """Scales down query dataset to the tile dataset"""

        querysize = dsquery.RasterXSize
        tilesize = dstile.RasterXSize
        tilebands = dstile.RasterCount

        if self.options.resampling == 'average':
            for i in range(1,tilebands+1):
                res = gdal.RegenerateOverview( dsquery.GetRasterBand(i),
                    dstile.GetRasterBand(i), 'average' )
                if res != 0:
                    self.error("RegenerateOverview() failed")
        else:
            # Other algorithms are implemented by gdal.ReprojectImage().
            dsquery.SetGeoTransform( (0.0, tilesize / float(querysize), 0.0, 0.0, 0.0, tilesize / float(querysize)) )
            dstile.SetGeoTransform( (0.0, 1.0, 0.0, 0.0, 0.0, 1.0) )
            
            res = gdal.ReprojectImage(dsquery, dstile, None, None, self.resampling)    
            if res != 0:
                self.error("ReprojectImage() failed on %s, error %d" % (tilefilename, res))

if __name__=='__main__':
    argv = gdal.GeneralCmdLineProcessor( sys.argv )
    if argv:
        gdal2cesium = GDAL2Cesium( argv[1:] )
        gdal2cesium.process()
