gdal2cesium
===========

Introduction
------------

A Python command line utility, based on GDAL and inspired by [gdal2tiles](http://www.gdal.org/gdal2tiles.html), to generate terrain heightmap tiles accordig to the native [Cesium format](http://cesiumjs.org/data-and-assets/terrain/formats/heightmap-1.0.html)

Any raster source supported by GDAL can be used. It can process multiple sources, with different size, resolution and covering area. The only requirement is that **sources must have the same data type** (Float32, Int16, etc.). You can transform the data type using [gdal_translate](http://www.gdal.org/gdal_translate.html).
The original data will be merged according to the following rules:

 - Lower resolution sources are used to generate lower zoom levels tiles
 - When two sources with different resolution overlaps, the lower one is used for lower zoom levels, the higher one is used for the highest zooms.
 - A tile for generated only if a source for the area covered by the tile exists and until source resolution is enough (less or equal) to the zoom level resolution.

This rules generate an optimal tiles coverage: tiles are generated only for those areas and zoom levels where a source meeting the required resolution is available. For each tile the data with the best resolution is choosen between the available sources.

Sources with different CRSs can be used but it's preferable and recomended to use data previously transformed to EPSG:4326 (WGS84). This will make the processing much more faster.

**NODATA** values inside the raster sources should be set to **0**.
An example of pre processing is:

<code>gdalwarp -s_srs epsg:3003 -t_srs epsg:4326 -srcnodata 9999 -dstnodata 0 source.tiff dest.tiff</code>

[gdalwarp](http://www.gdal.org/gdalwarp.html) will transform the source data to WGSS84 and set the destination NODATA value to 0.

Usage
-----

```
Usage: gdal2cesium.py [options] input_file(s)

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -s S_SRS, --s_srs=S_SRS
                        Define input raster CRS (eg EPSG:3003)
  -z ZOOM, --zoom=ZOOM  Zoom levels to render (format:'2-5' or '10').
  -r RESAMPLING, --resampling=RESAMPLING
                        Resampling method
                        (average,near,bilinear,cubic,cubicspline,lanczos) -
                        default 'average'
  -e, --resume          Resume mode. Generate only missing files.
  -v, --verbose         Print status messages to stdout
  -o OUTPUT, --o_dir=OUTPUT
                        Root output directory
  -i, --index           Create the shapefile of tiles index (True or False)
  -k, --keep            Keep temporary files reated by gdal2cesium
```

Example usage with tiff sources already in EPSG:4326:

<code>./gdal2cesium.py -o tiles *.tiff</code>

The ouput tiles directory will be saved into the /tiles subfolder
