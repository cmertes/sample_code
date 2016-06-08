import database_util as db
import sample_util as sm
import image_util as im 
import file_util as df
import numpy as np
from gdalconst import *
import gdal
import shutil
import ogr 
import os
import math
import sys
import collections
import glob
import csv
from final_conf_global_var import *


'''
create error probability surfaces from boosting algorithm. 
uses average confidence over multiple classifications for
landsat overlap segments.

input: class confidences for each landsat footprint and remap files.
output: statewide prediction probability for each land cover class.
'''

''' ______________________________________________ user defined variables'''

# classification level
level = 'level1'

# working/root directory
wd = 'S:/classification/maps/final_footprints/level1/'

# skip compositing for these footprints
skip_pr=['26027','23028']

# input directory containing probability surfaces
dir_control = wd + 'control/' + level +'/' 

# output directory for reprojected wtm files                                              
dir_wtm = wd + 'wtm/' + level +'/' 

# output directory for clipped subsets                                                             
dir_subsets = wd + 'subsets/' + level +'/'     

# output directory for composited subsets                                                               
dir_final_rasters = wd + 'final_rasters/' + level +'/' 

# output directory for final mosaics                                       
dir_final_mosaics = wd + 'final_mosaics/' + level + '/'            

'get subset_ids from database'
subsets = sm.getColumn(dbase,'subset_id','mosaic.mosaic_final')          

'get prs from database'
prs = sm.getColumn(dbase,'pr','boundary.boundary_wi_footprints')

'''_______________________________________________________ functions'''


def reprojectFiles(dirn, branch):
    '''reproject confidence file to wtm 3071'''
    
    'get all the maps and store them in list'
    file_list = glob.glob(dirn + '/*conf.tif')
                          
    for f in file_list:
        print f
        # pr from file name
        pr = f.split('/')[-1][:5]
        
        # utm zone from database
        source_epsg = im.getEPSG(pr,dbase)

        # use use def module to reproject
        out_dir = dir_wtm + branch + '/'
        im.reprojectRaster(f,pr,'30',source_epsg,out_dir,'conf') 


def bandCount(raster):
    '''Return count of bands in raster file'''
                          
    in_raster = gdal.Open(raster)
    bands = in_raster.RasterCount

    return bands


def createsubsetDir(label):
    '''create the subset directories to hold the clipped rasters'''
    
    # get subset_id column from database                          
    query = "SELECT subset_id FROM mosaic.mosaic_final"
    rows = db.fetch(query,dbase)
    for row in rows:

      # create directory 
      subset_dir = dir_subsets + "/" + label + "/subset" + str(row[0])
      df.createDir(subset_dir)

          
def getBoundingBoxList(subset_id):
    '''get bounding box to clip raster'''
    
    # get the image extent from database
    extent_query = "select st_extent(geom) as extent from mosaic.mosaic_final where subset_id = " + subset_id
    rows = db.fetch(extent_query,dbase)
    
    # process the rows object to get 4 seperate values from the box and store in a list
    for row in rows:

        # replace characters in the string
        box = row[0][4:-1].replace(",", " ")

        # split the string at space to separate coordinates
        cornerlist = box.split(' ')
                          
        return cornerlist        


def clipRaster(cond, label):
    '''Clip the map to the subsets that fall inside of it.'''

    # get list of maps that will be clipped from directories
    maplist = df.getFileListbywd(dir_wtm + '/' + label,[cond])
    
    # iterate throught the maps in maplist and clip to subset geometry if it falls inside of raster extent'
    for m in maplist:
        pr = m.split('/')[-1][2:7]
        subsetlist = db.fetch("SELECT subset_id FROM mosaic.mosaic_final where pr" + str(pr) + " = 'true'", dbase)

        for subset in subsetlist:
            print 'CLIPPING ', pr, 'TO subset: ', subset[0], '\n'
            subset_tif = dir_subsets + '/' + label + '\subset' + str(subset[0])  + '/pr' + str(pr) + '_subset' + str(subset[0]) + cond
            cornerlist = getBoundingBoxList(str(subset[0]))
            clip = 'gdalwarp -te ' + cornerlist[0] + ' ' + cornerlist[1] + ' ' + cornerlist[2] + ' ' + cornerlist[3] + \
                   ' -s_srs EPSG:3071 -cutline "PG:host=144.92.235.143 user=wiscland dbname=samples password=landcover" \
                    -csql "select * from mosaic.mosaic_final where subset_id = ' + str(subset[0]) + '" -dstnodata 255 ' + m + ' ' + subset_tif

            os.system(clip)
                         
                
def stackImages(suffix,label):
    '''stack images to multiband'''

    if not os.path.exists(dir_finalRasters + label): os.mkdir(dir_finalRasters + label)
    
    for subset in subsets:

        subset_dir = dir_subsets + '/' + label + '/subset' + str(subset[0]) + '/'
        subsetlist = [x for x in glob.glob(subset_dir + '*conf.tif')]

        if subsetlist:
            
            print 'ALIGNING IMAGES FOR: ', subset[0]
            if len(subsetlist)>1:
                
                for s in skip_pr: #files to skip compositing
                        skip=[m for m in subsetlist if s in m]
                        for path in skip:
                                os.remove(path)
            
            '''count again after removed'''
            subsetlist = [x for x in glob.glob(subset_dir + '*conf.tif')] 
                                                    
            '''move file to final raster location if no overlap in landsat scenes'''
            if len(subsetlist) == 1:
                print 'subset ', subset, 'already finalized.'
                in_file = subsetlist[0]
                out_file = dir_finalRasters + '/' + label + '/subset' + str(subset[0]) + '.tif'

                map_in = gdal.Open(in_file)
                rows = map_in.RasterYSize
                cols = map_in.RasterXSize
                bands = map_in.RasterCount
                out_arr = map_in.ReadAsArray() 
                
                'output raster'
                driver = map_in.GetDriver()
                outDs = driver.Create(out_file, cols, rows, 1, GDT_Byte)
                outBand = outDs.GetRasterBand(1)
                outBand.WriteArray(out_arr, 0, 0)
                outBand.FlushCache()
                outBand.SetNoDataValue(255)
          
                outDs.SetGeoTransform(map_in.GetGeoTransform())
                outDs.SetProjection(map_in.GetProjection())
                
                del map_in
                del outDs
            
            elif len(subsetlist) > 1:

                stacklist=[x for x in glob.glob(subset_dir + '*conf.tif')] 
                
                'name of the stacked image that will go into gdal_merge output argument'    
                stackimage = subset_dir + '/subset' + str(subset[0]) + 'stack.tif' 
                      
                merge = "gdal_merge.py -o " + stackimage + " -of GTiff -separate " + ' '.join(stacklist)
                os.system(merge)


def createAveragedRasters(branch):           
    '''use averaging technique to composite raster overlap areas'''
    
    for subset in subsets:
        
        subset_dir = dir_subsets + '/' + branch + '/subset' + str(subset[0]) + '/'
        out_file = dir_final_rasters + '/' + branch + '/' + 'subset' + str(subset[0]) + '.tif'

        if glob.glob(subset_dir + '*stack.tif'):

            conf = glob.glob(subset_dir + '*stack.tif')[0]
            map_in = gdal.Open(conf)
            rows = map_in.RasterYSize
            cols = map_in.RasterXSize
            in_arr = map_in.ReadAsArray()
            in_arr = in_arr.astype(float)
            in_arr[in_arr==0]=np.nan   #convert missing conf to nan
            
            'output raster'
            driver = map_in.GetDriver()
            outDs = driver.Create(out_file, cols, rows, 1, GDT_Byte)
            outBand = outDs.GetRasterBand(1)
            out_arr=np.rint(np.nanmean(in_arr, axis=0)).astype(int)
            outBand.WriteArray(out_arr, 0, 0)
            outBand.FlushCache()
            outBand.SetNoDataValue(255)
      
            outDs.SetGeoTransform(map_in.GetGeoTransform())
            outDs.SetProjection(map_in.GetProjection())

            del map_in
        
        
def remapDict(lookup):
        '''Read remap file of dict values'''

        dic = df.getDictfromCSV(lookup,'\t',1,0)

        return dic


def validLabels(branch):
    '''return valid children on classification branch'''

    if level=='level1':
        query="SELECT DISTINCT " + level + """
          FROM classify.holder_levels_c5
          WHERE level1<>2000 AND level1<>9999"""
        valid_labels=[str(x[0])for x in db.fetch(query,dbase)]

    else:
        query="SELECT DISTINCT " + level + """
          FROM classify.holder_levels_c5 
          WHERE """ + level + "<>9999 AND level" + str(int(level[-1])-1) + "='" + str(branch) + "' AND level" + str(int(level[-1])-1) + "<>" + level
        valid_labels=[str(x[0]) for x in db.fetch(query,dbase)]

    return sorted(valid_labels)

def branchName():
    '''return branch'''

    if level=='level1':
        branch='1000'
    else:
        fname=sorted(glob.glob(dir_control+'/*map.tif'))[0]
        branch=fname.split('_')[-2]
        
    return branch
        
def mosaicRasters(out_path, branch):
    '''Mosaics images within directory into single band image'''
    
    subsetlist = [x for x in glob.glob(dir_final_rasters + '/' + branch + '/*.tif')]
    out_file = out_path + 'wiscland2_' + level + '_' + branch + '_conf.tif'

    '''loop through the subsetList and append each of its elements into a string for input into gdal command'''
    filelist = ' '.join(subsetlist)
               
    task = "gdal_merge.py -o " + out_file + " -ps 30 30 -of GTiff -n 255 -a_nodata 255 " + filelist
    print '>>>FINAL STATEWIDE MOSAIC'
    os.system(task)

    edit = "gdal_edit.py -a_srs EPSG:3071 " + out_file
    os.system(edit)

def fileMeta(dirname, pattern):
    '''count files in dirname with pattern, return prs and file stem format'''

    files = glob.glob(dirname + '/*' + pattern + '*')
    prs = [f.split('/')[-1][3:8] for f in files]
    dbtem = files[0]
	
    return len(files), prs, dbtem

def getFootprint():
    '''get list of all prs from db'''
    
    query = "SELECT pr FROM boundary.boundary_wi_footprints"
    prs = [str(x[0]) for x in db.fetch(query, dbase)]
	
    return prs
                 
def fillMap(pr, out_file):
    '''create fill maps for footprints without confidences'''

    copy_map = glob.glob(dir_control[:-2] + '1/pr0' + str(pr) + '*confidence.tif')[0]
    copy_in = gdal.Open(copy_map)
    copy_arr = copy_in.ReadAsArray()

    rows = copy_in.RasterYSize
    cols = copy_in.RasterXSize
   
    out_driver = copy_in.GetDriver()
    out_ds = out_driver.Create(out_file, cols, rows, 1, gdal.GDT_Byte)
    out_ds.SetGeoTransform(copy_in.GetGeoTransform())
    out_ds.SetProjection(copy_in.GetProjection())
    
    out_band = out_ds.GetRasterBand(1)
    fill_arr = np.where(copy_arr[1]<>255, 0, 255)
    out_band.WriteArray(fill_arr, 0, 0)
    out_ds.SetGeoTransform(copy_in.GetGeoTransform())
    out_ds.SetProjection(copy_in.GetProjection())


def getExtent(filename):
    '''get image extent in native coordinates'''

    '''Store dataset information into variables and/or above lists'''
    dataset = gdal.Open(filename, GA_ReadOnly)
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize
    projection = dataset.GetProjection()  
    geotransform = dataset.GetGeoTransform()


    xmin = geotransform[0]
    xmax = geotransform[0] + (30 * cols)
    ymin = geotransform[3] - (30 * rows)
    ymax = geotransform[3]
 
    xres = abs(geotransform[1])
    yres = abs(geotransform[5])

    return xmin,xmax,ymin,ymax,xres,yres

    
def reprojectRaster(infile, outfile, epsg):
    '''repoject raster'''

    warp = 'gdalwarp -overwrite -t_srs EPSG:' + epsg + ' -of GTiff -co COMPRESS=LZW ' + infile + " " + outfile
    os.system(warp)


def getEPSG():
    '''retrieve epsg of footprint from table in db'''

    code_dict = dict()
    query =  "SELECT DISTINCT pr, utm_epsg FROM classify.utm_zones"
    codes= db.fetch(query,dbase)

    for c in codes:
        code_dict[str(c[0])] = str(c[1])

    return code_dict


def extractBand(in_file, band, out_file):
    '''extract band from multiband raster'''

    ds_in = gdal.Open(in_file)
    rows = ds_in.RasterYSize
    cols = ds_in.RasterXSize
    out_driver = ds_in.GetDriver()
    out_ds = out_driver.Create(out_file, cols, rows, 1, gdal.GDT_Byte)
    out_ds.SetGeoTransform(ds_in.GetGeoTransform())
    out_ds.SetProjection(ds_in.GetProjection())
    
    im_arr = ds_in.GetRasterBand(band).ReadAsArray()
    band = out_ds.GetRasterBand(1)
    band.WriteArray(im_arr)
    band.FlushCache()

    ds_in = ds_out = None

  
'''__________________________________________________________________________________________'''


subdirs = df.getSubdirs(dir_control)

for sd in subdirs:

    class_branch = sd.split('/')[-1]
    labels = validLabels(class_branch)
    fcount, prs, dbtem = fileMeta(sd, 'confidence.tif')
    allpr = getFootprint()

    for label in labels:
        out_path = sd + '/' + str(label)
        if not os.path.exists(out_path): os.mkdir(out_path)
        
    # extract each class confidence from multiband images

    for pr in prs:
        remap = remapDict(glob.glob(sd + '/pr0' + pr + '*remap_values.txt')[0])
        conf_file = glob.glob(sd + '/pr0' + pr + '*confidence.tif')[0]

        # all classes with predictions
        bands = [x for x in remap.iteritems() if x[1] in labels]
        values = sorted(int(x[0]) for x in bands)
        value_dic={}

        # determine band order from c5
        for i in range(1,len(values)+1):
            value_dic[str(values[i-1])]=i 
                
        for band in bands:
            out_file = sd + '/' + band[1] + '/' + pr + '_conf.tif'
            extractBand(conf_file, value_dic[band[0]], out_file)


        # all classes without predictions
        fill = [c for c in labels if c not in remap.itervalues()]

        for f in fill:
            out_file = sd + '/' + f + '/' + pr + '_conf.tif'
            print 'FILLING VALUE FOR {0:s} : {1:s}'.format(pr, f)
            fillMap(pr, out_file)

    # pr without branch classification
    missing_pr = [p for p in allpr if p not in prs]
    for mp in missing_pr:
        for label in labels:
            out_file = sd + '/' + label + '/' + mp + '_conf.tif'
            print 'FILLING VALUE FOR {0:s} : {1:s}'.format(mp, label)
            fillMap(mp, out_file)

    for label in labels:
        out_path = sd + '/' + str(label)
        if not os.path.exists(out_path): os.mkdir(out_path)
        reprojectFiles(out_path, label)
                          
        createsubsetDir(label)
        clipRaster('conf.tif', label)
        stackImages(['conf.tif'], label)
        createAveragedRasters(label)

        out_path = dir_final_mosaics + label + '/'
        if not os.path.exists(out_path): os.mkdir(out_path)
        mosaicRasters(out_path, label)

        


 
 

