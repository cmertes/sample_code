from osgeo import ogr
import numpy as np
from numpy import random, dstack
import random
import os
import sys
import csv
import subprocess
import glob
import gdal
from gdalconst import *
import dir_util as df
from apply_segments_global_var import *


'''
integrates naip image segmentation results with landsat land cover 
classification. uses simple majority rule on land cover over segment boundaries.
Segmented map is mosaicked, sieved, and remapped to final labels.

input: pixel-based classification results from c5 and segments (vector)
output: segmented map with 4-digit land cover codes
'''

wd='C:/classification/segmentation/segments/'
map_file='C:/classification/mosaicking/final_mosaics/level1/statewide_level1.tif'
out_file='C:/classification/maps/statewide/level1/statewide_level1_map.tif'
remap_file='C:/classification/maps/final_footprint/level1/pr026029_1305b_20151125_remap_values.txt'

'''_____________________________________________________________DEFINE FUNCTIONS'''

def resampleMap(map_file, divisor):
    '''downsample map resolution by divisor'''

    outmap = map_file.split('.')[0]+'_rs.tif'
    in_map = gdal.Open(map_file)
    geotransform = in_map.GetGeoTransform()
    x_res = str(geotransform[1] / divisor)
    y_res = str(geotransform[5] / divisor)
    del inmap

    resample = 'gdalwarp -tr ' + xres + ' ' + yres + ' ' + map_file + ' ' + outmap
    print 'resampling to x: {0:s} y:{1:s} '.format(x_res, y_res)
    os.system(resample)
    return outmap


def gme_isectpolyrst(polygon, raster):
        '''call GME isectpolyrst to intersect segments with map'''

        file_name = wd+'tmp.txt'
        newFile = 'isectpolyrst(in="' + polygon + '", raster="' + raster + '", prefix="lev", thematic=TRUE, proportion=TRUE, metrics="CNT", allowpartialoverlap=TRUE, medquant=FALSE);'
        newFileObj = open(file_name, 'w')
        newFileObj.write(newFile)
        newFileObj.close()
            
        print "Intersecting polygons with raster..."
        subprocess.call([r'C:\Program Files (x86)\SpatialEcology\GME\SEGME.exe', '-c', 'run(in="' + file_name + '");'])
        os.remove(file_name)

def gme_ptinpoly(polygon, point, where):
    '''call GME genpointinpoly to generate points in unlabeled polygons'''

    file_name = wd+'tmp.txt'
    newFile = 'genpointinpoly(in="' + polygon + '", out="' + point + '", position="LABEL", copyfields=TRUE, ' + where + ');'
    newFileObj = open(file_name, 'w')
    newFileObj.write(newFile)
    newFileObj.close()
            
    print "Points in polygons..."
    subprocess.call([r'C:\Program Files (x86)\SpatialEcology\GME\SEGME.exe', '-c', 'run(in="' + file_name + '");'])
    os.remove(file_name)
    
def gme_isectptrst(points, raster):
    '''intersect point with raster with gme'''
    file_name = wd+'tmp.txt'
    command = 'isectpntrst(in="' + points + '", raster="' + raster + '", field="label", update=TRUE);'
    newFileObj = open(file_name, 'w')
    newFileObj.write(command)
    newFileObj.close()
    
    print "Intersect points with raster..."
    subprocess.call([r'C:\Program Files (x86)\SpatialEcology\GME\SEGME.exe', '-c', 'run(in="' + file_name + '");'])
    os.remove(file_name)

def nearestPoint(points, compare_pt):
    '''Get nearest point'''

    file_name = wd+'tmp.txt'
    outf = points.split('.')[0]+'dist.csv'
    nn = points.split('.')[0]+'nn.csv'
    command = 'pointdistances(in="' + points + '", fld="RefID", out="'+outf+'", in2="' + compare_pt + '", fld2="RefId", nearest=1, nout="'+nn+'");'
    print command
    newFileObj = open(file_name, 'w')
    newFileObj.write(command)
    newFileObj.close()
    
    subprocess.call([r'C:\Program Files (x86)\SpatialEcology\GME\SEGME.exe', '-c', 'run(in="' + file_name + '");'])
    os.remove(file_name)

    return nn

def countFeatures(shapefile):
    '''count features in shapefile'''

    ds = ogr.Open(shapefile, 0)
    layer = ds.GetLayer()
    num_feat = layer.GetFeatureCount()
    
    del ds
    layer=None
    
    return num_feat

def nearestLabel(point, compare_pt, nn_file, dic):
    '''transcribe nearest label between points from csv of nearest neighbors'''
    
    ds = ogr.Open(point, update=0)
    pointf = ds.GetLayer()
    sql = 'SELECT RefID FROM ' +point.split('/')[-1].split('.')[0]+' WHERE label=254'
    replacef = ds.ExecuteSQL(sql)

    replace=[]
    for x,y in enumerate(replacef):
        replace.append(y.GetField(0))
    print 'replace', replace

    if replace:
        nns=[]                              #nearest neighbors for all points
        with open(nn_file, 'rU') as infile:
            in_nn=csv.reader(infile, delimiter=',')
            for row in in_nn:
                nns.append(row)
        del ds
        pointf=None

        ds=ogr.Open(point, update=0)
        layer=ds.GetLayer()
        num_feat=layer.GetFeatureCount()

        ds2=ogr.Open(allpt, update=0)
        layer2=ds2.GetLayer()
        
        for i in range(num_feat):
            feature=layer.GetFeature(i)
            fid=feature.GetFieldAsString("RefID")

            if int(fid) in replace:
                nn=[x[1] for x in nns if x[0]==fid]
                print nn
                nsql='SELECT RefId, label FROM ' + compare_pt.split('/')[-1].split('.')[0] + ' WHERE RefId='+nn[0]
                nn_label=ds2.ExecuteSQL(nsql)
                for i, j in enumerate(nn_label):
                    fid=int(fid)
                    dic[fid]=int(j.GetField(1))

    return dic
               
def update_dict(polygon, point):
    '''get update label values from point file'''

    update_dict={}
    ds = ogr.Open(point, update=0)
    pointf = ds.GetLayer()
    sql = 'SELECT RefID, label FROM '+ point.split('/')[-1].split('.')[0]+' WHERE label<255'
    att = ds.ExecuteSQL(sql)
    if att:
        for i, j in enumerate(att):
            update_dict[j.GetField(0)]=int(j.GetField(1))
    
    return update_dict

def updateLabel(layer, dic):
    '''Update shapefile 'label' field according to points'''
    
    keys=dic.keys()
    feature = layer.GetNextFeature()
    keys=[x for x in dic.keys()]
    while feature:
        refid = feature.GetFieldAsInteger("ID")
            
        if refid in keys:
            feature.SetField("label", dic[refid])
            layer.SetFeature(feature)
        
        feature = layer.GetNextFeature()
    print len(keys), ' features updated.'
        
def gid(layer, fieldname):
    '''serial another field'''

    num_features = layer.GetFeatureCount()
    print num_features
    for i in range(num_features):
            feature = layer.GetFeature(i)
            if feature:
                feature.SetField(fieldname, i)
                layer.SetFeature(feature)
   
def calcMajority(layer, classes):
        '''calcualte majority over values in column# start to column# end'''
        
        num_features = layer.GetFeatureCount()
        for i in range(num_features):
                feature = layer.GetFeature(i)
                zonal = [feature.GetField(x) for x in classes]
                max_cnt = max(zonal) 
                labels=[a for a,b in enumerate(zonal) if b==max_cnt]
        
                #resolve ties
                if len(labels)==1:
                        label=labels[0]
                elif len(labels)>1 and max_cnt>0:
                        label=labels[random.randrange(0,len(labels))]   #if tied, take random
                else: label=255                                         #no data
                feature.SetField('label', label)
                layer.SetFeature(feature)

        source = None
        print 'shapefile updated'

def overruleMaj(grids, mask):
    '''overrule majority if segment intersects another raster'''
    
    
def getFieldNames(layer):
        '''get field names from shapefile layer'''
        
        layer_defn = layer.GetLayerDefn()
        field_names = [layer_defn.GetFieldDefn(i).GetName() for i in range(layer_defn.GetFieldCount())]

        return field_names

def dissolvePolys(out_file, infile, fieldname):
        '''dissolve polygons on fieldname'''
        
        inlayer = infile.split('/')[-1].split('.')[0]
        dissolve = "ogr2ogr " + out_file + " " + infile + ' -dialect sqlite -sql \
        "SELECT ST_UNION(geometry), ' + fieldname + " FROM '" + inlayer + "' GROUP BY " + fieldname + '"'
        print 'Dissolving...'
        os.system(dissolve)

def addIntField(layer, fieldname):
        '''add integer type field to shapefile'''
        
        print "adding new field, '", fieldname, "'"
        new_field = ogr.FieldDefn(fieldname, ogr.OFTInteger)
        layer.CreateField(new_field)

def rasterizeSeg(infile, out_file, fieldname):
        '''Rasterize shapefile to 30m raster'''
        
        rasterize = 'gdal_rasterize -a ' + fieldname + ' -tr 30 30 -ot Byte -a_nodata 255 ' + infile + ' ' + out_file
        print 'Rasterizing segments...'
        os.system(rasterize)

def mosaicGrids(out_file, infiles):
        '''mosaic rasterized segment grids to statewide'''
        
        mosaic = 'gdal_merge.py -o ' + out_file + ' -of GTiff -n 255 -a_nodata 255 ' + ' '.join(infiles)
        print 'Mosaicking grids...'
        os.system(mosaic)

def sieveMap(infile, out_file):
        '''sieve classification results to reinforce 2 acre MMU'''
        sieve='gdal_sieve.py -st 9 -8 ' + infile + ' ' + out_file
        print sieve
        os.system(sieve)

def remapRaster(infile, out_file, lookup):
        '''remap raster values to those in lookup table'''
        inmap = gdal.Open(infile)
        rows = inmap.RasterYSize
        cols = inmap.RasterXSize
        map_arr = inmap.ReadAsArray()

        #remap values
        remap_dict = df.getDictfromCSV(lookup,'\t',1,0)
        remap_dict[0]=2000 #ag
        remap_dict[255]=32767 #nodata
        map_out = map_arr.astype(np.int16)
        print 'input map labels', np.unique(map_out)
        for r in remap_dict:
                print 'reclassifying', r, ': ', remap_dict[r]
                outval=int(remap_dict[r])
                temp=np.equal(map_out, int(r))
                np.putmask(map_out, temp, int(remap_dict[r]))
                temp=None
        print 'output map labels', np.unique(map_out)
        #output raster
        driver=inmap.GetDriver()
        outDs = driver.Create(out_file, cols, rows, 1, GDT_Int16)
        outDs.SetGeoTransform(inmap.GetGeoTransform())
        outDs.SetProjection(inmap.GetProjection())
        outband = outDs.GetRasterBand(1)

        outband.WriteArray(map_out, 0 ,0)
        outband.SetNoDataValue(32767)
        outband.FlushCache()
          
'''______________________________________________________________________________'''

shapefiles=glob.glob(wd+'/**/poly*.shp')
grid_dir=wd+'segmented_grids'
segmented_map=map_file.split('.')[0]+'_segmented.tif'

##skip this -- testing shows insignificant difference in results
##rsmap = map_file.split('.')[0]+'_rs.tif'
##rsmap = resampleMap(map_file, 4)

if not os.path.exists(grid_dir): os.mkdir(grid_dir)

for shapefile in shapefiles:
        
        grid=shapefile.split('/')[-1][4:-4] #grid number from filename
        point=wd+'point'+grid+'.shp'
        allpt=wd+'all'+grid+'.shp'
        print 'Processing Grid, ', grid
        
        '''_____________________________________________________ZONAL STATISTICS OVER MAP'''
        if not point in glob.glob(wd+'/*'):
            print 'points for ', grid
            gme_isectpolyrst(shapefile, map_file)

            '''__________________________________________________OPEN VECTOR FILE AND ADD FIELD'''

            source = ogr.Open(shapefile, 1)
            layer = source.GetLayer()
            addIntField(layer, "label")
            field_names = getFieldNames(layer)
            maj_class = [x for x in field_names if 'levV' in x][:-2] #labels for majority (skip 254, sum)


            '''_______________________________________________________CALCULATE MAJORITY LABEL'''

            if 'label_1' in field_names:
                print '"label" field already exists for grid ' + str(grid) 
                sys.exit() 
            else:
                print 'CLASS LABELS FOR MAJORITY: ', ' '.join([x for x in maj_class])
                calcMajority(layer, maj_class) #calc majority, skipping first and last two columns (ID, 254, sum)
        
            '''_____________________________________________GET LABEL FOR UNLABELED POLYGONS'''
        
            gme_ptinpoly(shapefile, point, 'where=""label"=255"')
        
            if countFeatures(point) > 0:
                gme_ptinpoly(shapefile, allpt, 'where=""label"<255"')
                gme_isectptrst(point, map_file)
                nearestPoint(point, allpt)    
            
                subdict = update_dict(shapefile, point)
                subdict = nearestLabel(point, allpt, point.split('.')[0]+'nn.csv', subdict)

                updateLabel(layer, subdict)
                 
            '''_____________________________________________________________________RASTERIZE'''
        
            outrast = grid_dir + '/grid' + grid + '.tif'
            rasterizeSeg(shapefile, outrast, 'label')
'''________________________________________________________________________________MOSAIC'''
mosaicGrids(segmented_map, glob.glob(grid_dir+'/*.tif'))
'''_________________________________________________________________________________SIEVE'''
sieved_map=segmented_map.split('.')[0]+'_sieved.tif'
sieveMap(segmented_map, sieved_map)
'''_________________________________________________________CONVERT DATA TYPE AND RELABEL'''
remap=remapRaster(sieved_map, out_file,remap_file)
edit = "gdal_edit.py -a_srs EPSG:3071 " + out_file       #add projection info to metadata
os.system(edit) 
