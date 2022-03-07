import arcpy as ap
from arcpy import AddMessage as write
from arcpy import AddFieldDelimiters as fieldDelim
import os
import math
import datetime as dt
import pandas as pd
import numpy as np
import re
import sys


INNER_DIST = "35 Meters"
OUTER_DIST = "65 Meters"
CONVERSION_FACTOR = 0.3048  # Multiply to feet to get meters, divide from meters to get feet
HGT_M = 46
HGT_FT = 46 / CONVERSION_FACTOR

DVOF_FC = {
    'Aeronautic' : (100456, 100167, 100443, 100177, 100130),
    'Agriculture' : (100691, 100133, 100044),
    'Culture' : (100082, 100303, 100108, 100116),
    'HydroAidNavigation' : (100253),
    'Hydrography' : (100279, 100272, 100337, 100330, 100295),
    'Industry' : (100012, 100013, 100029, 100032, 100033, 100016, 100003, 100004, 100007, 100001),
    'Military' : (100128, 132626),
    'Recreation' : (100053, 100060, 121747, 100072, 154703, 100074),
    'Storage' : (100134, 100139, 100142, 100131),
    'Structure' : (100111, 100083, 100087, 100105, 100084, 100122),
    'TransportationGround' : (100161, 100164, 100190, 100147, 100186),
    'TransportationWater' : (100174, 132749),
    'UtilityInfrastructure' : (100028, 100200, 100558, 100025, 100045, 100018, 100021, 100179, 100202, 100028, 100687)
    }


def check_set_spatial_reference(dvof_data, tds_data):
    write("Checking Spatial Reference...")
    dvof_desc = ap.Describe(dvof_data)
    tds_desc = ap.Describe(tds_data)
    dvof_sr = dvof_desc.spatialReference
    tds_sr = tds_desc.spatialReference

    if dvof_sr.name != tds_sr.name:
        write("Applying {} projection to DVOF...".format(tds_sr.name))
        ap.DefineProjection_management(dvof_data, tds_sr)
    else:
        write("DVOF and TDS projection match confirmed: {}".format(tds_sr.name))

def create_near_table(dvof, fc_list, search_dist):
    # Create near tables for each subclass and corresponding DVOF subclass
    '''
        Run the GenerateNearTable tool with the created DVOF feature classes and the correspoding TDS feature class.
    '''
    write("Generating Near Table...")
    dist_int = int(re.search(r'\d+', search_dist).group())
    
    in_feats = dvof
    
    table_name = 'in_memory\\{0}{1}_Table'.format("test", dist_int)
    #write("\nCreating {}\n".format(table_name))
    location = 'NO_LOCATION'
    angle = 'NO_ANGLE'
    closest = 'ALL'
    closest_count = 0
    method = 'GEODESIC'
    ap.GenerateNearTable_analysis(in_feats, fc_list, table_name, search_dist,
                                    location, angle, closest, closest_count, method)
    
    cols = [field.name for field in ap.ListFields(table_name)]
    dataframe = pd.DataFrame(ap.da.FeatureClassToNumPyArray(in_table=table_name, field_names=cols, skip_nulls=False, null_value=-999999))
    
    dataframe.rename(columns={'IN_FID':'FID'}, inplace=True)
    #write("\n")
    #write(dataframe.columns)

    write("Near Table Generated.")

    return dataframe
   
def create_dvof_dataframe(in_dvof):
    write("Converting Source DVOF to DataFrame...\n")
    cols = ['FID', 'FEATURETYP', 'HEIGHTAGL', 'TYPENAME']
    #cols.remove('Shape')
    write(cols)
    dataframe = pd.DataFrame(ap.da.FeatureClassToNumPyArray(in_table=in_dvof, field_names=cols, skip_nulls=False, null_value=-999999))
    dataframe.HEIGHTAGL = dataframe.HEIGHTAGL.astype(float)
    write("Source DVOF:\n\n")
    write(dataframe)

    return dataframe

def fix_fc_names(row):
    row.NEAR_FC = os.path.basename(row.NEAR_FC)
    return row

def meters_to_feet(row):
    row.hgt = row.hgt / CONVERSION_FACTOR
    return row

def create_feature_dict(dataframe):
    write("Identifying TDS Feature Classes of interest...\n")
    tds_list = dataframe.NEAR_FC.unique().tolist()
    write("Feature Classes Near DVOF Points:")
    write(tds_list)
    feat_dict = {}
    for i in tds_list:
        feat_dict[i] = dataframe.NEAR_FID.loc[dataframe.NEAR_FC == i].tolist()
    
    return feat_dict

def create_feature_layers(feature_dict):
    write("Collecting TDS Features of interest...")
    lyr_list = []
    for key in feature_dict.keys():
        lyr_name = key + "_dvof"
        if len(feature_dict[key]) == 1:
            query = """{0} = {1}""".format(fieldDelim(key, 'objectid'), feature_dict[key][0])
        else:
            query = """{0} IN {1}""".format(fieldDelim(key, 'objectid'), tuple(feature_dict[key]))
        ap.management.MakeFeatureLayer(key, lyr_name, query)
        lyr_list.append(lyr_name)
    
    return lyr_list

def create_layer_dfs(feat_lyr_list):
    write("Converting TDS Features to DataFrame...")
    cols = ['objectid', 'fcsubtype', 'hgt']
    lyr_df_list = []
    for lyr in feat_lyr_list:
        new_df = pd.DataFrame(ap.da.FeatureClassToNumPyArray(in_table=lyr, field_names=cols, skip_nulls=False, null_value=-999999))
        new_df.rename(columns={'objectid':'NEAR_FID'}, inplace=True)
        lyr_df_list.append(new_df)
    
    return lyr_df_list


def combine_df_tables(table_list):
    write("Compiling TDS Feature DataFrames...")
    cols = table_list[0].columns.values
    dataframe = pd.DataFrame(columns=cols)

    for i in table_list:
        #write(i)
        dataframe = pd.concat([dataframe, i], ignore_index=True)
    
    dataframe[['NEAR_FID', 'fcsubtype']] = dataframe[['NEAR_FID', 'fcsubtype']].astype(int)


    #write(dataframe)
    return dataframe
    #write("\n\nThe Column names: {}".format(col_names))


def main(*argv):
    dvof_source = argv[0]
    tds_db = argv[1]

    ap.env.workspace = tds_db

    fc_list = ap.ListFeatureClasses("*Pnt")

    check_set_spatial_reference(dvof_source, tds_db)

    near_df = create_near_table(dvof_source, fc_list, OUTER_DIST)
    dvof_df = create_dvof_dataframe(dvof_source)

    merged_df = pd.merge(near_df, dvof_df,on='FID', how='left')

    merged_df =  merged_df.apply(fix_fc_names, axis='columns')

    feature_dict = create_feature_dict(merged_df)

    tds_layer_list = create_feature_layers(feature_dict)
    tds_lyr_df_list = create_layer_dfs(tds_layer_list)
    
    merged_lyrs = combine_df_tables(tds_lyr_df_list)

    complete_df = pd.merge(merged_df, merged_lyrs,on='NEAR_FID', how='left')

    complete_df = complete_df[complete_df.hgt != -999999.0000]

    complete_df = complete_df.apply(meters_to_feet, axis='columns')

    complete_df['hgt_delta'] = complete_df.apply(lambda h: round(abs(h.HEIGHTAGL - h.hgt),3), axis='columns')

    # Identify DVOF points that match spatially and height
    matched_points = complete_df[(complete_df.NEAR_DIST < 35) & (complete_df.hgt_delta < 5)].FID.tolist()

    complete_df = complete_df[~complete_df.FID.isin(matched_points)]

    write(complete_df)

    write("\nDVOF Points matched:\n")
    write(matched_points)



if __name__=='__main__':
    ap.env.overwriteOutput = True
    argv = tuple(ap.GetParameterAsText(i) for i in range(ap.GetArgumentCount()))
    now = dt.datetime.now()
    main(*argv)
    write(dt.datetime.now() - now)