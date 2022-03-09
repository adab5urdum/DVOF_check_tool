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
HGT_FT = round(46 / CONVERSION_FACTOR, 4)

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


def create_output_dvof(folder_path, fc_template):
    write("Setting up output files...")
    out_name = "DVOF_Output"
    out_db = ap.management.CreateFileGDB(folder_path, out_name)
    
    dvof_fc_name = 'DVOF_Points'
    geometry_type = "POINT"
    dvof_out_path = ap.management.CreateFeatureclass(out_db, dvof_fc_name, geometry_type, fc_template)
    
    return dvof_out_path

def create_output_checks(out_path):
    '''
        Create shapefile to mark points that need to be manually inspected.
    '''
    out_name = "DVOF_Checks.shp"
    out_file = os.path.join(out_path, out_name)

    ap.CreateFeatureclass_management(out_path, out_name, 'POINT')

    ap.management.AddField(out_file, "typecode", "TEXT")
    ap.management.AddField(out_file, "typename", "TEXT")
    ap.management.AddField(out_file, "check", "TEXT")

    return out_file

def move_to_checkfile(check_db, row):
    in_fields = ['typecode', 'typename', 'check', 'SHAPE@']
    
    with ap.da.InsertCursor(check_db, in_fields) as iCur:
        iCur.insertRow(row)

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

def get_count(fc_layer):
    results = int(ap.GetCount_management(fc_layer).getOutput(0))

    return results

def check_tds_dvof_presence():
    feature_layer_dict = {"Points": [], "Surfaces": []}
    for fc in DVOF_FC.keys():
        srf_fc = fc + "Srf"
        pnt_fc = fc + "Pnt"
        srf_name = srf_fc + "_lyr"
        pnt_name = pnt_fc + "_lyr"
        if type(DVOF_FC[fc]) is int:
            where_clause = """{0} >= {1} AND {2} = {3}""".format(fieldDelim(pnt_fc, 'hgt'), HGT_M, fieldDelim(pnt_fc, 'fcsubtype'), DVOF_FC[fc])
        else:
            where_clause = """{0} >= {1} AND {2} IN {3}""".format(fieldDelim(pnt_fc, 'hgt'), HGT_M, fieldDelim(pnt_fc, 'fcsubtype'), DVOF_FC[fc])
        
        # Create Surface Layer and check count
        ap.management.MakeFeatureLayer(srf_fc, srf_name, where_clause)
        srf_count = get_count(srf_name)
        if srf_count == 0:
            write("No {} features above 46 meters present. Skipping...".format(srf_fc))
        else:
            write("{0} {1} features above 46 meters found.".format(srf_count, srf_fc))
            feature_layer_dict['Surfaces'].append(srf_name)
        
        # Create Point Layer and check count
        ap.management.MakeFeatureLayer(pnt_fc, pnt_name, where_clause)
        pnt_count = get_count(pnt_name)
        if pnt_count == 0:
            write("No {} features above 46 meters present. Skipping...".format(pnt_fc))
        else:
            write("{0} {1} features above 46 meters found.".format(pnt_count, pnt_fc))
            feature_layer_dict['Points'].append(pnt_name)

    return feature_layer_dict

def prepare_tds_data(layer_dict):
    working_list = layer_dict['Points']
    if len(layer_dict['Surfaces']) == 0:
        write("No qualifying surface features present.")
        return working_list
    else:
        write("Surface features with qualifying height present. Preparing for processing...")
        for lyr in layer_dict['Surfaces']:
            name = lyr+"_pnts"
            fc_path = "in_memory\\{}".format(name)
            pnt_type = "CENTROID"
            ap.management.FeatureToPoint(lyr, fc_path, pnt_type)
            working_list.append(fc_path)
            write("{} converted to Points.".format(lyr))
    return working_list

def check_source_dvof(dvof_data):
    total_count = get_count(dvof_data)
    lyr_name = "working_dvof_lyr"
    if get_field_type(dvof_data, 'HEIGHTAGL') == 'String':
        working_fids = find_dvof_by_hgt(dvof_data)
        if len(working_fids) == 1:
            where_clause = """{0} = {1}""".format(fieldDelim(dvof_data, 'FID'), working_fids[0])
        else:
            where_clause = """{0} IN {1}""".format(fieldDelim(dvof_data, 'FID'), tuple(working_fids))
        ap.management.MakeFeatureLayer(dvof_data, lyr_name, where_clause)
    else:
        where_clause = """{0} >= {1}""".format(fieldDelim(dvof_data, 'HEIGHTAGL'), HGT_FT)
        ap.management.MakeFeatureLayer(dvof_data, lyr_name, where_clause)
    work_count = get_count(lyr_name)
    write("{0} of {1} DVOF Points meet height threshold.\nChecking qualifying points against TDS.".format(work_count, total_count))

    return lyr_name

def get_field_type(data, field_name):
    fields = ap.ListFields(data)
    for i in fields:
        if i.name == field_name:
            return i.type
    else:
        return None

def find_dvof_by_hgt(dvof_data):
    fid_list = []
    fields = ['FID', 'HEIGHTAGL']
    with ap.da.SearchCursor(dvof_data, fields) as cursor:
        for row in cursor:
            if float(row[1]) >= HGT_FT:
                fid_list.append(row[0])
    
    return fid_list

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

    # For the case where there is only one working tds layer, return the layer df instead of a list
    if len(lyr_df_list) == 1:
        return lyr_df_list[0]
    
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


def identify_repeats(in_list):
    repeats_list = []
    for i in in_list:
        if in_list.count(i) > 1:
            if i not in repeats_list:
                repeats_list.append(i)
    
    return repeats_list

def add_checks(check_file, dvof_source, fid_list):
    fields = ['FID','FEATURETYP', 'TYPENAME', 'SHAPE@']
    with ap.da.SearchCursor(dvof_source, fields) as cursor:
        for row in cursor:
            if row[0] in fid_list:
                check_row = [row[1], row[2], 'More than one TDS feature match.', row[-1]]
                move_to_checkfile(check_file, check_row)


def main(*argv):
    # Bring in Parameters from ArcMap:
    #   1. DVOF source shapefile/feature class
    #   2. Delivery TDS feature dataset
    #   3. Output Folder
    dvof_source = argv[0]
    tds_db = argv[1]
    out_path = argv[2]
    fc_template = argv[3]

    # Set the workspace to the TDS feature dataset
    ap.env.workspace = tds_db

    # Create output files at designated output folder
    out_dvof = create_output_dvof(out_path, fc_template)
    write(out_dvof)

    out_checks = create_output_checks(out_path)
    write(out_checks)

    # Ensure that the DVOF input and TDS input have the same spatial reference
    # This is important to perform Geodesic measurements
    check_set_spatial_reference(dvof_source, tds_db)

    # Check TDS feature classes for qualifying features (hgt >= 46)
    tds_fc_dict = check_tds_dvof_presence()

    # Prepare working list of tds features for DVOF comparison
    working_tds_list = prepare_tds_data(tds_fc_dict)
    write(working_tds_list)

    # Check Source DVOF for points above 46 m threshold. Anything below is ignored.
    working_dvof_lyr = check_source_dvof(dvof_source)

    # Create the Near Tabel.
    # This is the primary table for this tool and will be added to and manipulated as needed.
    near_df = create_near_table(working_dvof_lyr, working_tds_list, OUTER_DIST)
    
    # Convert the necessary data from DVOF to Pandas DataFrame
    dvof_df = create_dvof_dataframe(working_dvof_lyr)

    # Merge the DVOF source data on to the Near Table based on DVOF IDs
    working_df = pd.merge(near_df, dvof_df,on='FID', how='left')

    # If there is only 1 TDS layer then the processing for multiple Near FCs is not necessary
    if len(working_tds_list) > 1:

        # Isolate the feature class names from their paths 
        working_df =  working_df.apply(fix_fc_names, axis='columns')

        # Create a dictionary of the OIDs found in each feature class
        feature_dict = create_feature_dict(working_df)

        # Using the identified OIDs create layers from each feature class
        tds_layer_list = create_feature_layers(feature_dict)

        # Pull the necessary information from the created layers into Pandas DataFrame
        tds_lyr_df_list = create_layer_dfs(tds_layer_list)

        # Combine the layer DataFrames into a single DataFrame
        combined_lyrs = combine_df_tables(tds_lyr_df_list)

        # Merge the layer data onto the Near table with the TDS OIDs as keys
        working_df = pd.merge(working_df, combined_lyrs,on='NEAR_FID', how='left')
    else:
        # Convert layer to Pandas DataFrame
        single_lyr_df = create_layer_dfs(working_tds_list)

        # Merge layer DataFrame with Working DataFrame
        working_df = pd.merge(working_df, single_lyr_df, on='NEAR_FID', how='left')


    # Drop all rows in which the TDS height is Null (-999999) --- Redundant step
    working_df = working_df[working_df.hgt != -999999.0000]

    # Convert the TDS heights from meters to feet for comparison
    working_df = working_df.apply(meters_to_feet, axis='columns')

    # Create a column to record the absolute difference between the DVOF height and TDS height
    working_df['hgt_delta'] = working_df.apply(lambda h: round(abs(h.HEIGHTAGL - h.hgt),3), axis='columns')

    # Identify DVOF points that match spatially (within 1 meter) and height (difference less than 5)
    matched_points = working_df[(working_df.NEAR_DIST < 1) & (working_df.hgt_delta < 5)].FID.tolist()
    if 'NEAR_FC' in working_df.columns:
        matched_tds_points = {fc: working_df[('NEAR_FC' == fc) & (working_df.FID.isin(matched_points))].NEAR_FID.tolist() for fc in working_df[working_df.FID.isin(matched_points)].NEAR_FC.tolist()}
    else:
        matched_tds_points = {working_tds_list[0]: working_df[(working_df.NEAR_DIST < 1) & (working_df.hgt_delta < 5)].NEAR_FID.tolist()}

    # Collect any DVOF points that match with more than one TDS feature
    check_matched = identify_repeats(matched_points)
    if len(check_matched) == 0:
        write("Matched DVOF Points Verified.")
    else:
        write("{} points identified for inspection.".format(len(check_matched)))
        # Remove the ids to be checked from the "No Change" list
        matched_points = [x for x in matched_points if x not in check_matched]
        # Add check_matched points to DVOF checkfile
        add_checks(out_checks, dvof_source, check_matched)

    # Remove the rows containing the matched DVOF points
    working_df = working_df[~working_df.FID.isin(matched_points)]

    # Identify DVOF points which are close enough for no change to be recorded
    no_change_points = working_df[(working_df.NEAR_DIST.le(35)) & (working_df.hgt_delta < 5)].FID.tolist()
    if 'NEAR_FC' in working_df.columns:
        no_change_tds_points = {fc: working_df["NEAR_FC" == fc & (working_df.FID.isin(no_change_points))].NEAR_FID.tolist() for fc in working_df[working_df.FID.isin(no_change_points)].NEAR_FC.tolist()}
    else:
        no_change_tds_points = {working_tds_list[0]: working_df[(working_df.NEAR_DIST.le(35)) & (working_df.hgt_delta < 5)].NEAR_FID.tolist()}

    # Collect any no_change_points that fit with more than one TDS point
    check_no_change = identify_repeats(no_change_points)
    if len(check_no_change) == 0:
        write("No Changed DVOF Points Verified.")
    else:
        write("{} points identified for inspection.".format(len(check_no_change)))
        # Remove the ids to be checked from the "No Change" list
        no_change_points = [x for x in no_change_points if x not in check_no_change]
        # Add check_no_change points to DVOF checkfile
        add_checks(out_checks, dvof_source, check_no_change)
    
    # Remove the rows containing the matched DVOF points
    working_df = working_df[~working_df.FID.isin(no_change_points)]

    remainder_dvof_points = working_df.FID.tolist()
    check_remainder = identify_repeats(remainder_dvof_points)
    
    if len(check_remainder) == 0:
        write("Remainder Points Verified.")
    else:
        write("{} points identified for inspection.".format(len(check_remainder)))
        remainder_dvof_points = [x for x in remainder_dvof_points if x not in check_remainder]
        # add check_remainder points to DVOF checkfile
        add_checks(out_checks, dvof_source, check_remainder)
    
    checks_remain_df = working_df[~working_df.FID.isin(remainder_dvof_points)]
    if 'NEAR_FC' in checks_remain_df.columns:
        checked_tds_points = {fc: checks_remain_df.NEAR_FID.tolist() for fc in checks_remain_df.NEAR_FC.tolist()}
    else:
        checked_tds_points = {working_tds_list[0]: checks_remain_df.NEAR_FID.tolist()}
    
    if 'Near_FC' in working_df.columns:
        change_points = {fc: working_df['NEAR_FC' == fc & (working_df.FID.isin(remainder_dvof_points))].NEAR_FID.tolist() for fc in working_df[(working_df.FID.isin(remainder_dvof_points))].NEAR_FC.tolist()}
    else:
        change_points = {working_tds_list[0]: working_df[working_df.FID.isin(remainder_dvof_points)].NEAR_FID.tolist()}


    # Populate change points in output dvof file

    # Look for adds in tds features (accounting for tds features already checked or added with change) and populate in output dvof file

    # Look for deletes in dvof points that did not appear in near table (thus is not close enough to a tds point) and populate in output dvof file
    
    write("\nBefore final filter:")
    write(working_df)
    write(remainder_dvof_points)
    write(change_points)

    write("\nThis is what remains:")
    write(checks_remain_df)
    write(checked_tds_points)

    write("\nDVOF Points matched:")
    write(matched_points)
    multiples_list = identify_repeats(matched_points)
    if len(multiples_list) > 0:
        write("\nMultiples Found:\n")
        write(multiples_list)


if __name__=='__main__':
    ap.env.overwriteOutput = True
    argv = tuple(ap.GetParameterAsText(i) for i in range(ap.GetArgumentCount()))
    now = dt.datetime.now()
    main(*argv)
    write(dt.datetime.now() - now)