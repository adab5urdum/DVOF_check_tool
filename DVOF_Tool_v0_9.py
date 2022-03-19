'''
    Name: DVOF Reporting Tool
    Description: This tool is designed to identify the Adds, Edits, and Deletes of DVOF for a TDS TPC.
                 Given a Source DVOF (either shapefile or feature class) and a TDS (feature dataset),
                 this tool will identify features within the TDS data that qualify for DVOF and compare
                 them to the existing DVOF points.
                 
                 If there is a change in position greater than 35 Meters or a change in height greater than 5 feet, 
                 then the feature will be marked as a 'Change'.
                 
                 If a TDS feature is not within 65 Meters of a DVOF point, then it will makred as an 'Add'.
                 
                 If a DVOF point does not have any qualifying TDS features within 65 Meters, then it will be
                 marked as a 'Delete'.

                 Any DVOF points that have multiple TDS features which qualify for a match will be flagged in a check file.
    
    Created by: John Jackson, Nat Cagle, and Kristen Hall
'''

import arcpy as ap
from arcpy import AddMessage as write
from arcpy import AddFieldDelimiters as fieldDelim
import os
import datetime as dt
import pandas as pd
import re
import sys

# Defining global constants
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

#-------------------------------------------------------------------------------------------------------------------------------------------------

def create_output_dvof(folder_path, fc_template):
    write("Setting up output files...")
    out_name = "DVOF_Output"
    out_db = ap.management.CreateFileGDB(folder_path, out_name)
    
    dvof_fc_name = 'DVOF_Points'
    geometry_type = "POINT"
    has_m = "DISABLED"
    has_z = "DISABLED"
    spatial_reference = ap.Describe(fc_template).spatialReference
    dvof_out_path = ap.management.CreateFeatureclass(out_db, dvof_fc_name, geometry_type, fc_template, has_m, has_z, spatial_reference)
    
    return dvof_out_path
#-------------------------------------------------------------------------------------------------------------------------------------------------

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
#-------------------------------------------------------------------------------------------------------------------------------------------------

def move_to_checkfile(check_db, row):
    in_fields = ['typecode', 'typename', 'check', 'SHAPE@']
    
    with ap.da.InsertCursor(check_db, in_fields) as iCur:
        iCur.insertRow(row)
#-------------------------------------------------------------------------------------------------------------------------------------------------

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
#-------------------------------------------------------------------------------------------------------------------------------------------------

def get_count(fc_layer):
    results = int(ap.GetCount_management(fc_layer).getOutput(0))

    return results
#-------------------------------------------------------------------------------------------------------------------------------------------------

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
#-------------------------------------------------------------------------------------------------------------------------------------------------

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
#-------------------------------------------------------------------------------------------------------------------------------------------------

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
#-------------------------------------------------------------------------------------------------------------------------------------------------

def get_field_type(data, field_name):
    fields = ap.ListFields(data)
    for i in fields:
        if i.name == field_name:
            return i.type
    else:
        return None
#-------------------------------------------------------------------------------------------------------------------------------------------------

def find_dvof_by_hgt(dvof_data):
    fid_list = []
    fields = ['FID', 'HEIGHTAGL']
    with ap.da.SearchCursor(dvof_data, fields) as cursor:
        for row in cursor:
            if float(row[1]) >= HGT_FT:
                fid_list.append(row[0])
    
    return fid_list
#-------------------------------------------------------------------------------------------------------------------------------------------------

def create_near_table(dvof, fc_list, search_dist):
    # Create near tables for each subclass and corresponding DVOF subclass
    '''
        Run the GenerateNearTable tool with the created DVOF feature classes and the correspoding TDS feature class.
    '''
    write("Generating Near Table...")
    dist_int = int(re.search(r'\d+', search_dist).group())
    
    in_feats = dvof
    
    table_name = 'in_memory\\{0}{1}_Table'.format("test", dist_int)
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

    write("Near Table Generated.")

    return dataframe
#-------------------------------------------------------------------------------------------------------------------------------------------------
#    
def create_dvof_dataframe(in_dvof):
    write("Converting Source DVOF to DataFrame...\n")
    cols = ['FID', 'FEATURETYP', 'HEIGHTAGL', 'TYPENAME']
    write(cols)
    dataframe = pd.DataFrame(ap.da.FeatureClassToNumPyArray(in_table=in_dvof, field_names=cols, skip_nulls=False, null_value=-999999))
    dataframe.HEIGHTAGL = dataframe.HEIGHTAGL.astype(float)
    write("Source DVOF:\n\n")
    write(dataframe)

    return dataframe
#-------------------------------------------------------------------------------------------------------------------------------------------------

def fix_fc_names(row):
    row.NEAR_FC = os.path.basename(row.NEAR_FC)
    return row
#-------------------------------------------------------------------------------------------------------------------------------------------------

def meters_to_feet(row):
    row.hgt = row.hgt / CONVERSION_FACTOR
    return row
#-------------------------------------------------------------------------------------------------------------------------------------------------

def check_table_for_tds_layer(tds_list, dataframe):
    near_fc_list = dataframe.NEAR_FC.unique().tolist()
    tds_layers_match = all(fc in tds_list for fc in near_fc_list)
    if tds_layers_match:
        write("\nNear table TDS references accepted!")
    else: 
        write("\nSomething has gone wrong....\nCall John\n")
        write("Near Table Entries: {0} \nWorking Layers: {1}".format(near_fc_list, tds_list))
    
    return tds_layers_match
#-------------------------------------------------------------------------------------------------------------------------------------------------

def create_layer_dfs(feat_lyr_list):
    write("\nConverting TDS Features to DataFrame...")
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

#-------------------------------------------------------------------------------------------------------------------------------------------------

def identify_checks(oids_list, checks_file, source_dvof, step_name):
    checks_list = identify_repeats(oids_list)
    if len(checks_list) == 0:
        write("\n{} DVOF Points Verified.".format(step_name))
    else:
        write("\n{} points identified for inspection.".format(len(checks_list)))
        # Remove the ids to be checked from the "No Change" list
        oids_list = [x for x in oids_list if x not in checks_list]
        # Add check_matched points to DVOF checkfile
        add_checks(checks_file, source_dvof, checks_list)
    
    return oids_list
#-------------------------------------------------------------------------------------------------------------------------------------------------

def combine_df_tables(table_list):
    write("Compiling TDS Feature DataFrames...")
    cols = table_list[0].columns.values
    dataframe = pd.DataFrame(columns=cols)

    for i in table_list:
        dataframe = pd.concat([dataframe, i], ignore_index=True)
    
    dataframe[['NEAR_FID', 'fcsubtype']] = dataframe[['NEAR_FID', 'fcsubtype']].astype(int)


    return dataframe

#-------------------------------------------------------------------------------------------------------------------------------------------------

def identify_repeats(in_list):
    repeats_list = []
    for i in in_list:
        if in_list.count(i) > 1:
            if i not in repeats_list:
                repeats_list.append(i)
    
    return repeats_list
#-------------------------------------------------------------------------------------------------------------------------------------------------

def add_checks(check_file, dvof_source, fid_list):
    fields = ['FID','FEATURETYP', 'TYPENAME', 'SHAPE@']
    with ap.da.SearchCursor(dvof_source, fields) as cursor:
        for row in cursor:
            if row[0] in fid_list:
                check_row = [row[1], row[2], 'More than one TDS feature match.', row[-1]]
                move_to_checkfile(check_file, check_row)
#-------------------------------------------------------------------------------------------------------------------------------------------------

def ignore_tds(tracker, list, dataframe):
    for key in tracker.keys():
        tds_list = dataframe[(dataframe.FID.isin(list)) & (dataframe.NEAR_FC == key)].NEAR_FID.tolist()
        for oid in tds_list:
            tracker[key].append(oid)
    
    return tracker
#-------------------------------------------------------------------------------------------------------------------------------------------------

def compile_change_points(dataframe):
    fc_list = dataframe.NEAR_FC.unique().tolist()
    change_dict = {fc: [] for fc in fc_list}
    for key in change_dict.keys():
        change_dict[key] = dataframe.NEAR_FID.tolist()
    
    return change_dict
#-------------------------------------------------------------------------------------------------------------------------------------------------

def populate_adds(tds_list, ignore_dict, dvof_fc):
    add_list = []
    tds_fields = ["OID@", "fcsubtype", "hgt", "SHAPE@"]
    otherremar = 'A'
    voidentifi = ''

    for lyr in tds_list:
        if len(ignore_dict[lyr]) == 1:
            where_clause = """{0} <> {1}""".format(fieldDelim(lyr, 'objectid'), ignore_dict[lyr][0])
        else:
            where_clause = """{0} NOT IN {1}""".format(fieldDelim(lyr, 'objectid'), tuple(ignore_dict[lyr]))
        with ap.da.SearchCursor(lyr, tds_fields, where_clause) as cursor:
            for row in cursor:
                geometry = row[-1]
                for pnt in geometry:
                    long_x = pnt.X
                    lat_y = pnt.Y
                featuretyp = '' # ***Conversion function in progress: update later***
                hgt_feet = round(row[2] / CONVERSION_FACTOR, 3)
                insert_row = [voidentifi, lat_y, long_x, featuretyp, hgt_feet, otherremar, geometry]
                add_list.append(row[0])
                move_to_delivery_dvof(dvof_fc, insert_row)
    write("\n{} Add Points populated:".format(len(add_list)))
#-------------------------------------------------------------------------------------------------------------------------------------------------

def populate_edits(edits_dict, dvof_fc):
    add_list = []
    tds_fields = ["OID@", "fcsubtype", "hgt", "SHAPE@"]
    otherremar = 'C'
    voidentifi = ''

    for lyr in edits_dict:
        if len(edits_dict[lyr]) == 1:
            where_clause = """{0} = {1}""".format(fieldDelim(lyr, 'objectid'), edits_dict[lyr][0])
        else:
            where_clause = """{0} IN {1}""".format(fieldDelim(lyr, 'objectid'), tuple(edits_dict[lyr]))
        with ap.da.SearchCursor(lyr, tds_fields, where_clause) as cursor:
            for row in cursor:
                geometry = row[-1]
                for pnt in geometry:
                    long_x = pnt.X
                    lat_y = pnt.Y
                hgt_feet = round(row[2] / CONVERSION_FACTOR, 3)
                featuretyp = 0 # ***Conversion function in progress: update later***
                insert_row = [voidentifi, lat_y, long_x, featuretyp, hgt_feet, otherremar, geometry]
                add_list.append(row[0])
                move_to_delivery_dvof(dvof_fc, insert_row)
    write("\n{} Change Points populated:".format(len(add_list)))

#-------------------------------------------------------------------------------------------------------------------------------------------------

def move_to_delivery_dvof(delivery_fc, row):
    delivery_fields = ['VOIDENTIFI', 'LATITUDE', 'LONGITUDE', 'FEATURETYP', 'HEIGHTAGL', 'OTHERREMAR', 'SHAPE@']
    with ap.da.InsertCursor(delivery_fc, delivery_fields) as iCur:
        iCur.insertRow(row)
#-------------------------------------------------------------------------------------------------------------------------------------------------

def fill_in_defaults_dvof(wac_file, country_bndry_file, imagery_footprint, delivery_fc):
    dvof_fields = ['COUNTRYCD', 'WAC', 'SOURCEDT', 'HEIGHTAMSL',
                    'SINGLEMULT', 'MULTIPLENO', 'VALIDATION', 'HORIZACC',
                    'HGTAGLACC', 'HGTAMSLACC', 'SECURITYCL', 'SURFACEMAT',
                    'HORIZDATUM', 'PRODUCER', 'PROCESSCD', 'LIGHTINGCD', 'SHAPE@']
    wac_fields = ['wac_index', 'SHAPE@']
    country_fields = ['CNTRY_NAME', 'SHAPE@']
    imagery_fields= ['IngestDate', 'SHAPE@']

    # Default values 
    tds_null = -999999
    sourcedt = ''               # 2
    heightamsl = tds_null       # 3
    singlemult = 'S'            # 4
    multipleno = 1              # 5
    validation = 3              # 6
    horizacc = tds_null         # 7
    hgtaglacc = tds_null        # 8
    hgtamslacc = tds_null       # 9
    securitycl = 'Unclassified' # 10
    surfacemat = 'Unknown'      # 11
    horizdatum = 'WGE'          # 12
    producer = 'E'              # 13
    processcd = 'OT'            # 14
    lightingcd = 'U'            # 15

    with ap.da.UpdateCursor(delivery_fc, dvof_fields) as uCur:
        for row in uCur:
            with ap.da.SearchCursor(country_bndry_file, country_fields) as c_cur:
                for c_row in c_cur:
                    if not row[-1].disjoint(c_row[-1]):
                        countrycd = c_row[0]

            with ap.da.SearchCursor(wac_file, wac_fields) as w_cur:
                for w_row in w_cur:
                    if not row[-1].disjoint(w_row[-1]):
                        wac_value = w_row[0]
            
            with ap.da.SearchCursor(imagery_footprint, imagery_fields) as f_cur:
                for f_row in f_cur:
                    if not row[-1].disjoint(f_row[-1]):
                        if sourcedt == '':
                            sourcedt = f_row[0]
                        elif sourcedt < f_row[0]:
                            sourcedt = f_row[0]
            row[0] = countrycd
            row[1] = wac_value
            row[2] = sourcedt
            row[3] = heightamsl
            row[4] = singlemult
            row[5] = multipleno
            row[6] = validation
            row[7] = horizacc
            row[8] = hgtaglacc
            row[9] = hgtamslacc
            row[10] = securitycl
            row[11] = surfacemat
            row[12] = horizdatum
            row[13] = producer
            row[14] = processcd
            row[15] = lightingcd

            uCur.updateRow(row)
    write("\nCoutnry Code, WAC, Source Date, and Defaults Populated.")

#-------------------------------------------------------------------------------------------------------------------------------------------------

def populate_deletes(deletes_list, source_dvof, delivery_dvof):
    if len(deletes_list) == 1:
        where_clause = """{0} = {1}""".format(fieldDelim(source_dvof, 'FID'), deletes_list[0])
    else:
        where_clause = """{0} IN {1}""".format(fieldDelim(source_dvof, 'FID'), tuple(deletes_list))
    
    source_fields = ['VOIDENTIFI', 'LATITUDE', 'LONGITUDE', 'FEATURETYP',
                'HEIGHTAGL', 'COUNTRYCD', 'WAC', 'SOURCEDT', 'HEIGHTAMSL',
                'SINGLEMULT', 'MULTIPLENO', 'VALIDATION', 'HORIZACC',
                'HGTAGLACC', 'HGTAMSLACC', 'SECURITYCL', 'SURFACEMAT', 'HORIZDATUM', 
                'PROCESSCD', 'LIGHTINGCD', 'SHAPE@']

    out_fields = ['VOIDENTIFI', 'LATITUDE', 'LONGITUDE', 'FEATURETYP',
                'HEIGHTAGL', 'COUNTRYCD', 'WAC', 'SOURCEDT', 'HEIGHTAMSL',
                'SINGLEMULT', 'MULTIPLENO', 'VALIDATION', 'HORIZACC',
                'HGTAGLACC', 'HGTAMSLACC', 'SECURITYCL', 'SURFACEMAT', 'HORIZDATUM', 
                'PROCESSCD', 'LIGHTINGCD', 'SHAPE@','OTHERREMAR']

    with ap.da.SearchCursor(source_dvof, source_fields, where_clause) as del_cur:
        for row in del_cur:
            with ap.da.InsertCursor(delivery_dvof, out_fields) as in_cur:
                new_row = [x for x in row]
                new_row.append('D')
                in_cur.insertRow(new_row)
    
   #with ap.da.UpdateCursor(delivery_dvof, ['OTHERREMAR']) as update:
   #    for update_row in update:
   #        if not update_row[0]:
   #            update_row[0] = 'D'
   #            update.updateRow(update_row)
    
    write("\n{} Delete points populated.".format(len(deletes_list)))

#-------------------------------------------------------------------------------------------------------------------------------------------------

def main(*argv):
    '''
        Set-up: To facilitate the processes in this tool the following steps are necessary:
                    1. Assign inputs to variables and define workspace
                    2. Create output files: Checks Shapefile, DVOF Feature Class
                    3. Ensure that the TDS and DVOF spatial references match
    '''
    # Bring in Parameters from ArcMap:
    #   1. DVOF source shapefile/feature class
    #   2. Delivery TDS feature dataset
    #   3. Output Folder
    dvof_source = argv[0]
    tds_db = argv[1]
    out_path = argv[2]
    fc_template = argv[3]
    wac = argv[4]
    country_boundaries = argv[5]
    imagery_footprint = argv[6]

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

    '''
        Pre-Processing: The data is imported and conditioned. First, layers are created
                        to filter out any unnecessary features. Next, the DVOF and TDS data
                        is used to create a Near Table, this provides the spatial distance 
                        information necessary. Finally, additional information from both the
                        DVOF and TDS data is pulled, conditoned, and added to the created table.
    '''

    # Check TDS feature classes for qualifying features (hgt >= 46)
    tds_fc_dict = check_tds_dvof_presence()

    # Prepare working list of tds features for DVOF comparison
    working_tds_list = prepare_tds_data(tds_fc_dict)
    write(working_tds_list)

    # Create means of tracking the features that need to be edited or added
    ignore_tracker = {lyr:[] for lyr in working_tds_list}
    #add_edit_tracker = {lyr:{} for lyr in working_tds_list}

    # Check Source DVOF for points above 46 m threshold. Anything below is ignored.
    working_dvof_lyr = check_source_dvof(dvof_source)

    # Create the Near Tabel.
    # This is the primary table for this tool and will be added to and manipulated as needed.
    near_df = create_near_table(working_dvof_lyr, working_tds_list, OUTER_DIST)
    
    # Convert the necessary data from DVOF to Pandas DataFrame
    dvof_df = create_dvof_dataframe(working_dvof_lyr)

    # Identify DVOF FIDs not in the Near table and flag as Deletes
    dvof_fids = dvof_df.FID.unique().tolist()
    near_dvof_fids = near_df.FID.unique()
    delete_list = [x for x in dvof_fids if x not in near_dvof_fids]
    write("\n\n{0} DVOF Points of {1} total marked as Delete".format(len(delete_list), len(dvof_fids)))


    # Merge the DVOF source data on to the Near Table based on DVOF IDs
    working_df = pd.merge(near_df, dvof_df,on='FID', how='left')

    # If there is only 1 TDS layer then the processing for multiple Near FCs is not necessary
    if "NEAR_FC" in working_df.columns:

        # Isolate the feature class names from their paths 
        working_df =  working_df.apply(fix_fc_names, axis='columns')

        # Ensure that the layers in the working_tds_list are matched to the Near_FC values
        table_check_approved = check_table_for_tds_layer(working_tds_list, working_df)
        if not table_check_approved:
            # Kill program
            sys.exit(0)

        # Pull the necessary information from the created layers into Pandas DataFrame
        tds_df_list = create_layer_dfs(working_tds_list)

        # Combine the layer DataFrames into a single DataFrame
        combined_lyrs = combine_df_tables(tds_df_list)

        # Merge the layer data onto the Near table with the TDS OIDs as keys
        working_df = pd.merge(working_df, combined_lyrs,on='NEAR_FID', how='left')
    else:
        # Convert layer to Pandas DataFrame
        single_lyr_df = create_layer_dfs(working_tds_list)

        # Merge layer DataFrame with Working DataFrame
        working_df = pd.merge(working_df, single_lyr_df, on='NEAR_FID', how='left')
        
        # Add NEAR_FC Column and populate with the single fc, simplifies code in the 'Processing' phase
        working_df["NEAR_FC"] = working_tds_list[0]
    
    # Convert the TDS heights from meters to feet for comparison
    working_df = working_df.apply(meters_to_feet, axis='columns')

    # Create a column to record the absolute difference between the DVOF height and TDS height
    working_df['hgt_delta'] = working_df.apply(lambda h: round(abs(h.HEIGHTAGL - h.hgt),3), axis='columns')

    '''
        Processing: Filter through the pre-processed table to identify any DVOF points that can be safely
                    ignored. Add any that are uncertain to the Check File. Whatever remains is populated
                    as a Change ('C') in post-processing.
    '''
    # Step 1 -- Identify any DVOF points that are 'near-perfect' matches to a TDS feature.
    #           Flag any that are close match for more than one TDS feature.
    step_1 = "Matched"

    # Select DVOF FIDs that are within 1 meter of a TDS feature with less than a 5 foot difference in height.
    matched_points = working_df[(working_df.NEAR_DIST < 1) & (working_df.hgt_delta < 5)].FID.tolist()

    # Since these are either matches or need to be checked, ignore the respective TDS features from now on.
    ignore_tracker = ignore_tds(ignore_tracker, matched_points, working_df)

    # The DVOF FIDs are scanned for repeats, meaning more than one TDS feature applied to the same DVOF point.
    # These are added to the Checks file.
    matched_points = identify_checks(matched_points, out_checks, dvof_source, step_1)
   
    # Remove the rows containing the matched DVOF points.
    working_df = working_df[~working_df.FID.isin(matched_points)]

    # Step 2 -- Identify any DVOF points that are within tolerances to not qualify for a Change entry.
    #           The individual steps are similar as those in Step 1, see explanations above.
    #           Difference between Step 1 and Step 2: the distance tolerance is expanded from 1 to 35 meters.
    step_2 = "No-Change"
    no_change_points = working_df[(working_df.NEAR_DIST.le(35)) & (working_df.hgt_delta < 5)].FID.tolist()
    ignore_tracker = ignore_tds(ignore_tracker, no_change_points, working_df)
    no_change_points = identify_checks(no_change_points, out_checks, dvof_source, step_2)
   
    # Remove the rows containing the 'no change' DVOF points.
    working_df = working_df[~working_df.FID.isin(no_change_points)]

    # Step 3 -- Filter any DVOF points that need to be checked and collect the rest as Change Points.
    step_3 = "Change"
    change_dvof_points = working_df.FID.tolist()
    ignore_tracker = ignore_tds(ignore_tracker, change_dvof_points, working_df)
    change_dvof_points = identify_checks(change_dvof_points, out_checks, dvof_source, step_3)

    # After the final set of checks is identified, anything that remains is considered a Change.
    # A change is any DVOF point that differs from its corresponding TDS feature by a distance of
    # more than 35 meters but less than 65 meters and/or a height difference greater than 5 feet.
    working_df = working_df[working_df.FID.isin(change_dvof_points)]
    changes_tracker = compile_change_points(working_df)

    # Ensure that there were no DVOF points unaccounted for in processing.    
    checks_remain_df = working_df[~working_df.FID.isin(change_dvof_points)]
    if not checks_remain_df.empty:
        remains_points = checks_remain_df.FID.unique().tolist()
        ignore_tracker = ignore_tds(ignore_tracker, remains_points, checks_remain_df)
        write("Unknown points added to check file.")

    '''
        Post-Processing: All points identified as Change are added to the output. Then any TDS features
                         not identified as a match or check (i.e. not present in the ignore_tracker) are
                         populated as an Add ('A'). The regional codes, source dates, and defaults are then
                         filled in. Finally, the DVOF points not present in the processing table, which are
                         not close enough to a qualifying TDS feature, are populated as a Delete ('D').

                         After the output has been written, a basic report is printed for the user.
    '''
    write("\nProcessing Complete. Writing Output DVOF...")
    
    # Populate change points in output dvof file
    populate_edits(changes_tracker, out_dvof)
    
    # Look for adds in tds features (accounting for tds features already checked or added with change) and populate in output dvof file
    populate_adds(working_tds_list, ignore_tracker, out_dvof)

    # Populate the default and extra values in dvof after Adds and Edits
    fill_in_defaults_dvof(wac, country_boundaries, imagery_footprint, out_dvof)

    # Look for deletes in dvof points that did not appear in near table (thus is not close enough to a tds point) and populate in output dvof file
    populate_deletes(delete_list, dvof_source, out_dvof)


    write("\n\nAccounted for all DVOF Points? {}".format(checks_remain_df.empty))

    check_count = get_count(out_checks)
    write("\n{} checks written to Check File.\n\n".format(check_count))


if __name__=='__main__':
    ap.env.overwriteOutput = True
    argv = tuple(ap.GetParameterAsText(i) for i in range(ap.GetArgumentCount()))
    now = dt.datetime.now()
    main(*argv)
    write(dt.datetime.now() - now)