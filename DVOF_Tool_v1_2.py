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
    feature_layer_list = []
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
            feature_layer_list.append(srf_name)

        # Create Point Layer and check count
        ap.management.MakeFeatureLayer(pnt_fc, pnt_name, where_clause)
        pnt_count = get_count(pnt_name)
        if pnt_count == 0:
            write("No {} features above 46 meters present. Skipping...".format(pnt_fc))
        else:
            write("{0} {1} features above 46 meters found.".format(pnt_count, pnt_fc))
            feature_layer_list.append(pnt_name)

    return feature_layer_list
#-------------------------------------------------------------------------------------------------------------------------------------------------

'''
    This is an unnecessary step but leaving the function here for the moment as it may be useful in other programs.
'''
def convert_tds_srf_to_pnt(layer_dict):
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

def get_source_oid(source_dvof, oid_name=[]):
    if oid_name == []:
        field_list = ap.ListFields(source_dvof, field_type='OID')
        name = field_list[0].name
        oid_name.append(name)
        return name
    else:
        return oid_name[0]
#-------------------------------------------------------------------------------------------------------------------------------------------------

def limit_dvof(source_dvof):
    desc_tds = ap.Describe(ap.env.workspace)
    tds_extent = desc_tds.extent
    tds_extent_poly = tds_extent.polygon

    extent_check_name = "DVOF_extent_check_lyr"

    ap.management.MakeFeatureLayer(source_dvof, extent_check_name)

    ap.management.SelectLayerByLocation(extent_check_name, "WITHIN", tds_extent_poly)

    within_ext_count = get_count(extent_check_name)
    if within_ext_count == 0:
        return None
    
    return_lyr_name = "extent_dvof_lyr"
    ap.management.MakeFeatureLayer(extent_check_name, return_lyr_name)

    ap.management.Delete(extent_check_name)

    return return_lyr_name
#-------------------------------------------------------------------------------------------------------------------------------------------------

def check_source_dvof(dvof_data):
    oid_field = get_source_oid(dvof_data)
    total_count = get_count(dvof_data)
    lyr_name = "working_dvof_lyr"
    if get_field_type(dvof_data, 'HEIGHTAGL') == 'String':
        working_fids = find_dvof_by_hgt(dvof_data)
        if len(working_fids) == 1:
            where_clause = """{0} = {1}""".format(fieldDelim(dvof_data, oid_field), working_fids[0])
        else:
            where_clause = """{0} IN {1}""".format(fieldDelim(dvof_data, oid_field), tuple(working_fids))
        ap.management.MakeFeatureLayer(dvof_data, lyr_name, where_clause)
    else:
        where_clause = """{0} >= {1}""".format(fieldDelim(dvof_data, 'HEIGHTAGL'), HGT_FT)
        ap.management.MakeFeatureLayer(dvof_data, lyr_name, where_clause)
    work_count = get_count(lyr_name)
    if work_count == 0:
        return None
    write("{0} of {1} DVOF Points meet height threshold.\nChecking qualifying points against TDS.".format(work_count, total_count))

    return lyr_name
#-------------------------------------------------------------------------------------------------------------------------------------------------

def get_field_type(data, field_name):
    fields = ap.ListFields(data)
    for i in fields:
        if i.name == field_name:
            return i.type
    return None
#-------------------------------------------------------------------------------------------------------------------------------------------------

def find_dvof_by_hgt(dvof_data):
    fid_list = []
    fields = ['OID@', 'HEIGHTAGL']
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
    oid_field = get_source_oid(in_dvof)
    cols = [oid_field, 'FEATURETYP', 'HEIGHTAGL', 'TYPENAME']
    write(cols)
    dataframe = pd.DataFrame(ap.da.FeatureClassToNumPyArray(in_table=in_dvof, field_names=cols, skip_nulls=False, null_value=-999999))
    dataframe.HEIGHTAGL = dataframe.HEIGHTAGL.astype(float)
    if oid_field != 'FID':
        dataframe.rename(columns={oid_field: 'FID'}, inplace=True)
    write("Source DVOF:\n\n")
    write(dataframe)

    return dataframe
#-------------------------------------------------------------------------------------------------------------------------------------------------

def fix_fc_names(row):
    row.NEAR_FC = os.path.basename(row.NEAR_FC)
    #if row.NEAR_FC.count('in_memory//') > 0:
    #    row.NEAR_FC = row.NEAR_FC.strip('in_memory//')
    return row
#-------------------------------------------------------------------------------------------------------------------------------------------------

def meters_to_feet(row):
    row.hgt = row.hgt / CONVERSION_FACTOR
    return row
#-------------------------------------------------------------------------------------------------------------------------------------------------

def check_table_for_tds_layer(tds_list, dataframe):
    near_fc_list = dataframe.NEAR_FC.unique().tolist()
    check_list = []
    for i in tds_list:
        if i.count('in_memory//') > 0:
            corrected = i.strip('in_memory//')
            check_list.append(corrected)
        else:
            check_list.append(i)

    tds_layers_match = all(fc in check_list for fc in near_fc_list)
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
        write("\n{0} {1} points identified for inspection.".format(len(checks_list), step_name))
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
    fields = ['OID@','FEATURETYP', 'TYPENAME', 'SHAPE@']
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

    oid_match_dict = pd.Series(dataframe.FID.values, index=dataframe.NEAR_FID).to_dict()

    return change_dict, oid_match_dict
#-------------------------------------------------------------------------------------------------------------------------------------------------

def populate_adds(tds_list, ignore_dict, dvof_fc, check_db):
    add_list = []
    tds_fields = ["OID@", "fcsubtype", "hgt", "SHAPE@XY"]
    otherremar = 'A'
    voidentifi = -999999
    horizacc = '25'

    for lyr in tds_list:
        
        working_layer = lyr
        if len(ignore_dict[lyr]) == 1:
            where_clause = """{0} <> {1}""".format(fieldDelim(working_layer, 'objectid'), ignore_dict[lyr][0])
        elif len(ignore_dict[lyr]) > 1:
            where_clause = """{0} NOT IN {1}""".format(fieldDelim(working_layer, 'objectid'), tuple(ignore_dict[lyr]))
        else:
            where_clause = ''
        with ap.da.SearchCursor(working_layer, tds_fields, where_clause) as cursor:
            for row in cursor:
                geometry = row[-1]
                long_x = geometry[0]
                lat_y = geometry[1]
               
                featuretyp = get_featuretyp_code(row[1])
                # Check that there is a featuretyp
                if featuretyp == 0:
                    check_row = ['', '', 'Fix DVOF Output: {} does not have FEATURETYP assigned.'.format(row[1]), row[-1]]
                    move_to_checkfile(check_db, check_row)
                hgt_feet = round(row[2] / CONVERSION_FACTOR, 3)
                insert_row = [voidentifi, lat_y, long_x, featuretyp, hgt_feet, horizacc, otherremar, geometry]
                add_list.append(row[0])
                move_to_delivery_dvof(dvof_fc, insert_row)
    write("\n{} Add Points populated:".format(len(add_list)))
#-------------------------------------------------------------------------------------------------------------------------------------------------

def populate_edits(edits_dict, dvof_fc, check_db, oid_match_dict, source_dvof):
    add_list = []
    tds_fields = ["OID@", "fcsubtype", "hgt", "SHAPE@XY"]
    otherremar = 'C'
    


    for lyr in edits_dict:
       
        working_layer = lyr
        if len(edits_dict[lyr]) == 1:
            where_clause = """{0} = {1}""".format(fieldDelim(working_layer, 'objectid'), edits_dict[lyr][0])
        else:
            where_clause = """{0} IN {1}""".format(fieldDelim(working_layer, 'objectid'), tuple(edits_dict[lyr]))
        with ap.da.SearchCursor(working_layer, tds_fields, where_clause) as cursor:
            for row in cursor:
                geometry = row[-1]
                long_x = geometry[0]
                lat_y = geometry[1]
                source_oid = oid_match_dict[row[0]]
                voidentifi = get_source_value(source_oid, "VOIDENTIFI", source_dvof)
                horizacc = get_source_value(source_oid, "HORIZACC", source_dvof)
                hgt_feet = round(row[2] / CONVERSION_FACTOR, 3)
                featuretyp = get_featuretyp_code(row[1])
                # Check that there is a featuretyp
                if featuretyp == 0:
                    check_row = ['', '', 'Fix DVOF Output: {} does not have FEATURETYP assigned.'.format(row[1]), row[-1]]
                    move_to_checkfile(check_db, check_row)
                insert_row = [voidentifi, lat_y, long_x, featuretyp, hgt_feet, horizacc, otherremar, geometry]
                add_list.append(row[0])
                move_to_delivery_dvof(dvof_fc, insert_row)
    write("\n{} Change Points populated:".format(len(add_list)))

#-------------------------------------------------------------------------------------------------------------------------------------------------
# Function will take in a table, a field, and an oid and return the field value
def get_source_value(source_oid, field_name, source_dvof_table):
    oid_field = get_source_oid(source_dvof_table)
    where_clause = """{0} = {1}""".format(fieldDelim(source_dvof_table, oid_field), source_oid)
    with ap.da.SearchCursor(source_dvof_table, [field_name], where_clause) as cursor:
        for row in cursor:
            field_value = row[0]
    
    return field_value


#-------------------------------------------------------------------------------------------------------------------------------------------------

def move_to_delivery_dvof(delivery_fc, row):
    delivery_fields = ['VOIDENTIFI', 'LATITUDE', 'LONGITUDE', 'FEATURETYP', 'HEIGHTAGL', 'HORIZACC', 'OTHERREMAR', 'SHAPE@']
    with ap.da.InsertCursor(delivery_fc, delivery_fields) as iCur:
        iCur.insertRow(row)
#-------------------------------------------------------------------------------------------------------------------------------------------------

def fill_in_defaults_dvof(wac_file, country_bndry_file, imagery_footprint, delivery_fc, check_db):
    dvof_fields = ['COUNTRYCD', 'WAC', 'SOURCEDT', 'HEIGHTAMSL',
                    'SINGLEMULT', 'MULTIPLENO', 'VALIDATION', 'HGTAGLACC', 'HGTAMSLACC', 'SECURITYCL', 'SURFACEMAT',
                    'HORIZDATUM', 'PRODUCER', 'PROCESSCD', 'LIGHTINGCD', 'FEATURETYP', 'SHAPE@']
    wac_fields = ['wac_index', 'SHAPE@']
    country_fields = ['CNTRY_CD', 'SHAPE@']
    imagery_fields= ['IngestDate', 'SHAPE@']

    # Default values
    tds_null = -999999
    countrycd = ''              # 0

    sourcedt = ''               # 2
    heightamsl = tds_null       # 3
    singlemult = 'S'            # 4
    multipleno = 1              # 5
    validation = 3              # 6
    #horizacc = tds_null         # 7     -- Defunct
    hgtaglacc = tds_null        # 8
    hgtamslacc = tds_null       # 9
    securitycl = 'U'            # 10
    surfacemat = 'U'            # 11
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
            # Check for country code. Add to checkfile if not present
            if countrycd == '':
                feat_typ = row[-2]
                if feat_typ == 122:
                    feat_typ = 123
                    row[-2] = feat_typ
                elif feat_typ == 437:
                    feat_typ = 438
                    row[-2] = feat_typ
                check_row = [feat_typ, '', "Point requires country code.", row[-1]]
                move_to_checkfile(check_db, check_row)

            # Assign and update row values
            row[0] = countrycd
            row[1] = wac_value
            row[2] = first_of_month(sourcedt)
            row[3] = heightamsl
            row[4] = singlemult
            row[5] = multipleno
            row[6] = validation
            #row[7] = horizacc
            row[7] = hgtaglacc
            row[8] = hgtamslacc
            row[9] = securitycl
            row[10] = surfacemat
            row[11] = horizdatum
            row[12] = producer
            row[13] = processcd
            row[14] = lightingcd

            uCur.updateRow(row)
    write("\nCoutnry Code, WAC, Source Date, and Defaults Populated.")

#-------------------------------------------------------------------------------------------------------------------------------------------------

def populate_deletes(deletes_list, source_dvof, delivery_dvof, imagery_footprint):
    fixed_deletes = [int(x) for x in deletes_list]
    
    oid_field = get_source_oid(source_dvof)

    if len(deletes_list) == 1:
        where_clause = """{0} = {1}""".format(fieldDelim(source_dvof, oid_field), fixed_deletes[0])
    else:
        where_clause = """{0} IN {1}""".format(fieldDelim(source_dvof, oid_field), tuple(fixed_deletes))

    source_fields = ['VOIDENTIFI', 'LATITUDE', 'LONGITUDE', 'FEATURETYP',
                'HEIGHTAGL', 'COUNTRYCD', 'WAC', 'HEIGHTAMSL',
                'SINGLEMULT', 'MULTIPLENO', 'VALIDATION', 'HORIZACC',
                'HGTAGLACC', 'HGTAMSLACC', 'SECURITYCL', 'SURFACEMAT', 'HORIZDATUM',
                'PROCESSCD', 'LIGHTINGCD', 'SHAPE@']

    out_fields = ['VOIDENTIFI', 'LATITUDE', 'LONGITUDE', 'FEATURETYP',
                'HEIGHTAGL', 'COUNTRYCD', 'WAC', 'HEIGHTAMSL',
                'SINGLEMULT', 'MULTIPLENO', 'VALIDATION', 'HORIZACC',
                'HGTAGLACC', 'HGTAMSLACC', 'SECURITYCL', 'SURFACEMAT', 'HORIZDATUM',
                'PROCESSCD', 'LIGHTINGCD', 'SHAPE@','OTHERREMAR', 'SOURCEDT']

    # UPDATE: Include SOURCEDT fix
    imagery_fields= ['IngestDate', 'SHAPE@']
    with ap.da.SearchCursor(source_dvof, source_fields, where_clause) as del_cur:
        for row in del_cur:
            with ap.da.SearchCursor(imagery_footprint, imagery_fields) as f_cur:
                for f_row in f_cur:
                    if not row[-1].disjoint(f_row[-1]):
                        if sourcedt == '':
                            sourcedt = f_row[0]
                        elif sourcedt < f_row[0]:
                            sourcedt = f_row[0]
                fixed_date = first_of_month(sourcedt)
            with ap.da.InsertCursor(delivery_dvof, out_fields) as in_cur:
                new_row = [x for x in row]
                new_row.append('D')
                new_row.append(fixed_date)
                in_cur.insertRow(new_row)

    write("\n{} Delete points populated.".format(len(deletes_list)))

#-------------------------------------------------------------------------------------------------------------------------------------------------

def get_featuretyp_code(tds_subtype):
    fcsubtype_featuretyp_dict = {
        100053: 350, 100295: 282, 100012: 116, 100161: 260, 100163: 260, 100164: 260,
        100083: 99, 100013: 121, 100167: 702, 100026: 183, 100028: 188, 100029: 186,
        100330: 925, 100200: 721, 100007: 170, 100018: 130, 100080: 111, 132626: 607,
        100134: 822, 100133: 821, 100033: 189, 100443: 533, 100105: 535, 100253: 765,
        100108: 680, 100177: 785, 100021: 137, 100558: 540, 100207: 717, 100004: 103,
        100072: 335, 100132: 122,    # Change 122 to 123 if it does not have country code
        100025: 182, 154703: 324, 100131: 820, 100139: 801, 100122: 501, 100337: 926, 100142: 824, 100045: 437,    # Change 437 to 438 if it does not have country code
        100044: 437    # Change 437 to 438 if it does not have country code
    }
    if tds_subtype in fcsubtype_featuretyp_dict.keys():
        return fcsubtype_featuretyp_dict[tds_subtype]
    else:
        return 0

#-------------------------------------------------------------------------------------------------------------------------------------------------
# Nat's date fix
def first_of_month(in_date):
    date_field = str(in_date)
    feat_date = dt.strptime(date_field, r"%m/%d/%Y")
    dumb_feat_date = feat_date.replace(day=1)
    return dumb_feat_date

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
    #project_aoi = argv[7]

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
    # Check booleans to skip functions if there are no TDS or DVOF features present,
    #   or stop the program if neither are present.
    tds_features_present = True
    dvof_features_present = True

    # Check TDS feature classes for qualifying features (hgt >= 46)
    tds_fc_list = check_tds_dvof_presence()
    if len(tds_fc_list) == 0:
        write("No qualifying TDS features detected.")
        tds_features_present = False  

    # Create means of tracking the features that need to be edited or added
    ignore_tracker = {lyr:[] for lyr in tds_fc_list}

    # Limit Source DVOF to the extent of the Delivery TDS
    dvof_limit_lyr = limit_dvof(dvof_source)
    if len(dvof_limit_lyr) is None:
        write("No Source DVOF Features detected within TDS Extent.")
        if not tds_features_present:
            write("This isn't the DVOF you're looking for. Move along.")
            sys.exit(0)
        dvof_features_present = False

    # Check Source DVOF for points above 46 m threshold. Anything below is ignored.
    working_dvof_lyr = check_source_dvof(dvof_limit_lyr)
    if working_dvof_lyr is None:
        write("No Source DVOF above threshold present.")
        if not tds_features_present:
            write("You're Done.")
            sys.exit(0)
        dvof_features_present = False
    
    # ******Temp Stop if either TDS features aren't present or DVOF features aren't, then stop the program******
    if not dvof_features_present or not tds_features_present:
        write("Unable to complete processing. Both DVOF and TDS features are required.")
        sys.exit(0)

    # Create the Near Table.
    # This is the primary table for this tool and will be added to and manipulated as needed.
    if dvof_features_present and tds_features_present:
        near_df = create_near_table(working_dvof_lyr, tds_fc_list, OUTER_DIST)


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
        table_check_approved = check_table_for_tds_layer(tds_fc_list, working_df)   # working_tds_list changed here to tds_fc_list, if surface test fails change back.
        if not table_check_approved:
            # Kill program
            #sys.exit(0)
            write("\nHold on to your butts...\n\n")

        # Pull the necessary information from the created layers into Pandas DataFrame
        tds_df_list = create_layer_dfs(tds_fc_list) # working_tds_list changed here to tds_fc_list, if surface test fails change back.

        # Combine the layer DataFrames into a single DataFrame
        combined_lyrs = combine_df_tables(tds_df_list)

        # Merge the layer data onto the Near table with the TDS OIDs as keys
        working_df = pd.merge(working_df, combined_lyrs,on='NEAR_FID', how='left')
    else:
        # Convert layer to Pandas DataFrame
        single_lyr_df = create_layer_dfs(tds_fc_list)   # working_tds_list changed here to tds_fc_list, if surface test fails change back.

        # Merge layer DataFrame with Working DataFrame
        working_df = pd.merge(working_df, single_lyr_df, on='NEAR_FID', how='left')

        # Add NEAR_FC Column and populate with the single fc, simplifies code in the 'Processing' phase
        working_df["NEAR_FC"] = tds_fc_list[0]  # working_tds_list changed here to tds_fc_list, if surface test fails change back.

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
    write("Looking for {} features.".format(step_1))

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
    write("Looking for {} features.".format(step_2))
    no_change_points = working_df[(working_df.NEAR_DIST.le(35)) & (working_df.hgt_delta < 5)].FID.tolist()
    ignore_tracker = ignore_tds(ignore_tracker, no_change_points, working_df)
    no_change_points = identify_checks(no_change_points, out_checks, dvof_source, step_2)

    # Remove the rows containing the 'no change' DVOF points.
    working_df = working_df[~working_df.FID.isin(no_change_points)]

    # Step 3 -- Filter any DVOF points that need to be checked and collect the rest as Change Points.
    step_3 = "Change"
    write("Looking for {} features.".format(step_3))
    change_dvof_points = working_df.FID.tolist()
    ignore_tracker = ignore_tds(ignore_tracker, change_dvof_points, working_df)
    change_dvof_points = identify_checks(change_dvof_points, out_checks, dvof_source, step_3)

    # After the final set of checks is identified, anything that remains is considered a Change.
    # A change is any DVOF point that differs from its corresponding TDS feature by a distance of
    # more than 35 meters but less than 65 meters and/or a height difference greater than 5 feet.
    working_df = working_df[working_df.FID.isin(change_dvof_points)]
    changes_tracker, oid_match_dict = compile_change_points(working_df)

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
    populate_edits(changes_tracker, out_dvof, out_checks, oid_match_dict, dvof_source)

    # Look for adds in tds features (accounting for tds features already checked or added with change) and populate in output dvof file
    populate_adds(tds_fc_list, ignore_tracker, out_dvof, out_checks)    # working_tds_list changed here to tds_fc_list, if surface test fails change back.

    # Populate the default and extra values in dvof after Adds and Edits
    fill_in_defaults_dvof(wac, country_boundaries, imagery_footprint, out_dvof, out_checks)

    # Look for deletes in dvof points that did not appear in near table (thus is not close enough to a tds point) and populate in output dvof file
    populate_deletes(delete_list, dvof_source, out_dvof, imagery_footprint)


    write("\n\nAccounted for all DVOF Points? {}".format(checks_remain_df.empty))

    check_count = get_count(out_checks)
    write("\n{} checks written to Check File.\n\n".format(check_count))


if __name__=='__main__':
    ap.env.overwriteOutput = True
    argv = tuple(ap.GetParameterAsText(i) for i in range(ap.GetArgumentCount()))
    now = dt.datetime.now()
    main(*argv)
    write(dt.datetime.now() - now)
