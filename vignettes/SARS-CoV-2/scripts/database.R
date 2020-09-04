library(plyr)
library(tidyverse)
library(RMySQL)
library(magrittr)
library(tictoc)


get_test_database_connection <- function(){
    con <- DBI::dbConnect(
        RMySQL::MySQL(),
        host = "covid19cp.cgymeokijgns.us-east-1.rds.amazonaws.com",
        port = 3306,
        user = "covid19cp",
        password = "Genes-brett-Flip-9Bottling")
    con %>% DBI::dbSendQuery("USE covid19cp")
    con
}

set_schema <- function(con, schema) {
    available_schemas <- con %>% DBI::dbGetQuery("SHOW SCHEMAS;")
    if (schema %in% available_schemas$Database) {
        con %>% DBI::dbSendQuery(paste0("USE ", schema))
    } else{
        cat(
            "Unrecognized schema, available schemas are:\n  ",
            available_schemas$Database %>% paste0(collapse = "\n  "), "\n",
            sep = "")
    }
    con
}

# TODO: gather database arguments to be parsed by argparse
# then use these arguments to initialize a database connection
add_mysql_argparse_parameters <- function(parser){
    parser$add_argument(
        "--database_driver",
        action = "store",
        default = "RMySQL",
        dest = "database_driver",
        type = "character")
    parser$add_argument()
        
}

        


get_primary_database_connection <- function(schema="covid19primary") {
    con <- DBI::dbConnect(
        RMySQL::MySQL(),
        host = "covid19primary.cgymeokijgns.us-east-1.rds.amazonaws.com",
        port = 3306,
        user = "covid19primary",
        password = "atop-9most-5Inn-Dandruff9")
#        password = "frightful-bootlace-Cats-flaring8")
    con %>% set_schema(schema)
    con
}


build_indices <- function(plate_id) {
    # if a table is on the LHS of a join
    # an index over the join columns makes the query faster
    con %>% DBI::dbSendQuery(paste0("
       CREATE UNIQUE INDEX SARS_", plate_id, "_Per_Image_index
       ON ", plate_prefix, "_Per_Image (ImageNumber ASC);"))
    
    con %>% DBI::dbSendQuery(paste0("
       CREATE UNIQUE INDEX SARS_", plate_id, "_Per_syncytia_index
       ON ", plate_id, "_Per_syncytia (
            ImageNumber ASC,
            syncytia_Number_Object_Number ASC);"))     
    # 27,299 rows
    
    con %>% DBI::dbSendQuery(paste0("
       CREATE UNIQUE INDEX", plate_id, "_Per_Nuclei_index
       ON ", plate_id, "_Per_Nuclei (
            ImageNumber ASC,
            Nuclei_Number_Object_Number ASC);"))
    # <MySQLResult:38546520,0,19>
    
    con %>% DBI::dbSendQuery(paste0("
       CREATE UNIQUE INDEX ", plate_id, "_Per_Cells_index
       ON ", plate_id, "_Per_Cells (
            ImageNumber ASC,
            Cells_Number_Object_Number ASC);"))
    # <MySQLResult:135600456,0,20>
}

#######################


get_acas_database <- function(stage=FALSE){
    con <- DBI::dbConnect(
        RPostgres::Postgres(),
        host=paste0("umsexton", ifelse(stage, "-stage", ""), ".onacaslims.com"),
        port=5432,
        user="acas",
        password="acas")
    con %>% DBI::dbSendQuery("SET search_path TO acas;")
    con
}



"lsthing_tag"
"ls_thing_state"
"ls_tag"
"ls_seq_trt_grp"
"ls_seq_protocol"
"ls_seq_itx_subj_cntr"
"ls_thing_label"
"ls_seq_itx_protocol_protocol"
"ls_interaction"
"ls_transaction"
"ls_seq_subject"
"ls_seq_itx_cntr_cntr"
"ls_thing"
"ls_seq_itx_expt_expt"
"ls_seq_expt"
"ls_role"

"role_kind"
"role_type"
"subject_state"
"treatment_group_value"
"operator"
"subject_label"
"vw_experiment"
"v_api_dv_experiment"
"vw_api_dv_experiment"
"vw_dv_expt_ag_kinds"
"v_api_dv_protocol"
"itx_subject_container"
"experiment_analysisgroup"
"container_value"
"file_type"
"experiment"
"analysis_group_label"
"analysisgroup_treatmentgroup"
"itx_subject_container_value"
"analysis_group_state"
"salt_form"
"analysis_group_value"
"protocol_type"
"protocol_tag"
"parent_annotation"
"physical_state"
"state_kind"
"solution_unit"
"state_type"
"stereo_category"
"salt_loader"
"salt_form_alias_type"
"salt_form_alias"
"salt_mol_idx_shadow_hash"
"dry_run_compound_mol_idx_shadow"
"salt"
"isotope"
"code_kind"
"bulk_load_template"
"author_label"
"author_value"
"author"
"author_state"
"iso_salt"
"bulk_load_file"
"author_role"
"compound"
"itx_container_container_state"
"itx_expt_expt"
"interaction_type"
"interaction_kind"
"file_thing"
"experiment_type"
"dry_run_compound_mol_idx_shadow_hash"
"dry_run_compound"
"itx_expt_expt_state"
"itx_container_container_value"
"unit_kind"
"parent_mol_idx_shadow"
"parent_mol_idx_shadow_hash"
"analysis_group"
"schema_version"
"salt_form_alias_kind"
"pre_def_corp_name"
"parent"
"itx_protocol_protocol_state"
"itx_protocol_protocol"
"itx_expt_expt_value"
"thing_kind"
"uncertainty_kind"
"unit"
"treatment_group_label"
"thing_type"
"thing_page_archive"
"thing_page"
"itx_ls_thing_ls_thing_value"
"itx_ls_thing_ls_thing_state"
"itx_ls_thing_ls_thing"
"ddict_value"
"experiment_tag"
"experiment_kind"
"experiment_label"
"corp_name"

"value_type"
"unit_type"
"treatment_group"
"code_type"
"ddict_type"
"ddict_kind"

"cron_job"
"container_type"
"container_kind"
"compound_type"
"code_origin"
"file_list"
"itx_container_container"
"container_label"
"vendor"
"salt_form_mol_idx_shadow_hash"
"lot"
"cmpd_reg_app_setting"
"treatmentgroup_subject"
"update_log"
"salt_form_mol_idx_shadow"

"label_sequence"
"label_kind"
"itx_subject_container_state"
"itx_protocol_protocol_value"
"subject_value"
"subject"
"container_state"
"container"
"application_setting"
"purity_measured_by"
"protocol_value"
"protocol_state"
"protocol_label"
"protocol"
"treatment_group_state"
"experiment_value"
"experiment_state"
"temp_select_table"
"ls_seq_container"
"ls_seq_anl_grp"
"qc_compound"
"parent_alias_type"
"parent_alias_kind"
"parent_alias"
"operator_type"
"v_api_dv_ag_results"
"dv_api_all_ag_results"
"vw_expt_tested_lots"
"vw_expt_columns"
"vw_protocol_conditions"
"vw_dv_ag_results"
"api_system_statistics"
"api_subject_container_results"
"vw_assay_tree"
"batch_code_experiment_links"
"api_all_data"
"api_batch_cmpd_reg_links"
"application_paths"
"api_analysis_group_results"
"api_protocol"
"api_experiment"
"p_api_analysis_group_results"
"ls_thing_value"
"lot_alias_type"
"lot_alias_kind"
"lot_alias"
"label_type"
"label_sequence_ls_role"
"value_kind"
"protocol_kind"
"operator_kind"
"salt_mol_idx_shadow"
"api_dose_response"
"api_curve_params"
"api_subject_results"
"api_salt_iso_salt"
"api_isotope_iso_salt"
"api_hts_treatment_results"
"api_container_contents"
"api_batch_parentname"
"api_file_list"
"api_experiment_results"
"api_experiment_approved"

