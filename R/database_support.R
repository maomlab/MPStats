


#' Construct a database from database connection parameters
#'
#'    It's recommended to store database connection parameters in the
#'    <project>/parameters.R file 
#'
#' @param database_parameters List with the following elements
#'    database_type, host, port, user, password, schema
#'
#' @return DBI database connect with the given schema specified
#' 
#' @export
get_database_connection <- function(database_parameters){
    if(database_parameters$database_type == "mysql"){
        if(!requireNamespace("RMySQL")){
            stop(paste0(
              "Requesting to load a mysql database, ",
              "please install the 'RMySQL' package with ",
              "'install.packages(\"RMySQL\")'."))
            return(invisible())
        }
        con <- DBI::dbConnect(
            RMySQL::MySQL(),
            host = database_parameters$host,
            port = database_parameters$port,
            user = database_parameters$user,
            password = database_parameters$password) %>%
            set_schema(database_parameters$schema)
    } else {
        stop(paste0(
          "Unrecognized or missing database_type ",
          database_parameters$database_type))
    }
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
