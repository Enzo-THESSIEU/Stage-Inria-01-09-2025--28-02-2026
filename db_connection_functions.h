#pragma once

// Necessary Libraries
#include <string>
#include <windows.h>
#include <sqlext.h>

// Structures
#include "structures.h"

// Different classes in Header

class odbc_functions {
public:
    // Build: "Driver=...;Server=...;Database=...;"
    std::string produce_conn_str(const std::string& server,
                                 const std::string& database) const;

    // Alloc ENV/DBC/STMT + connect, returns handles
    sql_params connect_to_db(const std::string& connStr, sql_params sql_param = {}) const;

    // Free STMT + disconnect + free DBC/ENV
    sql_params cleanup_sql_query(sql_params sql_param) const;

    // Reset the statement so it can be reused for a different query
    sql_params free_up_statement_for_new_query_type(sql_params sql_param) const;

private:
    static std::string odbc_diag(SQLSMALLINT handleType, SQLHANDLE handle);
    static void check_rc(SQLRETURN rc, const char* where);
    static void check_rc_diag(SQLRETURN rc, const char* where, SQLSMALLINT handleType, SQLHANDLE handle);
};