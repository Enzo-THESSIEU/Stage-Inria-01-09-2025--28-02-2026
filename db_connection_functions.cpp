// Inclusions
#include "db_connection_functions.h"

#include <windows.h>
#include <sqlext.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <sstream>
#include <cstdio>


// Checker functions
std::string odbc_functions::odbc_diag(SQLSMALLINT handleType, SQLHANDLE handle) {
    SQLCHAR state[6]{};
    SQLCHAR msg[1024]{};
    SQLINTEGER nativeErr = 0;
    SQLSMALLINT msgLen = 0;

    if (SQLGetDiagRecA(handleType, handle, 1, state, &nativeErr, msg, sizeof(msg), &msgLen) == SQL_SUCCESS) {
        return std::string((char*)state) + " | " + std::string((char*)msg);
    }
    return "No diagnostics";
}

void odbc_functions::check_rc(SQLRETURN rc, const char* where) {
    if (rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO) return;
    throw std::runtime_error(std::string("ODBC error at: ") + where);
}

void odbc_functions::check_rc_diag(SQLRETURN rc, const char* where, SQLSMALLINT handleType, SQLHANDLE handle) {
    if (rc == SQL_SUCCESS || rc == SQL_SUCCESS_WITH_INFO) return;
    throw std::runtime_error(std::string(where) + " | " + odbc_diag(handleType, handle));
}



// ---------------- ODBC Functions -------------------------
std::string odbc_functions::produce_conn_str(const std::string& server, const std::string& database) const {
    // Note: Encrypt=no for local dev; if you want Encrypt=yes, add TrustServerCertificate=yes
    return "Driver={ODBC Driver 18 for SQL Server};"
           "Server=" + server + ";"
           "Database=" + database + ";"
           "Trusted_Connection=yes;"
           "Encrypt=no;";
};


sql_params odbc_functions::connect_to_db(const std::string& connStr, sql_params sql_param) const {
    // Allocate ENV
    SQLRETURN rc = SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &sql_param.env);
    check_rc(rc, "SQLAllocHandle(SQL_HANDLE_ENV)");

    rc = SQLSetEnvAttr(sql_param.env, SQL_ATTR_ODBC_VERSION, (SQLPOINTER)SQL_OV_ODBC3, 0);
    check_rc(rc, "SQLSetEnvAttr(SQL_ATTR_ODBC_VERSION)");

    // Allocate DBC
    rc = SQLAllocHandle(SQL_HANDLE_DBC, sql_param.env, &sql_param.dbc);
    check_rc(rc, "SQLAllocHandle(SQL_HANDLE_DBC)");

    // Connect
    rc = SQLDriverConnectA(
        sql_param.dbc, nullptr,
        (SQLCHAR*)connStr.c_str(), SQL_NTS,
        nullptr, 0, nullptr,
        SQL_DRIVER_NOPROMPT
    );
    // If connect fails, diagnostics are on the DBC handle
    check_rc_diag(rc, "SQLDriverConnectA", SQL_HANDLE_DBC, sql_param.dbc);

    // Allocate STMT
    rc = SQLAllocHandle(SQL_HANDLE_STMT, sql_param.dbc, &sql_param.stmt);
    check_rc(rc, "SQLAllocHandle(SQL_HANDLE_STMT)");

    return sql_param;
}


sql_params odbc_functions::cleanup_sql_query(sql_params sql_param) const {
    if (sql_param.stmt != SQL_NULL_HSTMT) {
        SQLFreeHandle(SQL_HANDLE_STMT, sql_param.stmt);
        sql_param.stmt = SQL_NULL_HSTMT;
    }
    if (sql_param.dbc != SQL_NULL_HDBC) {
        SQLDisconnect(sql_param.dbc);
        SQLFreeHandle(SQL_HANDLE_DBC, sql_param.dbc);
        sql_param.dbc = SQL_NULL_HDBC;
    }
    if (sql_param.env != SQL_NULL_HENV) {
        SQLFreeHandle(SQL_HANDLE_ENV, sql_param.env);
        sql_param.env = SQL_NULL_HENV;
    }
    return sql_param;
}


sql_params odbc_functions::free_up_statement_for_new_query_type(sql_params sql_param) const {
    if (sql_param.stmt == SQL_NULL_HSTMT) return sql_param;

    SQLFreeStmt(sql_param.stmt, SQL_CLOSE);        // close any open cursor/results
    SQLFreeStmt(sql_param.stmt, SQL_UNBIND);       // unbind result columns
    SQLFreeStmt(sql_param.stmt, SQL_RESET_PARAMS); // unbind parameters

    return sql_param;
}