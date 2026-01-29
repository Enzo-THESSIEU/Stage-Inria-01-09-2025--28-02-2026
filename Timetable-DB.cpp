#include <windows.h>
#include <sqlext.h>
#include <iostream>
#include <string>
#include <cstring>
#include <cstdio>
#include "small_instance_generator.h"

int extract_from_db() {
    // Initialize ODBC environment
    SQLHENV env = nullptr;
    SQLHDBC dbc = nullptr;
    SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &env);
    SQLSetEnvAttr(env, SQL_ATTR_ODBC_VERSION, (void*)SQL_OV_ODBC3, 0);
    SQLAllocHandle(SQL_HANDLE_DBC, env, &dbc);

    // Option A: Windows auth (local SQL Server)
    // Example: SERVER=localhost;DATABASE=YourDb;Trusted_Connection=Yes;
    SQLWCHAR connStr[] =
        L"DRIVER={ODBC Driver 18 for SQL Server};SERVER=localhost;DATABASE=train_timetable;"
        L"Trusted_Connection=Yes;TrustServerCertificate=Yes;";

    // Option B: SQL Server auth
    SQLWCHAR outStr[1024];
    SQLSMALLINT outLen;
    SQLRETURN ret = SQLDriverConnectW(
        dbc, NULL, connStr, SQL_NTS, outStr, 1024, &outLen, SQL_DRIVER_COMPLETE
    );

    // Check the connection status
    if (!(ret == SQL_SUCCESS || ret == SQL_SUCCESS_WITH_INFO)) {
        std::cerr << "Connection failed\n";
        return 1;
    }

    // Allocate statement handle
    SQLHSTMT stmt = nullptr;
    SQLAllocHandle(SQL_HANDLE_STMT, dbc, &stmt);

    // Simple query
    SQLWCHAR sql[] = L"SELECT TOP 10 TrainId, TripName, DepartureTime, ArrivalTime FROM dbo.Timetable ORDER BY DepartureTime;";
    SQLExecDirectW(stmt, sql, SQL_NTS);

    // Fetch and display results
    SQLINTEGER id;
    SQLWCHAR trip[256];
    TIMESTAMP_STRUCT dep, arr;

    while (SQLFetch(stmt) == SQL_SUCCESS) {
        SQLGetData(stmt, 1, SQL_C_SLONG, &id, 0, NULL);
        SQLGetData(stmt, 2, SQL_C_WCHAR, trip, sizeof(trip), NULL);
        SQLGetData(stmt, 3, SQL_C_TYPE_TIMESTAMP, &dep, 0, NULL);
        SQLGetData(stmt, 4, SQL_C_TYPE_TIMESTAMP, &arr, 0, NULL);

        std::wcout << L"#" << id << L" " << trip
                   << L" dep " << dep.hour << L":" << dep.minute
                   << L" arr " << arr.hour << L":" << arr.minute << L"\n";
    }

    // Clean up
    SQLFreeHandle(SQL_HANDLE_STMT, stmt);
    SQLDisconnect(dbc);
    SQLFreeHandle(SQL_HANDLE_DBC, dbc);
    SQLFreeHandle(SQL_HANDLE_ENV, env);
    return 0;
}


// Helper function to print ODBC errors
void printOdbcError(SQLSMALLINT handleType, SQLHANDLE handle, const char* msg) {
    std::cerr << msg << "\n";
    SQLCHAR sqlState[6] = {0};
    SQLCHAR message[1024] = {0};
    SQLINTEGER nativeError = 0;
    SQLSMALLINT textLength = 0;

    SQLSMALLINT i = 1;
    while (SQLGetDiagRecA(handleType, handle, i, sqlState, &nativeError,
                         message, sizeof(message), &textLength) == SQL_SUCCESS) {
        std::cerr << "  [" << sqlState << "] (" << nativeError << ") " << message << "\n";
        ++i;
    }
}

std::string minutes_to_datetime(double minutes, const std::string& date = "2026-01-01") {
    int total = (int)minutes;
    int hh = total / 60;
    int mm = total % 60;

    char buf[32];
    std::snprintf(buf, sizeof(buf), "%s %02d:%02d:00", date.c_str(), hh, mm);
    return std::string(buf);
}

int insert_train_into_db(const std::vector<Timetabled_Train_Generator::Train>& all_trains) {
    const std::string connStr =
        "DRIVER={ODBC Driver 18 for SQL Server};"
        "SERVER=.\\SQLEXPRESS;"
        "DATABASE=darp;"
        "Trusted_Connection=Yes;"
        "Encrypt=Yes;"
        "TrustServerCertificate=Yes;";

    // Server=localhost\SQLEXPRESS;Database=master;Trusted_Connection=True;

    SQLHENV env = SQL_NULL_HENV;
    SQLHDBC dbc = SQL_NULL_HDBC;
    SQLHSTMT stmt = SQL_NULL_HSTMT;

    // Init
    if (SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &env) != SQL_SUCCESS) return 1;
    if (SQLSetEnvAttr(env, SQL_ATTR_ODBC_VERSION, (SQLPOINTER)SQL_OV_ODBC3, 0) != SQL_SUCCESS) return 1;
    if (SQLAllocHandle(SQL_HANDLE_DBC, env, &dbc) != SQL_SUCCESS) return 1;

    // Connect
    SQLCHAR outConn[1024];
    SQLSMALLINT outLen = 0;
    SQLRETURN ret = SQLDriverConnectA(
        dbc, NULL,
        (SQLCHAR*)connStr.c_str(), SQL_NTS,
        outConn, sizeof(outConn), &outLen,
        SQL_DRIVER_NOPROMPT
    );

    if (!(ret == SQL_SUCCESS || ret == SQL_SUCCESS_WITH_INFO)) {
        printOdbcError(SQL_HANDLE_DBC, dbc, "Connection failed");
        return 1;
    }

    // // Start transaction
    // SQLExecDirectA((SQLHSTMT)dbc, (SQLCHAR*)"BEGIN TRAN;", SQL_NTS); // quick&dirty
    // // Better: allocate stmt and exec "BEGIN TRAN;" there (shown below)

    // Allocate statement
    if (SQLAllocHandle(SQL_HANDLE_STMT, dbc, &stmt) != SQL_SUCCESS) return 1;

    SQLExecDirectA(stmt, (SQLCHAR*)"TRUNCATE TABLE dbo.Timetable;", SQL_NTS);

    // Begin tran properly on stmt
    SQLExecDirectA(stmt, (SQLCHAR*)"BEGIN TRAN;", SQL_NTS);

    // Prepare insert
    const char* sql =
        "INSERT INTO dbo.Timetable (TripName, DepartureTime, ArrivalTime) "
        "VALUES (?, ?, ?);";

    if (SQLPrepareA(stmt, (SQLCHAR*)sql, SQL_NTS) != SQL_SUCCESS) {
        printOdbcError(SQL_HANDLE_STMT, stmt, "SQLPrepare failed");
        return 1;
    }

    // Buffers for parameters (reused each loop)
    char tripBuf[256];
    char depBuf[32];
    char arrBuf[32];
    SQLLEN tripLen = SQL_NTS, depLen = SQL_NTS, arrLen = SQL_NTS;

    // Bind once (buffers updated per row)
    SQLBindParameter(stmt, 1, SQL_PARAM_INPUT, SQL_C_CHAR, SQL_VARCHAR, 200, 0, tripBuf, 0, &tripLen);
    SQLBindParameter(stmt, 2, SQL_PARAM_INPUT, SQL_C_CHAR, SQL_VARCHAR, 30,  0, depBuf,  0, &depLen);
    SQLBindParameter(stmt, 3, SQL_PARAM_INPUT, SQL_C_CHAR, SQL_VARCHAR, 30,  0, arrBuf,  0, &arrLen);

    for (const auto& train : all_trains){
        // TripName
        std::snprintf(tripBuf, sizeof(tripBuf), "%s", train.stations.c_str());

        // Times
        std::string dep = minutes_to_datetime(train.departure);
        std::string arr = minutes_to_datetime(train.arrival);

        std::snprintf(depBuf, sizeof(depBuf), "%s", dep.c_str());
        std::snprintf(arrBuf, sizeof(arrBuf), "%s", arr.c_str());

        ret = SQLExecute(stmt);
        if (!(ret == SQL_SUCCESS || ret == SQL_SUCCESS_WITH_INFO)) {
            printOdbcError(SQL_HANDLE_STMT, stmt, "SQLExecute failed (insert)");
            SQLExecDirectA(stmt, (SQLCHAR*)"ROLLBACK;", SQL_NTS);
            return 1;
        }
    }


    // Commit
    SQLExecDirectA(stmt, (SQLCHAR*)"COMMIT;", SQL_NTS);

    // Cleanup
    SQLFreeHandle(SQL_HANDLE_STMT, stmt);
    SQLDisconnect(dbc);
    SQLFreeHandle(SQL_HANDLE_DBC, dbc);
    SQLFreeHandle(SQL_HANDLE_ENV, env);
    return 0;
}

 