// define header for timetable database extraction code
#ifndef Timetable_DB_h
#define Timetable_DB_h

// Include Libraries
#include <windows.h>
#include <sqlext.h>
#include <iostream>
#include <string>
#include "small_instance_generator.h"

// Include the different functions used
int extract_from_db();
void printOdbcError(SQLSMALLINT handleType, SQLHANDLE handle, const char* msg);
std::string minutes_to_datetime(double minutes, const std::string& date = "2026-01-01");
int insert_train_into_db(const std::vector<Timetabled_Train_Generator::Train>& all_trains);

#endif