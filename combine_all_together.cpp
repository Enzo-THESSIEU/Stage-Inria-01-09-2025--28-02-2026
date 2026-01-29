#include "small_instance_generator.h"
#include "Timetable-DB.h"
#include <vector>

int main() {
    Timetabled_Train_Generator gen;
    auto all_trains = gen.get_all_trains();
    insert_train_into_db(all_trains);
    
    return 0;
    }