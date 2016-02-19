#pragma once

#include <stack>

#define RECT_IDS    1024
#define CONFIG_IDS  1024

namespace som {

// List of ID's free for use
extern std::stack<int> rectangle_ids;
extern std::stack<int> config_ids;

int get_rectangle_id();
int get_config_id();

void release_rectangle_id(int id);
void release_config_id(int id);

}
