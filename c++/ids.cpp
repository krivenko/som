#include "ids.hpp"

#include <stack>

namespace som {

std::stack<int> rectangle_ids = [](int size) {
 std::stack<int> tmp;
 for(int i = size - 1; i >= 0; --i) tmp.push(i);
 return tmp;
}(RECT_IDS);

std::stack<int> config_ids = [](int size) {
 std::stack<int> tmp;
 for(int i = size - 1; i >= 0; --i) tmp.push(i);
 return tmp;
}(CONFIG_IDS);

int get_rectangle_id() { int id = rectangle_ids.top(); rectangle_ids.pop(); return id; }
int get_config_id() { int id = config_ids.top(); config_ids.pop(); return id; }

void release_rectangle_id(int id) { rectangle_ids.push(id); }
void release_config_id(int id) { config_ids.push(id); }

}
