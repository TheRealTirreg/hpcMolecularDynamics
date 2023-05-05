#include <iostream>
#include <memory>

void set_to_one(int& value) {
    value = 1;
}

int main() {
    // int * pointer = new int;
    std::unique_ptr<int> pointer = std::make_unique<int>(0);  // guarantees delete before return
    set_to_one(*pointer);
    // if (pointer = nullptr)
    if (pointer == nullptr)
        std::cout << "pointer is nullptr\n";
    set_to_one(*pointer);
    return 0;
}