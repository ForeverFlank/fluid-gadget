#include <iostream>
#include "fixed16.cpp"

using namespace std;

int main()
{
    fixed16 a(2.5f);
    fixed16 b(1.5f);

    cout << "a = " << a.to_float() << "\n";
    cout << "b = " << b.to_float() << "\n";

    cout << "a + b = " << (a + b).to_float() << "\n";
    cout << "a - b = " << (a - b).to_float() << "\n";
    cout << "a * b = " << (a * b).to_float() << "\n";
    cout << "a / b = " << (a / b).to_float() << "\n";

    return 0;
}