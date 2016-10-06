#include <iostream>

using namespace std;

#ifdef _MSC_VER
    typedef _int64 Long;
    typedef _int32 Int;
#elif __GNUC__
    typedef long long Long;
    typedef long Int;
#elif __BCPLUSPLUS__
    typedef __int64 Long;
    typedef __int32 Int;
#else
    typedef long Long;
    typedef int Int;
#endif

int main()
{
    cout << "Compiler ";
    #ifdef _MSC_VER
        cout << "MSVC (" << _MSC_VER << ")";
    #elif __GNUC__
        cout << "GCC (" << __GNUC__ << ")";
    #elif __BCPLUSPLUS__
        cout << "BCC (" << __BCPLUSPLUS__ << ")";
    #else
        cout << "Unknown";
    #endif
    cout << endl;

    Int a;
    Long b;
    cout << "Int : " << sizeof(a)*8 << " bit" << endl;
    cout << "Long : " << sizeof(b)*8 << " bit" << endl;

    return 0;
}
