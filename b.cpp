#include <iostream>

using namespace std ;

int main() {
    int a[3] = {1, 2, 3} ;
    int i = 0 ;
    int b = a[i++] ;
    cout << b << endl ;
    cout << i << endl ;
    return  0 ;
}
