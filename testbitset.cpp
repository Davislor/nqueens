#include <iostream>
#include <bitset>

using namespace std;

int main(void)
{
  cout << "A bitfield<" << 1 << "> is " << sizeof(bitset<1>);
  cout << " bytes long." << endl;

  cout << "A bitfield<" << 2 << "> is " << sizeof(bitset<2>);
  cout << " bytes long." << endl;

  cout << "A bitfield<" << 3 << "> is " << sizeof(bitset<3>);
  cout << " bytes long." << endl;

  cout << "A bitfield<" << 4 << "> is " << sizeof(bitset<4>);
  cout << " bytes long." << endl;

  cout << "A bitfield<" << 8 << "> is " << sizeof(bitset<8>);
  cout << " bytes long." << endl;

  cout << "A bitfield<" << 12 << "> is " << sizeof(bitset<12>);
  cout << " bytes long." << endl;

  cout << "A bitfield<" << 16 << "> is " << sizeof(bitset<16>);
  cout << " bytes long." << endl;

  cout << "A bitfield<" << 24 << "> is " << sizeof(bitset<24>);
  cout << " bytes long." << endl;

  cout << "A bitfield<" << 32 << "> is " << sizeof(bitset<32>);
  cout << " bytes long." << endl;

  cout << "A bitfield<" << 40 << "> is " << sizeof(bitset<40>);
  cout << " bytes long." << endl;

  cout << "A bitfield<" << 48 << "> is " << sizeof(bitset<48>);
  cout << " bytes long." << endl;

  cout << "A bitfield<" << 56 << "> is " << sizeof(bitset<56>);
  cout << " bytes long." << endl;

  cout << "A bitfield<" << 64 << "> is " << sizeof(bitset<64>);
  cout << " bytes long." << endl;

  return 0;
}
