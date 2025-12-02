#include <iostream>
class Atom {
public:
  double x;
};

Atom atom;

void setvalue()
{
  atom.x=2.0;
}


//double Atom::x = 1.0;

int main()
{
  atom.x=1.0;
  std::cout << atom.x << std::endl;
  setvalue();
  std::cout << atom.x << std::endl;
}
