// Copyright (C) 2017 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef DotProd_hh
#define DotProd_hh

//--------------------------------------------------------------------------//

/// This class defines a DotProd, whose task is to compute
/// the cross product between two vectors

template <class TYPE>
class DotProd {
public:

  /// Constructor
  DotProd() {}

  /// Desctructor
  ~DotProd() {}

  void computeDotProd(std::vector <TYPE> a, std::vector <TYPE> b)
  {
    c=0.0;
    if(a.size()!=b.size()) {
     std::cout << "DotProd::error in dimension size\n";
    }
    for (unsigned int i=0; i<a.size(); i++){
      c = c + a.at(i) * b.at(i);
    }
  }

  double getDotProd() { return c; }

private:

  /// vector of the cross product
  double c;
};

//--------------------------------------------------------------------------//

#endif // DotProd_hh 
