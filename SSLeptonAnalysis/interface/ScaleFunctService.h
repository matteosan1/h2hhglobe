#ifndef SCALEFUNCTIONSERVICE_H
#define SCALEFUNCTIONSERVICE_H

#include "SSLeptonAnalysis/interface/ScaleFunct.h"

scaleFunctBase<double * > * scaleFunctService( const int identifier ){
  switch ( identifier ) {
  case ( 50 ): return ( new scaleFunct50<double * > ); break;
  default: std::cout << "scaleFunctService error: wrong identifier = " << identifier << std::endl; exit(1);
  }
}

#endif
