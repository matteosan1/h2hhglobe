#ifndef RESOLFUNCTSERVICE_H
#define RESOLFUNCTSERVICE_H

/// Service to build the scale functor corresponding to the passed identifier                                                                               
resolFunctBase<double * > * resolutionFunctService( const int identifier ){
  switch ( identifier ) {
  case ( 45 ): return ( new resolFunct45<double * > ); break;
  default: std::cout << "resolutionFunctService error: wrong identifier = " << identifier << std::endl; exit(1);
  }
}

#endif
