/*! \file TC_utils.c
 *  \brief Various utilities used by other functions
 */ 
/** 
 * \brief Much faster version of pow for small integers 
 *  (considering that reactions are generally no more than 3rd order).
 *  The C version of pow requires a double exponent; pow is the slowest part of the 
 *  code when used instead of fastIntPow. 
 */
static double fastIntPow(double val, int exponent)
{
  switch(exponent)
  {
  case 1:
    return ( val );
  case 2:
    return ( val*val );
  case 0:
    return ( 1.0 );
  case 3:
    return ( val*val*val );
  default:
    return ( pow(val,exponent) );
  }
}

