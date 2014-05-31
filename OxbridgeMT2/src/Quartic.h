// Source file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include <complex>

namespace Mt2 {
  namespace quartic {

    //******************************************************************************
    std::complex<double> c8_sqrt (const std::complex<double> x );
    double c8_argument (const std::complex<double> & x);
    double c8_magnitude (const std::complex<double> & x);
    double r8_sign ( double x );
    void r8poly3_root (const double a, const double b, const double c, const double d, std::complex<double> & r1, std::complex<double> & r2, std::complex<double> & r3 );

    void r8poly4_root (const double a,
		       const double b,
		       const double c,
		       const double d,
		       const double e,
		       std::complex<double> & r1,
		       std::complex<double> & r2,
		       std::complex<double> & r3,
		       std::complex<double> & r4 )

    //******************************************************************************
    //
    //  Purpose:
    //
    //    R8POLY4_ROOT returns the four roots of a quartic polynomial.
    //
    //  Discussion:
    //
    //    The polynomial has the form:
    //
    //      A * X**4 + B * X**3 + C * X**2 + D * X + E = 0
    //
    //  Modified:
    //
    //    27 October 2005
    //
    //  Parameters:
    //
    //    Input, double A, B, C, D, the coefficients of the polynomial.
    //    A must not be zero.
    //
    //    Output, complex *R1, *R2, *R3, *R4, the roots of the polynomial.
    //
    {
      double a3;
      double a4;
      double b3;
      double b4;
      double c3;
      double c4;
      double d3;
      double d4;
      std::complex<double> p;
      std::complex<double> q;
      std::complex<double> r;
      const std::complex<double> zero(0,0);

      if ( a == 0.0 )
	{
	  std::cout << "\n";
	  std::cout << "R8POLY4_ROOT - Fatal error!\n";
	  std::cout << "  A must not be zero!\n";
	  throw "R8POLY4_ROOT - Fatal error!  A must not be zero!";
	}

      a4 = b / a;
      b4 = c / a;
      c4 = d / a;
      d4 = e / a;
      //
      //  Set the coefficients of the resolvent cubic equation.
      //
      a3 = 1.0;
      b3 = -b4;
      c3 = a4 * c4 - 4.0 * d4;
      d3 = -a4 * a4 * d4 + 4.0 * b4 * d4 - c4 * c4;
      //
      //  Find the roots of the resolvent cubic.
      //
      r8poly3_root ( a3, b3, c3, d3, r1, r2, r3 );
      //
      //  Choose one root of the cubic, here R1.
      //
      //  Set R = sqrt ( 0.25 * A4**2 - B4 + R1 )
      //
  
      // if r1 is real, and if r1 + 0.25 * a4 * a4 - b4 is positive, then r that is about to follow will be exactly real.

      r = c8_sqrt ( r1 + 0.25 * a4 * a4 - b4 );
      //  if ( r != zero )
      if ((r.real()!=0.0) || (r.imag()!=0.0))
	{
	  p = c8_sqrt ( zero + 0.75 * a4 * a4 - r * r - 2.0 * b4 
			+ 0.25 * ( 4.0 * a4 * b4 - 8.0 * c4 - a4 * a4 * a4 ) / r );

	  q = c8_sqrt ( zero + 0.75 * a4 * a4 - r * r - 2.0 * b4 
			- 0.25 * ( 4.0 * a4 * b4 - 8.0 * c4 - a4 * a4 * a4 ) / r );
	}
      else
	{
	  p = c8_sqrt ( zero + 0.75 * a4 * a4 - 2.0 * b4 
			+ 2.0 * c8_sqrt ( r1 * r1 - 4.0 * d4 ) );

	  q = c8_sqrt ( zero + 0.75 * a4 * a4 - 2.0 * b4 
			- 2.0 * c8_sqrt ( r1 * r1 - 4.0 * d4 ) );
	}
      //
      //  Set the roots.
      //
      r1 = zero  -0.25 * a4 + 0.5 * r + 0.5 * p;
      r2 = zero  -0.25 * a4 + 0.5 * r - 0.5 * p;
      r3 = zero  -0.25 * a4 - 0.5 * r + 0.5 * q;
      r4 = zero  -0.25 * a4 - 0.5 * r - 0.5 * q;
      return;
    }


    //******************************************************************************

    void r8poly3_root (const double a, const double b, const double c, const double d, std::complex<double> & r1, std::complex<double> & r2, std::complex<double> & r3 )

    //******************************************************************************
    //
    //  Purpose:
    //
    //    R8POLY3_ROOT returns the three roots of a cubic polynomial.
    //
    //  Discussion:
    //
    //    The polynomial has the form
    //
    //      A * X**3 + B * X**2 + C * X + D = 0
    //
    //  Modified:
    //
    //    25 October 2005
    //
    //  Parameters:
    //
    //    Input, double A, B, C, D, the coefficients of the polynomial.
    //    A must not be zero.
    //
    //    Output, complex *R1, *R2, *R3, the roots of the polynomial, which
    //    will include at least one real root.
    //
    {
      std::complex<double> i;
      double pi = 3.141592653589793;
      double q;
      double r;
      double s1;
      double s2;
      double temp;
      double theta;
      const std::complex<double> zero(0,0);

      if ( a == 0.0 )
	{
	  std::cout << "\n";
	  std::cout << "R8POLY3_ROOT - Fatal error!\n";
	  std::cout << "  A must not be zero!\n";
	  throw "R8POLY3_ROOT - Fatal error! A must not be zero!";
	}

      i = std::complex<double> ( 0.0, 1.0 );

      q = ( pow ( b / a, 2 ) - 3.0 * ( c / a ) ) / 9.0;

      r = ( 2.0 * pow ( b / a, 3 ) - 9.0 * ( b / a ) * ( c / a ) 
	    + 27.0 * ( d / a ) ) / 54.0;

      if ( r * r < q * q * q )
	{
	  theta = acos ( r / sqrt ( pow ( q, 3 ) ) );
	  r1 = -2.0 * sqrt ( q ) * cos (   theta              / 3.0 );
	  r2 = -2.0 * sqrt ( q ) * cos ( ( theta + 2.0 * pi ) / 3.0 );
	  r3 = -2.0 * sqrt ( q ) * cos ( ( theta + 4.0 * pi ) / 3.0 );
	}
      else if ( q * q * q <= r * r )
	{
	  temp = -r + sqrt ( r * r - q * q * q );
	  s1 = r8_sign ( temp ) * pow ( fabs ( temp ), 1.0 / 3.0 );

	  temp = -r - sqrt ( r * r - q * q * q );
	  s2 = r8_sign ( temp ) * pow ( fabs ( temp ), 1.0 / 3.0 );

	  r1 = s1 + s2;
	  r2 = zero -0.5 * ( s1 + s2 ) + i * 0.5 * sqrt ( 3.0 ) * ( s1 - s2 );
	  r3 = zero -0.5 * ( s1 + s2 ) - i * 0.5 * sqrt ( 3.0 ) * ( s1 - s2 );
	}

      r1 = r1 - b / ( 3.0 * a );
      r2 = r2 - b / ( 3.0 * a );
      r3 = r3 - b / ( 3.0 * a );

      return;
    }


    //******************************************************************************

    std::complex<double> c8_sqrt (const std::complex<double> x )

      //******************************************************************************
      //
      //  Purpose:
      //
      //    C8_SQRT returns the principal square root of a C8.
      //
      //  Discussion:
      //
      //    A C8 is a double precision complex value.
      //
      //  Modified:
      //
      //    27 October 2005
      //
      //  Author:
      //
      //    John Burkardt
      //
      //  Parameters:
      //
      //    Input, complex X, the number whose square root is desired.
      //
      //    Output, complex C8_SQRT, the square root of X.
      //
      {
	double argument;
	double magnitude;
	std::complex<double> value;
  
	argument = c8_argument ( x );
	magnitude = c8_magnitude ( x );

	if ( magnitude == 0.0 )
	  {
	    value = std::complex<double> ( 0.0, 0.0 );
	  }
	else
	  {
	    value = pow ( magnitude, 1.0 / 2.0 ) 
	      * std::complex<double> ( cos ( argument / 2.0 ), sin ( argument / 2.0 ) );
	  }

	return value;
      }

    //*********************************************************************

    double r8_sign ( double x )

    //*********************************************************************
    //
    //  Purpose:
    //
    //    R8_SIGN returns the sign of a double precision number.
    //
    //  Modified:
    //
    //    18 October 2004
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Parameters:
    //
    //    Input, double X, the number whose sign is desired.
    //
    //    Output, double R8_SIGN, the sign of X.
    //
    {
      if ( x < 0.0 )
	{
	  return ( -1.0 );
	} 
      else
	{
	  return ( 1.0 );
	}
    }

    //******************************************************************************

    double c8_argument (const std::complex<double> & x)

    //******************************************************************************
    //
    //  Purpose:
    //
    //    C8_ARGUMENT returns the argument of a C8.
    //
    //  Discussion:
    //
    //    A C8 is a double precision std::complex<double> value.
    //
    //  Modified:
    //
    //    21 October 2005
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Parameters:
    //
    //    Input, complex X, the value whose argument is desired.
    //
    //    Output, double C8_ARGUMENT, the argument of X.
    //
    {
      double argument;

      argument = atan2(x.imag(),x.real());

      return argument;
    }

    //******************************************************************************

    double c8_magnitude (const std::complex<double> & x)

    //******************************************************************************
    //
    //  Purpose:
    //
    //    C8_MAGNITUDE returns the magnitude of a C8.
    //
    //  Discussion:
    //
    //    A C8 is a double precision complex value.
    //
    //  Modified:
    //
    //    22 October 2005
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Parameters:
    //
    //    Input, complex X, the value whose norm is desired.
    //
    //    Output, double C8_MAGNITUDE, the magnitude of X.
    //
    {
      return abs(x);
      //double magnitude;
      //
      //magnitude = sqrt(pow(x.real(),2) + pow(x.imag(),2));
      //
      //return magnitude;
    }


  } // namespace quartic
} //  namespace Mt2 

