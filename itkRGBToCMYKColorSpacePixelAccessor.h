/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRGBToCMYKColorSpacePixelAccessor.h,v $
  Language:  C++
  Date:      $Date: 2009-03-03 15:08:46 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRGBToCMYKColorSpacePixelAccessor_h
#define __itkRGBToCMYKColorSpacePixelAccessor_h


#include "itkRGBPixel.h"
#include "itkVector.h"
#include "vnl/vnl_math.h"
#include "itkNumericTraits.h"

namespace itk
{
namespace Accessor
{
/**
 * \class RGBToCMYKColorSpacePixelAccessor
 * \brief Give access to a RGBPixel as if it were in CMYK Color Space as a Vector type.
 *
 * This class is intended to be used as parameter of 
 * an ImageAdaptor to make an RGBPixel image appear as being
 * an image of Vector pixel type in CMYK Color Space.
 *
 * \sa ImageAdaptor
 * \ingroup ImageAdaptors
 *
 */

template <class TInput, class TOutput>
class ITK_EXPORT RGBToCMYKColorSpacePixelAccessor
{
public:
  /** Standard class typedefs. */
  typedef   RGBToCMYKColorSpacePixelAccessor        Self;

 /** External typedef. It defines the external aspect
   * that this class will exhibit */
  typedef  Vector<TOutput,4>     ExternalType;

  /** Internal typedef. It defines the internal real
   * representation of data */
  typedef   RGBPixel<TInput>    InternalType;

  /** Write access to the RGBToCMYKColorSpace component */
  inline void Set( InternalType & output, const ExternalType & input ) const
    { 
    // Normalize RGB values.
    double r = (double)input[0] / (double)NumericTraits<TInput>::max();
    double g = (double)input[1] / (double)NumericTraits<TInput>::max();
    double b = (double)input[2] / (double)NumericTraits<TInput>::max();
    
    double C = 1 - r;
    double M = 1 - g;
    double Y = 1 - b;
    
    double k = 1;
    if(C < k)
      k = C;
    if(M < k)
      k = M;
    if(Y < k)
      k = Y;
      
    if(k == 1)
      {
	  C = 0;
	  M = 0;
	  Y = 0;    
      }
    else
      {
	  C = (C - k) / (1 - k);    
	  M = (M - k) / (1 - k);
	  Y = (Y - k) / (1 - k);
	  }

    output[0] = static_cast<TInput>(C); // C
    output[1] = static_cast<TInput>(M); // M
    output[2] = static_cast<TInput>(Y); // Y
    output[3] = static_cast<TInput>(k); // K
    
    return output;
    }

  /** Read access to the RGBToCMYKColorSpace component */
  inline ExternalType Get( const InternalType & input ) const
    {
    // Normalize RGB values.
    double r = (double)input[0] / (double)NumericTraits<TInput>::max();
    double g = (double)input[1] / (double)NumericTraits<TInput>::max();
    double b = (double)input[2] / (double)NumericTraits<TInput>::max();
    
    double C = 1 - r;
    double M = 1 - g;
    double Y = 1 - b;
    
    double k = 1;
    if(C < k)
      k = C;
    if(M < k)
      k = M;
    if(Y < k)
      k = Y;
      
    if(k == 1)
      {
	  C = 0;
	  M = 0;
	  Y = 0;    
      }
    else
      {
	  C = (C - k) / (1 - k);    
	  M = (M - k) / (1 - k);
	  Y = (Y - k) / (1 - k);
	  }
      
    ExternalType output;
    output[0] = static_cast<TOutput>(C); // C
    output[1] = static_cast<TOutput>(M); // M
    output[2] = static_cast<TOutput>(Y); // Y
    output[3] = static_cast<TOutput>(k); // K
    
    return output;
    }

private:
};
  
}  // end namespace Accessor
}  // end namespace itk

#endif
