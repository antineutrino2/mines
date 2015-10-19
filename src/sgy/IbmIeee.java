package sgy;

public class IbmIeee {
  public static void shortToFloat(short[] s, float[] f) {
    int n = s.length;
    for (int i=0; i<n; ++i)
      f[i] = s[i];
  }
  public static void ieeeToFloat(int[] ieee, float[] f) {
    int n = ieee.length;
    for (int i=0; i<n; ++i)
      f[i] = ieeeToFloat(ieee[i]);
  }
  public static void ibmToFloat(int[] ibm, float[] f) {
    int n = ibm.length;
    for (int i=0; i<n; ++i)
      f[i] = ibmToFloat(ibm[i]);
  }
  public static float ieeeToFloat(int ieee) {
    return Float.intBitsToFloat(ieee);
  }
  public static int floatToIeee(float f) {
    return Float.floatToIntBits(f);
  }
  public static void floatToIeee(int[] ieee, float[] f) {
	    int n = f.length;
	    for (int i=0; i<n; ++i)
	      ieee[i] = floatToIeee(f[i]);
  }
//  public static void floatToIBM(int[] ibm, float[] f) {
//		int n = f.length;
//	    for (int i=0; i<n; ++i)
//	    	ibm[i] = floatToIBM(f[i]);
//  }
  public static float ibmToFloat(int ibm) {
    // 1) Extract sign bit, exponent, and mantissa.
    // 2) Convert exponent: subtract 64, multiply by 4, subtract 1, add 127
    // 3) Convert mantissa:
    //    while high mantissa bit is zero {
    //      shift mantissa left (low order bits are zeroed)
    //      decrement exponent
    // 4) Put sign and exponent bits back in number
    // 5) Reverse order of bytes?
    /*
    ibm = ((ibm    )     )<<24 |
          ((ibm>> 8)&0xff)<<16 |
          ((ibm>>16)&0xff)<< 8 |
          ((ibm>>24)&0xff);
    */
    int s = 0x80000000&ibm; // sign bit
    int e = 0x7f000000&ibm; // exponent
    int m = 0x00ffffff&ibm; // mantissa
    int ieee = 0;
    if (m!=0) {
      e = (e>>22)-130; // = ((e>>24)-64)*4 - 1 + 127
      while ((m&0x00800000)==0) {
        m <<= 1;
        --e;
      }
      if (e<=0) {
        ieee = 0;
      } else if (e>=255) {
        ieee = s|0x7f7fffff;
      } else {
        ieee = s|(e<<23)|(m&0x007fffff);
      }
    }
    /*
    ieee = ((ieee    )     )<<24 |
           ((ieee>> 8)&0xff)<<16 |
           ((ieee>>16)&0xff)<< 8 |
           ((ieee>>24)&0xff);
    */
    return Float.intBitsToFloat(ieee);
  }
//  public static int floatToIBM(float f) {
//	    // 1) Extract sign bit, exponent, and mantissa.
//	    // 2) Convert exponent: subtract 64, multiply by 4, subtract 1, add 127
//	    // 3) Convert mantissa:
//	    //    while high mantissa bit is zero {
//	    //      shift mantissa left (low order bits are zeroed)
//	    //      decrement exponent
//	    // 4) Put sign and exponent bits back in number
//	    // 5) Reverse order of bytes?
//	    /*
//	    ibm = ((ibm    )     )<<24 |
//	          ((ibm>> 8)&0xff)<<16 |
//	          ((ibm>>16)&0xff)<< 8 |
//	          ((ibm>>24)&0xff);
//	    */
//	  	int ibm = 0;
//	  	int ieee = Float.floatToIntBits(f);
//	    int s = 0x80000000; // sign bit
//	    int e = 0x7f800000; // exponent
//	    int m = 0x7fffff;   // mantissa
//	    int signbit = ieee&s;
//	    int exponent = ((ieee&e) >> 23) - 127;
//	    // The IBM 7-bit exponent is to the base 16 and the mantissa is 
//	    // presumed to be entirely to the right of the radix point. In 
//	    // contrast, the IEEE exponent is to the base 2 and there is an 
//	    // assumed 1-bit to the left of the radix point.
//	    double exp16 = Math.floor((exponent+1) / 4);
//	    int exp_remainder = (exponent+1) % 4;
//	    int downshift = 0;
//	    if (exp_remainder != 0) {
//	    	exp16 += exp_remainder;
//	    	downshift = 4-exp_remainder;
//	    }
//	    int ibmexp = exp16 + 64;
//	    return ibm;
//	  }
  private static float ibmToFloatSu(int ibm) {
    int fconv = ((ibm    )     )<<24 |
                ((ibm>> 8)&0xff)<<16 |
                ((ibm>>16)&0xff)<< 8 |
                ((ibm>>24)&0xff);
    int fmant = 0x00ffffff&fconv;
    if (fconv!=0 && fmant!=0) {
      int t = ((0x7f000000&fconv)>>>22)-130;
      while ((fmant&0x00800000)==0) {
        --t;
        fmant <<= 1;
      }
      if (t>254) {
        fconv = (0x80000000&fconv)|0x7f7fffff;
      } else if (t<=0) {
        fconv = 0;
      } else {
        fconv = (0x80000000&fconv)|(t<<23)|(0x007fffff&fmant);
      }
    }
    return Float.intBitsToFloat(fconv);
  }
}
