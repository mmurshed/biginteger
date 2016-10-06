/**
 * BigInteger C-Code
 * Version 6.6.25
 *
 * Copyright (c) 2004
 * Mahbub Murshed Suman (suman@bttb.net.bd)
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation. Mahbub Murshed Suman makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 */

struct BigInteger
{
    // The integer array to hold the number
    int *TheNumber;
    // Start of the location of the number in the array
    unsigned int Start;
    // End  of the location of the number in the array
    unsigned int End;
    // True if the number is negative
    bool isNegative;
};


  // Character array constructor
  void BigInteger(BigInteger *x, char const* n)
  {
    if(n[0]=='-') { x->isNegative = true; n++; }
    else if(n[0]=='+') { x->isNegative = false; n++; }
    else x->isNegative = false;

    while(*n=='0') n++;

    int l = strlen(n);
    if(l==0)
    {
      x->TheNumber = 0;
	  x->Start = x->End = 0;
      return;
    }
    x->Start = 0;
    x->End = (unsigned int)(l/LOG10BASE + l%LOG10BASE - 1);
    x->TheNumber = new int [Digits(x)];

    int cur = l - 1;
    for(unsigned int i = x->End; i>=x->Start;i--)
    {
      if(cur<0) break;
        int r = 0;
      int TEN = 1;
      for(l=0;l<4;l++)
      {
        if(cur<0) break;
        r = r + TEN*(n[cur]-'0');
        TEN *= 10;
        cur--;
      }
      TheNumber[i] = r;
    }
    TrimZeros(x);
    if(isZero()) x->isNegative = false;
  }

  // Copy constructor
  CopyBigInteger(BigInteger *from, BigInteger* to)
  {
    to->Start = 0;
    to->End = from->Digits() - 1;
    to->TheNumber = new int [from.Digits()];

    datacopy(from,to, from.Digits());
    to->isNegative = from.isNegative;
  }

  // Copies data from `a' to `b'
  void datacopy(BigInteger* a,BigInteger* b, SizeT size)
  {
    for(SizeT i=0;i<size;i++)
      b->TheNumber[Start+i] = a->TheNumber[a->Start+i];
  }

  // Returns number of digits in a BigInteger
  unsigned int Digits(BigInteger *x) const
  { return x->End - x->Start+1; }

  // true if 'this' is zero
  bool isZero(BigInteger *x) const
  { return x->End==x->Start && x->TheNumber[x->Start]==0; }


  // Trims Leading Zeros
  void TrimZeros(BigInteger *x)
  {
    while(x->TheNumber[x->Start]==0 && x->Start<x->End)
      x->Start++;
  }

  // Compares this with `with' irrespective of sign
  // Returns
  // 0 if equal
  // 1 if this>with
  // -1 if this<with
  int UnsignedCompareTo(BigInteger *x, BigInteger* with)
  {
    if(isZero(x) && isZero(with)) return 0;
    if(!isZero(x) && isZero(with)) return 1;
    if(isZero(x) && !isZero(with)) return -1;

    long temp = Digits(x) - Digits(with);
    // Case 3: First One got more digits
    // Case 4: First One got less digits
    if(temp!=0) return temp<0?-1:1;

    // Now I know Both have same number of digits
    // Case 5,6,7:
    /*
    Compares two array of data considering the
    case that both of them have the same number
    of digits
    */

    temp = 0;
    for(unsigned int i=0;i<Digits(x);i++)
    {
      temp = x->TheNumber[i+x->Start] - with->TheNumber[i+with->Start];
      if(temp!=0)
        return temp<0?-1:1;
    }

    return 0;
  }

  // Compares this with `with'
  // Returns
  // 0 if equal
  // 1 if this>with
  // -1 if this<with
  int CompareTo(BigInteger* x, BigInteger* with)
  {
    int cmp = UnsignedCompareTo(x,with);

    // Case 1: Positive , Negative
    if(x->isNegative==false && with->isNegative==true) return 1;
    // Case 2: Negative, Positive
    if(x->isNegative==true && with->isNegative==false) return -1;
    // Now, Both have Same Sign
    int both_neg = 1;
    if(x->isNegative==true && with->isNegative==true) both_neg = -1;
    return cmp*both_neg;
  }

  // Implentation of addition by paper-pencil method
void Add(BigInteger *Big, BigInteger* Small,BigInteger* Result) const
  {
    Result->TheNumber = *new int[Digits(Big)+1];

    long Carry=0,Plus;
    long i=Digits(Big) - 1,
      j= Digits(Small) - 1;

    for(; i>=0 ;i--,j--){
      Plus = Big->TheNumber[i+Big->Start] + Carry;
      if(j>=0) Plus += Small->TheNumber[j+Small->Start] ;

      Result->TheNumber[i+1] = Plus%BASE;
      Carry = Plus/BASE;
    }
    i++;

    if(Carry) Result->TheNumber[i--] = Carry;

    TrimZeros(Result);
  }

  // Addition operator
/*
  operator+(BigInteger const& N1, BigInteger const& N2)
  {
    if(N1.isNegative && !N2.isNegative)
    {
      // First is negaive and second is not
      BigInteger& Res = *new BigInteger(N1);
      Res.isNegative = false;
      Res = N2-Res;
      return Res;
    }
    if(!N1.isNegative && N2.isNegative)
    {
      // Second is negaive and first is not
      BigInteger& Res = *new BigInteger(N2);
      Res.isNegative = false;
      Res = N1-Res;
      return Res;
    }

    BigInteger& ret = *new BigInteger();

    if(N1.UnsignedCompareTo(N2)<0)
      ret = N2.Add(N1);
    else
      ret = N1.Add(N2);
    if(N1.isNegative==true && N2.isNegative==true)
      ret.isNegative = true;
    return ret;
  }
*/
  // Implentation of subtraction by paper-pencil method
  void Subtract(BigInteger* Big, BigInteger* Small,BigInteger* Result)const
  {
    Result->TheNumber = *new int[Digits(Big)+1];

    long Carry=0,Minus;

    long i = Big->Digits() - 1,
      j= Small->Digits() - 1;

    for( ; i>=0 ;i--,j--)
    {
      Minus = Big->TheNumber[i+Big->Start] - Carry;
      if(j>=0) Minus -= Small->TheNumber[j+Small->Start];

      if(Minus < 0)
      {
        Result->TheNumber[i+1] = Minus + BASE;
        Carry = 1;
      }
      else
      {
        Result->TheNumber[i+1]  = Minus;
        Carry = 0;
      }
    }
    TrimZeros(Result);
  }

/*
  // Subtruction operator
  BigInteger& operator-(BigInteger const& N1, BigInteger const& N2)
  {
    if(N1.isNegative && !N2.isNegative)
    {
      BigInteger& res = *new BigInteger(N1);
      res.isNegative = false;
      res = res+N2;
      res.isNegative = true;
      return res;
    }
    if(!N1.isNegative && N2.isNegative)
    {
      BigInteger& res = *new BigInteger(N2);
      res.isNegative = false;
      res = res+N1;
      return res;
    }


    BigInteger& ret = *new BigInteger();
    int cmp = N1.UnsignedCompareTo(N2);
    if(cmp==0)
    {
      ret = *new BigInteger();
    }
    if(cmp<0)
    {
      ret = N2.Subtract(N1);
      ret.isNegative = true;
    }
    else
    {
      ret = N1.Subtract(N2);
      ret.isNegative = false;
    }
    return ret;
  }
*/
  // Implentation of multiplication by paper-pencil method
  BigInteger& BigInteger::Multiply(BigInteger const& Small) const
  {
    BigInteger const& Big = *this;
    BigInteger& Result = *new BigInteger(Big.Digits()+Small.Digits(),0);

    long Carry,Multiply;

    SizeT i;
    SizeT j;

    for(i = 0 ; i< Small.Digits() ; i++)
    {
      Carry = 0;
      for(j = 0 ; j< Big.Digits() ; j++)
      {
        Multiply = ( (long)Small.TheNumber[Small.End-i] * (long)Big.TheNumber[Big.End-j] )
          + Carry + Result.TheNumber[Result.End-i-j];
        Result.TheNumber[Result.End-i-j] = Multiply%BASE;
        Carry = Multiply/BASE ;
      }
      Result.TheNumber[Result.End-i-j] = Carry;
    }

    Result.TrimZeros();

    return Result;
  }

  // Implentation of multiplication by paper-pencil method
  // See: D. E. Knuth 4.3.1
  BigInteger& BigInteger::Multiply(DATATYPE const& with) const
  {
    BigInteger& Result = *new BigInteger(Digits()+1,0);

    long Carry,Multiply;

    SizeT i;

    Carry = 0;
    for(i = 0 ; i< Digits() ; i++)
    {
      Multiply = Carry + (long)TheNumber[End-i] * (long)with;
      Carry = Multiply / BASE;
      Result.TheNumber[Result.End-i] = Multiply % BASE;
    }
    Result.TheNumber[Result.End-i] = Carry;
    Result.TrimZeros();

    return Result;
  }

  // Multiplication operator
  BigInteger& operator*(BigInteger const& N1, BigInteger const& N2)
  {
    if(N1.isZero() || N2.isZero()) return *new BigInteger();
    BigInteger& ret = *new BigInteger();
    if(N1.UnsignedCompareTo(N2)<0)
      ret = N2.Multiply(N1);
    else
      ret = N1.Multiply(N2);
    if(N1.isNegative!=N2.isNegative)
      ret. isNegative = true;
    return ret;
  }

  // Multiplication operator
  BigInteger& operator*(BigInteger const& N1, DATATYPE const& N2)
  {
    if(N1.isZero() || N2==0) return *new BigInteger();
    BigInteger& ret = N1.Multiply(N2);
    // if(N1.isNegative!=(N2<0)) ret. isNegative = true;
    return ret;
  }

  // Multiplication operator
  BigInteger& operator*(DATATYPE const& N1, BigInteger const& N2)
  {
    if(N2.isZero() || N1==0) return *new BigInteger();
    BigInteger& ret = N2.Multiply(N1);
    // if(N2.isNegative!=(N1<0)) ret. isNegative = true;
    return ret;
  }

  // Implentation of division by paper-pencil method
  // Here `this' is divided by 'V'
  BigInteger& BigInteger::DivideAndRemainder(DATATYPE const& V,DATATYPE& _R,bool skipRemainder=false) const
  {
    BigInteger& W = *new BigInteger(Digits(),0,false);
    DATATYPE R = 0;
    SizeT j;
    for(j=0;j<=W.End;j++)
    {
      W.TheNumber[j] = (R*BASE+TheNumber[Start+j])/V;
      R = (R*BASE+TheNumber[Start+j])%V;
    }

    W.TrimZeros();
    if(skipRemainder==false)
      _R = R;

    return W;
  }

  // Implentation of division by paper-pencil method
  // Does not perform the validity and sign check
  // It is assumed that this > `_V'
  // See: D.E.Knuth 4.3.1
  BigInteger& BigInteger::DivideAndRemainder(BigInteger const& _V,BigInteger& R,bool skipRemainder=false) const
  {
    SizeT m = this->Digits()-_V.Digits();
    SizeT n = _V.Digits();
    BigInteger& Q = *new BigInteger(m+1,0,false);

    DATATYPE d, qhat, rhat;
    long temp,x,y;
    SizeT i;
    int j;

    d = (BASE-1)/_V.TheNumber[_V.Start];

    BigInteger& U = this->Multiply(d);
    BigInteger& V = _V.Multiply(d);

    for(j = m; j>=0; j--)
    {
      temp = (long)U.TheNumber[U.End-j-n]*(long)BASE + (long)U.TheNumber[U.End-j-n+1];
      x = temp / (long)V.TheNumber[V.Start];
      y = temp % (long)V.TheNumber[V.Start];
      if(x>(long)BASE) x /= BASE;
      if(y>(long)BASE) y %= BASE;
      qhat = (DATATYPE) x;
      rhat = (DATATYPE) y;

	  bool badRhat = false;
	  do
	  {
		x = (long)qhat * (long)V.TheNumber[V.Start+1];
		y = (long)BASE*(long)rhat + (long)U.TheNumber[U.End-j-n+1];

		if(qhat==BASE || x > y)
		{
		  qhat --;
		  rhat += V.TheNumber[V.Start];
		  if(rhat<BASE) badRhat = true;
		  else badRhat = false;
		}
		else break;
	  }while(badRhat);

      // `x' Will act as Carry in the next loop
      x = 0;
      temp = 0;
      for(i=0;i<=n;i++)
      {
        if(V.End>=i) temp = (long)qhat*(long)V.TheNumber[V.End-i] + temp;
        y = (long)U.TheNumber[U.End-j-i] - temp%BASE - x;
        temp /= BASE;
        if(y < 0)
        {
          U.TheNumber[U.End-j-i] = (DATATYPE)(y+BASE);
          x = 1;
        }
        else
        {
          U.TheNumber[U.End-j-i]  = (DATATYPE)y;
          x = 0;
        }
      }
      // if(x) U.TheNumber[U.Start+j+i] --;

      Q.TheNumber[Q.End-j] = qhat;

      if(x)
      {
        Q.TheNumber[Q.End-j] --;
        // `x' Will act as Carry in the next loop
        x = 0;
        for(i=0;i<=n;i++)
        {
          y = (long)U.TheNumber[U.End-j-i] + x;
          if(V.End>=i) y += (long)V.TheNumber[V.End-i];
          U.TheNumber[U.End-j-i] = (DATATYPE)(y % BASE);
          x = y / BASE;
        }
        U.TheNumber[U.End-j-i] = (DATATYPE)x;
      }
    }

    U.TrimZeros();
    DATATYPE _t;
    if(skipRemainder==false)
      R = U.DivideAndRemainder(d,_t,true);
    Q.TrimZeros();

    return Q;
  }

  // Front end for Divide and Remainder
  // Performs the validity and sign check
  BigInteger& DivideAndRemainder(BigInteger const& U, BigInteger const& V,BigInteger& R,bool skipRemainder=false)
  {
    if(V.isZero())
      throw logic_error (DumpString ("DivideAndRemainder (BigInteger/BigInteger)",BigMathDIVIDEBYZERO));

    if(U.isZero())
    {
      R = *new BigInteger();
      return *new BigInteger();
    }

    int cmp = U.UnsignedCompareTo(V);
    if(cmp==0)
    {
      R = *new BigInteger();
      BigInteger& W = *new BigInteger(1l);
      if(U.isNegative!=V.isNegative) W.isNegative = true;
      return W;
    }
    else if(cmp<0)
    {
      R = *new BigInteger(U);
      if(U.isNegative!=V.isNegative) R.isNegative = true;
      return *new BigInteger();
    }
    BigInteger& ret = U.DivideAndRemainder(V,R,skipRemainder);
    if(U.isNegative!=V.isNegative) ret.isNegative = true;
    return ret;
  }



void main()
{
	return 0;
}
