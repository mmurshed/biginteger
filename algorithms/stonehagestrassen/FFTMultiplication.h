#ifndef FFT_MULTIPLICATION
#define FFT_MULTIPLICATION

#include <vector>
#include <complex>
#include <cmath>

#include "../BigIntegerUtil.h"

using namespace std;

namespace BigMath
{
    const long double PI = acos(-1.0L);

    class FFTMultiplication
    {
        // In-place iterative Cooley-Tukey FFT.
        // If invert is true, computes the inverse FFT.
        static void fft(vector<complex<long double>> &a, bool invert)
        {
            Int n = (Int)a.size();
            for (Int i = 1, j = 0; i < n; ++i)
            {
                Int bit = n >> 1;
                while (j & bit)
                {
                    j ^= bit;
                    bit >>= 1;
                }
                j ^= bit;
                if (i < j)
                    swap(a[i], a[j]);
            }

            for (Int len = 2; len <= n; len <<= 1)
            {
                long double angle = 2 * PI / len * (invert ? -1 : 1);
                complex<long double> wlen(cos(angle), sin(angle));
                for (Int i = 0; i < n; i += len)
                {
                    complex<long double> w(1);
                    for (int j = 0; j < len / 2; j++)
                    {
                        complex<long double> u = a[i + j];
                        complex<long double> v = a[i + j + len / 2] * w;
                        a[i + j] = u + v;
                        a[i + j + len / 2] = u - v;
                        w *= wlen;
                    }
                }
            }
            if (invert)
            {
                for (Int i = 0; i < n; i++)
                    a[i] /= n;
            }
        }

        // Computes convolution (polynomial multiplication) of two coefficient vectors A and B.
        // Returns a vector C such that C[k] = sum_{i+j=k} A[i]*B[j].
        static vector<DataT> convolution(const vector<DataT> &A, const vector<DataT> &B)
        {
            DataT n = 1;
            while (n < (Int)A.size() + (Int)B.size() - 1)
                n <<= 1;

            vector<complex<long double>> fa(n), fb(n);
            for (SizeT i = 0; i < A.size(); i++)
                fa[i] = A[i];
            for (SizeT i = A.size(); i < (size_t)n; i++)
                fa[i] = 0;
            for (SizeT i = 0; i < B.size(); i++)
                fb[i] = B[i];
            for (SizeT i = B.size(); i < (size_t)n; i++)
                fb[i] = 0;

            fft(fa, false);
            fft(fb, false);
            for (Int i = 0; i < n; i++)
                fa[i] *= fb[i];
            fft(fa, true);

            vector<DataT> C(n, 0);
            for (Int i = 0; i < n; i++)
            {
                // Round the real part to the nearest integer
                C[i] = (DataT)round(fa[i].real());
            }
            return C;
        }

    public:
        // Multiply two vectors of digits using FFT-based convolution.
        // Here each vector element is assumed to be a digit (0–9).
        static vector<DataT> Multiply(const vector<DataT> &a, const vector<DataT> &b, BaseT base)
        {
            // Compute the convolution which gives the raw coefficient vector.
            vector<DataT> c = convolution(a, b);
            DataT carry = 0;
            for (SizeT i = 0; i < c.size(); i++)
            {
                DataT total = c[i] + carry;
                c[i] = total % base;  // keep the current digit in base B
                carry = total / base; // propagate the carry
            }
            // Append any remaining carry.
            while (carry > 0)
            {
                c.push_back(carry % base);
                carry /= base;
            }
            BigIntegerUtil::TrimZeros(c);
            return c;
        }

        // Multiply two vectors of digits using FFT-based convolution.
        // Here each vector element is assumed to be a digit (0–9).
        static vector<DataT> Multiply(const vector<DataT> &a, DataT b, BaseT base)
        {
            if (b == 0 || BigIntegerUtil::IsZero(a))
            {
                return vector<DataT>(1, 0);
            }

            vector<DataT> c(1, b);
            return Multiply(a, c, base);
        }

        // Multiply two vectors of digits using FFT-based convolution.
        // Here each vector element is assumed to be a digit (0–9).
        static void MultiplyTo(vector<DataT> &a, DataT b, BaseT base)
        {
            if (b == 0 || BigIntegerUtil::IsZero(a))
            {
                a = vector<DataT>(1, 0);
                return;
            }
            vector<DataT> c(1, b);
            a = Multiply(a, c, base);
        }
    };
}
#endif // TOOMCOOK_MULTIPLICATION_2a