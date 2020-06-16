//Note: the original version of this code comes from https://github.com/antiduh/ErrorCorrection.git

using System;

namespace ReedSolomonOnBytes
{
    public readonly struct GaloisField256
    { 
	    public readonly int Size; // Fixed at 2^8 since we want to operate only on bytes

        private readonly byte[] m_logarithms;

        private readonly byte[,] m_multTable;

        public GaloisField256(int fieldGenPoly, int size = 256)
        {
	        Size = size;
	        m_multTable = default;
            Field = default;
            Inverses = default;
            m_logarithms = default;

	        Field = BuildField(fieldGenPoly);
            m_logarithms = BuildLogarithms();
            m_multTable = BuildMultTable();
            Inverses = BuildInverses();
        }

        public byte[] Field { get; }

        public byte[] Inverses { get; }

        public byte Multiply(byte left, byte right)
        {
            return m_multTable[left, right];
        }

        public byte Divide(byte dividend, byte divisor)
        {
            return m_multTable[dividend, Inverses[divisor]];
        }

        public Span<byte> PolyMult(in ReadOnlySpan<byte> left, in ReadOnlySpan<byte> right)
        {
	        Span<byte> result = new byte[left.Length + right.Length - 1];

            for(var leftIndex = 0; leftIndex < left.Length; ++leftIndex)
            {
                for(var rightIndex = 0; rightIndex < right.Length; ++rightIndex)
                {
	                result[leftIndex + rightIndex] = (byte)(result[leftIndex + rightIndex] ^ Multiply(left[leftIndex], right[rightIndex]));
                }
            }

            return result;
        }

        public byte PolyEval(in ReadOnlySpan<byte> poly, byte x)
        {
	        var sum = poly[0];

            for(var i = 1; i < poly.Length; ++i)
            {
	            if (poly[i] == 0)
	            {
		            continue;
	            }

                var power = (m_logarithms[poly[i]] + m_logarithms[x] * i) % (Size - 1);
                
                sum ^= Field[power + 1];
            }

            return sum;
        }

		private byte[] BuildField(int fieldGenPoly)
		{
			var field = new byte[Size];
	        
            field[1] = 1;

            var f = 1;

            for (var i = 2; i < Size; ++i)
            {
                f <<= 1;

                if (f >= Size)
                {
                    f ^= fieldGenPoly;
                }

                field[i] = (byte)f;
            }

            return field;
		}

        private byte[] BuildLogarithms()
        {
	        var logarithms = new byte[Size];

            for(var i = 0; i < Size; ++i)
            {
	            logarithms[Field[i]] = (byte)(i - 1);
            }

            return logarithms;
        }

        private byte[,] BuildMultTable()
        {
            var table = new byte[Size, Size];

            for(var left = 0; left < Size; ++left)
            {
                for(var right = 0; right < Size; ++right)
                {
	                table[left, right] = InternalMult(left, right);
                }
            }

            return table;
        }

        private byte[] BuildInverses()
        {
	        var inverses = new byte[Size];
	        
            for(var i = 1; i < inverses.Length; ++i)
            {
	            inverses[Field[i]] = InternalDivide(1, Field[i]);
            }

            return inverses;
        }

        private byte InternalMult(int left, int right)
        {
	        if (left == 0 || right == 0)
	        {
		        return 0;
	        }

            return Field[((m_logarithms[left] + m_logarithms[right]) % (Size - 1)) + 1];
        }

        private byte InternalDivide(int dividend, int divisor)
        {
	        if (dividend == 0)
	        {
		        return 0;
	        }

            return Field[((m_logarithms[dividend] - m_logarithms[divisor] + (Size - 1)) % (Size - 1)) + 1 ];
        }
    }
}
