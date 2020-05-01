//Note: the original version of this code comes from https://github.com/antiduh/ErrorCorrection.git

using System;

namespace ReedSolomonOnBytes
{
    public readonly ref struct Rs256EncoderLittleEndian // Little-Endian because byte message[0] is the highest mathematical message element
    {
	    private readonly int m_parityLength;
        
	    private readonly GaloisField256 m_gf;
        private readonly ReadOnlySpan<byte> m_codeGenPoly;

		public Rs256EncoderLittleEndian(int codewordLength, int messageLength, in GaloisField256 gf) 
        {
            m_parityLength = codewordLength - messageLength;

            m_gf = gf;

            m_codeGenPoly = default;
            m_codeGenPoly = BuildCodeGenPoly();
        }

        public void CalculateParity(in ReadOnlySpan<byte> message, Span<byte> parityOut)
        {
            parityOut.Fill(0);

	        byte r;
	        byte z_0 = 0;
	        
            for(var i = 0; i < message.Length - 1; ++i)
            {
                r = (byte)(z_0 ^ message[i]);

                for (var zIter = 0; zIter < parityOut.Length; ++zIter)
                {
	                parityOut[parityOut.Length - zIter - 1] ^= m_gf.Multiply(m_codeGenPoly[zIter], r);
                }

                z_0 = parityOut[0];

                for (var zIter = 0; zIter < parityOut.Length - 1; ++zIter)
                {
	                parityOut[zIter] = parityOut[zIter + 1];
                }

                parityOut[parityOut.Length - 1] = 0;
            }

            r = (byte)(z_0 ^ message[message.Length - 1]);

            for (var zIter = 0; zIter < parityOut.Length; ++zIter)
            {
	            parityOut[parityOut.Length - zIter - 1] ^= m_gf.Multiply(m_codeGenPoly[zIter], r);
            }
        }

        private Span<byte> BuildCodeGenPoly()
        {
	        var polys = new byte[m_parityLength][];

            // Build the degree-1 polynomials (we need 2t = numElements of them).
            // Eg 2t = 4, need four of them:
            //   (x + 1) is {1, 1}
            //   (x + 2) is {2, 1}
            //   (x + 4) is {4, 1}
            //   (x + 8) is {8, 1}

            // Remember that field[0] is 0, field[1] is a^0.
            for (var i = 0; i < m_parityLength; ++i)
            {
                polys[i] = new byte[] { m_gf.Field[i + 1], 1 };
            }

            // Multiply them one at a time to produce the field generator poly.
            var current = polys[0].AsSpan();
            for (var i = 1; i < m_parityLength; ++i)
            {
                current = m_gf.PolyMult(current, polys[i]);
            }

            return current;
        }
    }
}
