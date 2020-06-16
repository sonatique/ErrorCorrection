//Note: the original version of this code comes from https://github.com/antiduh/ErrorCorrection.git
using System;
using System.Runtime.CompilerServices;

namespace ReedSolomonOnBytes
{
    public readonly struct Rs256DecoderLittleEndian // Little-Endian because byte message[0] is the highest mathematical message element
    {
        private readonly GaloisField256 m_gf;
		
        private readonly byte[] m_syndroms;

        private readonly byte[] m_lambda;
        
        private readonly byte[] m_lambdaPrime;

        private readonly byte[] m_omega;

        private readonly byte[] m_errorIndexes;
        
        private readonly byte[] m_chienCache;

        public Rs256DecoderLittleEndian(int codewordLength, int messageLength, in GaloisField256 gf)
        {
	        CodewordLength = codewordLength;
	        MessageLength = messageLength;

	        ParityLength = codewordLength - messageLength;
			
            m_gf = gf;

            // Syndrom calculation buffers
            m_syndroms = new byte[ParityLength];

            // Lamda calculation buffers
            m_lambda = new byte[ParityLength - 1];
            
            // LambdaPrime calculation buffers
            m_lambdaPrime = new byte[ParityLength - 2];

            // Omega calculation buffers
            m_omega = new byte[ParityLength - 2];
            
            // Error position calculation
            m_errorIndexes = new byte[codewordLength];

            // Cache of the lookup used in the ChienSearch process.
            m_chienCache = new byte[codewordLength];

            for(var i = 0; i < m_chienCache.Length; ++i)
            {
                m_chienCache[i] = m_gf.Inverses[m_gf.Field[i + 1]];
            }
        }

        public int ParityLength { get; }
        public int CodewordLength { get; }
        public int MessageLength { get; }

		public int CorrectInPlace(in Span<byte> messageAndParity)
		{
			if (messageAndParity.Length == 0)
			{
				return 0;
			}

			var errorIndexes = GetErrorIndexes(messageAndParity);

            return errorIndexes.Length > 0 ? RepairErrors(messageAndParity, errorIndexes, m_omega, m_lambdaPrime) : 0;
        }

		private int RepairErrors(in Span<byte> messageAndParity, in ReadOnlySpan<byte> errorIndexes, in ReadOnlySpan<byte> omega, in ReadOnlySpan<byte> lp)
        {
			var messageLen = messageAndParity.Length;

	        var repairedBitCount = 0;

	        var corrections = new Rs256Algorithms.FixedList<(byte errorIndex, byte xorValue)>(stackalloc (byte errorIndex, byte xorValue)[CodewordLength]);

            GetCorrections(messageAndParity.Length, errorIndexes, omega, lp, ref corrections);

            foreach (var error in corrections)
            {
	            var repairingValue = error.xorValue;

	            repairedBitCount += PopulationCount(repairingValue);
                    
	            messageAndParity[messageLen - 1 - error.errorIndex] ^= repairingValue;
            }

            return repairedBitCount;
        }

		private void GetCorrections(int messageLen, in ReadOnlySpan<byte> errorIndexes, in ReadOnlySpan<byte> omega, in ReadOnlySpan<byte> lp, ref Rs256Algorithms.FixedList<(byte errorIndex, byte xorValue)> corrections)
		{
			foreach (var errorIndex in errorIndexes)
			{
				if (errorIndex >= messageLen)
				{
					continue; // repair would occur in leading 0 bytes
				}

				var x = m_gf.Field[errorIndex + 1];

				var xInverse = m_gf.Inverses[x];

				var top = m_gf.PolyEval(omega, xInverse);

				top = m_gf.Multiply(top, x);
	            
				var bottom = m_gf.PolyEval(lp, xInverse);

				var repairingValue = m_gf.Divide(top, bottom);

				corrections.Add((errorIndex, repairingValue));
			}
		}

		private ReadOnlySpan<byte> GetErrorIndexes(in ReadOnlySpan<byte> messageAndParity)
		{
			CalcSyndromPoly(messageAndParity);
			CalcLambda();
			CalcLambdaPrime();
			CalcOmega();

			return ChienSearch();
		}

        private void CalcLambda()
        {
            // Need to clear lambda 
            Span<byte> lambda = m_lambda;
            lambda.Fill(0);
            lambda[0] = 1;
            
            Span<byte> corrPoly = stackalloc byte[ParityLength - 1];
            corrPoly[1] = 1;
			
            var l = 0;
            
            Span<byte> lambdaStar = stackalloc byte[ParityLength - 1];
            
            for(var k = 1; k <= ParityLength; ++k)
            {            
                // --- Calculate e ---
                var e = m_syndroms[k - 1];

                for(var i = 1; i <= l; ++i)
                {
                    e ^= m_gf.Multiply(lambda[i], m_syndroms[k - 1 - i]);
                }

                // --- Update estimate if e != 0 ---
                if(e != 0)
                {
                    // D*(x) = D(x) + e * C(x);
                    for(var i = 0; i < lambdaStar.Length; ++i)
                    {
                        lambdaStar[i] = (byte)(lambda[i] ^ m_gf.Multiply(e, corrPoly[i]));
                    }

                    if(2 * l < k)
                    {
	                    // L = K - L;
                        l = k - l;

                        if (l >= lambda.Length)
                        {
	                        l = lambda.Length - 1;
                        }

                        // C(x) = D(x) * e^(-1);
                        var eInv = m_gf.Inverses[e]; // temp to store calculation of 1 / e aka e^(-1)
                        for(var i = 0; i < corrPoly.Length; ++i)
                        {
                            corrPoly[i] = m_gf.Multiply(lambda[i], eInv);
                        }
                    }
                }

                // --- Advance C(x) ---

                // C(x) = C(x) * x
                for(var i = corrPoly.Length - 1; i >= 1; i--)
                {
                    corrPoly[i] = corrPoly[i - 1];
                }

                corrPoly[0] = 0;

                if(e != 0)
                {
                    // D(x) = D*(x);
                    lambdaStar.CopyTo(lambda);
                }
            }
        }

        private void CalcLambdaPrime()
        {
            // Forney's says that we can just set even powers to 0 and then take the rest and 
            // divide it by x (shift it down one). 
            
            // No need to clear this.lambdaPrime between calls; full assignment is done every call.

            for(var i = 0; i < m_lambdaPrime.Length; ++i)
            {
                if((i & 0x1) == 0)
                {
                    m_lambdaPrime[i] = m_lambda[i + 1];
                }
                else
                {
                    m_lambdaPrime[i] = 0;
                }
            }
        }

        private void CalcOmega()
        {
            for (var i = 0; i < m_omega.Length; ++i)
            {
                m_omega[i] = m_syndroms[i];

                for (var lIter = 1; lIter <= i; lIter++)
                {
                    m_omega[i] ^= m_gf.Multiply(m_syndroms[i - lIter], m_lambda[lIter]);
                }
            }
        }

        private ReadOnlySpan<byte> ChienSearch()
        {
	        var errorIndexCount = 0;

	        ReadOnlySpan<byte> chienCache = m_chienCache;

            for(var i = 0; i < chienCache.Length; ++i)
            {
	            var error = m_gf.PolyEval(
		            m_lambda,
		            chienCache[i]
	            );

	            if (error != 0)
	            {
                    continue;
	            }

	            m_errorIndexes[errorIndexCount++] = (byte)i;
            }

            return new ReadOnlySpan<byte>(m_errorIndexes, 0, errorIndexCount);
        }
        
        private void CalcSyndromPoly(in ReadOnlySpan<byte> message)
        {
	        // Don't need to zero this.syndromes first - it's not used before its assigned to.

            for(var synIndex = 0; synIndex < m_syndroms.Length; ++synIndex)
            {
                // EG, if g(x) = (x+a^0)(x+a^1)(x+a^2)(x+a^3) 
                //             = (x+1)(x+2)(x+4)(x+8),
                // Then for the first syndrom S_0, we would provide root = a^0 = 1.
                // S_1 --> root = a^1 = 2 etc.

                var root = m_gf.Field[synIndex + 1];
                byte syndrome = 0;

                for(var i = 0; i < message.Length -1; ++i)
                {
                    syndrome = m_gf.Multiply((byte)(syndrome ^ message[i]), root);
                }

                m_syndroms[synIndex] = (byte)(syndrome ^ message[message.Length - 1]);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static int PopulationCount(byte x)
        {
	        var x1 = x - (x >> 1 & 0b1010101);
	        var x2 = (x1 >> 2 & 0b110011) + (x1 & 0b110011);
	        return (x2 >> 4) + (x2 & 0b1111);
        }
    }
}
