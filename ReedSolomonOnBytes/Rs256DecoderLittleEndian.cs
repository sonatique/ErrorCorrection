//Note: the original version of this code comes from https://github.com/antiduh/ErrorCorrection.git

using System;

namespace ReedSolomonOnBytes
{
    public readonly ref struct Rs256DecoderLittleEndian // Little-Endian because byte message[0] is the highest mathematical message element
    {
        private readonly GaloisField256 m_gf;
		private readonly int m_parityLength;

        private readonly Span<byte> m_syndroms;

        private readonly Span<byte> m_lambda;
        private readonly Span<byte> m_corrPoly;
        private readonly Span<byte> m_lambdaStar;

        private readonly Span<byte> m_lambdaPrime;

        private readonly Span<byte> m_omega;

        private readonly Span<byte> m_errorIndexes;

        private readonly Span<byte> m_chienCache;

        public Rs256DecoderLittleEndian(int codewordLength, int messageLength, in GaloisField256 gf)
        {
			m_parityLength = codewordLength - messageLength;

            m_gf = gf;

            // Syndrom calculation buffers
            m_syndroms = new byte[m_parityLength];

            // Lamda calculation buffers
            m_lambda = new byte[m_parityLength - 1];
            m_corrPoly = new byte[m_parityLength - 1];
            m_lambdaStar = new byte[m_parityLength - 1];

            // LambdaPrime calculation buffers
            m_lambdaPrime = new byte[m_parityLength - 2];

            // Omega calculation buffers
            m_omega = new byte[m_parityLength - 2];
            
            // Error position calculation
            m_errorIndexes = new byte[codewordLength];

            // Cache of the lookup used in the ChienSearch process.
            m_chienCache = new byte[codewordLength];

            for(var i = 0; i < m_chienCache.Length; i++)
            {
                m_chienCache[i] = m_gf.Inverses[m_gf.Field[i + 1]];
            }
        }

		public void CorrectInPlace(in Span<byte> messageAndParity)
        {
            CalcSyndromPoly(messageAndParity);
            CalcLambda();
            CalcLambdaPrime();
            CalcOmega();

            ChienSearch();

            RepairErrors(messageAndParity, m_errorIndexes, m_omega, m_lambdaPrime);
        }

        private void RepairErrors(Span<byte> messageAndParity, in ReadOnlySpan<byte> errorIndexes, in ReadOnlySpan<byte> omega, in ReadOnlySpan<byte> lp)
        {
	        var messageLen = messageAndParity.Length;

            for(var i = 0; i < messageLen; i++)
            {
	            if (errorIndexes[i] != 0)
	            {
		            continue;
	            }

	            var x = m_gf.Field[i + 1];

	            var xInverse = m_gf.Inverses[x];

	            var top = m_gf.PolyEval(omega, xInverse);
	            top = m_gf.Multiply(top, x);
	            var bottom = m_gf.PolyEval(lp, xInverse);
                    
	            messageAndParity[messageLen - 1 - i] ^= m_gf.Divide(top, bottom);
            }
        }

        private void CalcLambda()
        {
	        // --- Initial conditions ----
            // Need to clear lambda and corrPoly, but not lambdaStar. lambda and corrPoly 
            // are used and initialized iteratively in the algorithm, whereas lambdaStar isn't.
            m_corrPoly.Fill(0);
            m_lambda.Fill(0);
            
            var k = 1;
            var l = 0;
            m_corrPoly[1] = 1;
            m_lambda[0] = 1;
            
            while(k <= m_parityLength)
            {            
                // --- Calculate e ---
                var e = m_syndroms[k - 1];

                for(var i = 1; i <= l; i++)
                {
                    e ^= m_gf.Multiply(m_lambda[i], m_syndroms[k - 1 - i]);
                }

                // --- Update estimate if e != 0 ---
                if(e != 0)
                {
                    // D*(x) = D(x) + e * C(x);
                    for(var i = 0; i < m_lambdaStar.Length; i++)
                    {
                        m_lambdaStar[i] = (byte)(m_lambda[i] ^ m_gf.Multiply(e, m_corrPoly[i]));
                    }

                    if(2 * l < k)
                    {
                        // L = K - L;
                        l = k - l;

                        // C(x) = D(x) * e^(-1);
                        var eInv = m_gf.Inverses[e]; // temp to store calculation of 1 / e aka e^(-1)
                        for(var i = 0; i < m_corrPoly.Length; i++)
                        {
                            m_corrPoly[i] = m_gf.Multiply(m_lambda[i], eInv);
                        }
                    }
                }

                // --- Advance C(x) ---

                // C(x) = C(x) * x
                for(var i = m_corrPoly.Length - 1; i >= 1; i--)
                {
                    m_corrPoly[i] = m_corrPoly[i - 1];
                }

                m_corrPoly[0] = 0;

                if(e != 0)
                {
                    // D(x) = D*(x);
                    m_lambdaStar.CopyTo(m_lambda);
                }

                k += 1;
            }
        }

        private void CalcLambdaPrime()
        {
            // Forney's says that we can just set even powers to 0 and then take the rest and 
            // divide it by x (shift it down one). 
            
            // No need to clear this.lambdaPrime between calls; full assignment is done every call.

            for(var i = 0; i < m_lambdaPrime.Length; i++)
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
            for (var i = 0; i < m_omega.Length; i++)
            {
                m_omega[i] = m_syndroms[i];

                for (var lIter = 1; lIter <= i; lIter++)
                {
                    m_omega[i] ^= m_gf.Multiply(m_syndroms[i - lIter], m_lambda[lIter]);
                }
            }
        }

        private void ChienSearch()
        {
            for(var i = 0; i < m_errorIndexes.Length; i++)
            {
                m_errorIndexes[i] = m_gf.PolyEval(
                    m_lambda,
                    m_chienCache[i]
               );
            }
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
    }
}
