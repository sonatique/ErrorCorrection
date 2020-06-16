using System;

namespace ReedSolomonOnBytes
{
	public static class Rs256Algorithms
	{
		public static void CorrectFecInPlace(in Span<byte> messageAndParityWithUncertainLength, 
			in Rs256EncoderLittleEndian encoder, in Rs256DecoderLittleEndian decoder, 
			int probableTotalMessageLength, int possibleNumberOfAdditionalMessageBytes,
			out int totalFecedLength,
			out int totalCorrectedBits)
		{
			totalFecedLength = GetTotalFecedLength(messageAndParityWithUncertainLength, encoder, decoder, probableTotalMessageLength, possibleNumberOfAdditionalMessageBytes);

			var lastBlockLength = totalFecedLength % decoder.CodewordLength;
			
			totalCorrectedBits = 0;

			var currentStartIndex = totalFecedLength - lastBlockLength; 

			totalCorrectedBits += decoder.CorrectInPlace(messageAndParityWithUncertainLength.Slice(currentStartIndex, lastBlockLength));

			currentStartIndex -= decoder.CodewordLength;
			
			while (currentStartIndex >= 0)
			{
				totalCorrectedBits += decoder.CorrectInPlace(messageAndParityWithUncertainLength.Slice(currentStartIndex, decoder.CodewordLength));
				currentStartIndex -= decoder.CodewordLength;
			}
		}

		private static int GetTotalFecedLength(in ReadOnlySpan<byte> messageAndParityWithUncertainLength, 
			in Rs256EncoderLittleEndian encoder, in Rs256DecoderLittleEndian decoder, 
			int probableTotalMessageLength, int possibleNumberOfAdditionalMessageBytes)
		{
			var chunkMaxLength = messageAndParityWithUncertainLength.Length % encoder.CodewordLength;

			var currentStartIndex = messageAndParityWithUncertainLength.Length - chunkMaxLength;

			var lengthList = new FixedList<int>(stackalloc int[encoder.CodewordLength]);

			while (currentStartIndex >= 0)
			{
				SearchLengths(
					messageAndParityWithUncertainLength.Slice(currentStartIndex, chunkMaxLength),
					encoder, decoder,
					ref lengthList);

				if (lengthList.Length > 0) // there was not only junk in the last supposed part
				{
					break;
				}

				chunkMaxLength = encoder.CodewordLength;
				currentStartIndex -= chunkMaxLength;

				lengthList.Clear();
			}

			if (currentStartIndex < 0)
			{
				return messageAndParityWithUncertainLength.Length;
			}

			var lastBlockLength = lengthList[0];

			if (lengthList.Length > 1)
			{
				var smallestDif = encoder.CodewordLength;

				foreach (var candidateLengthOfLastBlock in lengthList)
				{
					var candidateTotalLength = currentStartIndex + candidateLengthOfLastBlock;

					if (probableTotalMessageLength <= candidateTotalLength)
					{
						var diff = candidateTotalLength - probableTotalMessageLength;
						if (diff < smallestDif)
						{
							smallestDif = diff;
							lastBlockLength = candidateLengthOfLastBlock;
						}
					}
					else if(probableTotalMessageLength + possibleNumberOfAdditionalMessageBytes <= candidateTotalLength)
					{
						var diff = candidateTotalLength - (probableTotalMessageLength + possibleNumberOfAdditionalMessageBytes);
						if (diff < smallestDif)
						{
							smallestDif = diff;
							lastBlockLength = candidateLengthOfLastBlock;
						}
					}

				}
			}

			return currentStartIndex + lastBlockLength;
		}

		public static void RemoveAndCorrectFecInPlace(ref Span<byte> messageAndParityWithUncertainLength, in Rs256EncoderLittleEndian encoder, in Rs256DecoderLittleEndian decoder, int probableTotalMessageLength, int possibleNumberOfAdditionalMessageBytes, out int totalCorrectedBits)
		{
			CorrectFecInPlace(in messageAndParityWithUncertainLength, encoder, decoder, probableTotalMessageLength, possibleNumberOfAdditionalMessageBytes, out var totalFecedLength, out totalCorrectedBits);

			var fullBlockCount = totalFecedLength / decoder.CodewordLength;
			var lastBlockLength = totalFecedLength % decoder.CodewordLength;

			for (var i = 0; i < fullBlockCount - 1; ++i)
			{
				messageAndParityWithUncertainLength.Slice((i+1) * decoder.CodewordLength, decoder.MessageLength)
					.CopyTo(messageAndParityWithUncertainLength.Slice((i+1) * decoder.MessageLength, decoder.MessageLength));
			}

			messageAndParityWithUncertainLength.Slice(fullBlockCount * decoder.CodewordLength, lastBlockLength)
				.CopyTo(messageAndParityWithUncertainLength.Slice(fullBlockCount * decoder.MessageLength, lastBlockLength));

			messageAndParityWithUncertainLength = messageAndParityWithUncertainLength.Slice(0, fullBlockCount * decoder.MessageLength + lastBlockLength - decoder.ParityLength);
		}

		private static int ProbableTotalMessageBlockLength(int allFecedDataLength, int codewordLength, int parityLength)
		{
			var numberOfBlock = allFecedDataLength / codewordLength;
			if (allFecedDataLength % codewordLength > 0)
			{
				++numberOfBlock;
			}

			return allFecedDataLength - numberOfBlock * parityLength;
		}
	
		public static void SearchLengths(in ReadOnlySpan<byte> messageAndParityWithUncertainLength, in Rs256EncoderLittleEndian encoder, in Rs256DecoderLittleEndian decoder, ref FixedList<int> foundLengthsOutput)
		{
			if (messageAndParityWithUncertainLength.Length <= encoder.ParityLength)
			{
				if (messageAndParityWithUncertainLength.Length == encoder.ParityLength && messageAndParityWithUncertainLength.IsAllZero())
				{
					foundLengthsOutput.Add(0);
				}
				return;
			}

			Span<byte> expectedParity = stackalloc byte[encoder.ParityLength];

			Span<byte> copy = stackalloc byte[messageAndParityWithUncertainLength.Length];

			messageAndParityWithUncertainLength.CopyTo(copy);

			var iterationCount = messageAndParityWithUncertainLength.Length - encoder.ParityLength; // we could limit to a smaller number of possible extra byte

			for (var i = 0; i < iterationCount; ++i)
			{
				var candidateData = copy.Slice(0, copy.Length - i);

				var correctedBitCount = decoder.CorrectInPlace(candidateData);

				encoder.CalculateParity(candidateData.Slice(0, candidateData.Length - encoder.ParityLength), expectedParity);

				var candidateParity = candidateData.Slice(candidateData.Length - encoder.ParityLength, encoder.ParityLength);

				if (candidateParity.SequenceEqual(expectedParity))
				{
					var originalParity =
						messageAndParityWithUncertainLength.Slice(candidateData.Length - encoder.ParityLength, encoder.ParityLength);

					var checkingSpecialCaseInvalidateCandidate = false;
					for (var p = encoder.ParityLength - 1; p >= 0; --p)
					{
						if (candidateParity[p] != 0 || originalParity[p] == 0)
						{
							continue;
						}

						checkingSpecialCaseInvalidateCandidate = true;
						break;
					}

					if (!checkingSpecialCaseInvalidateCandidate)
					{
						if(!foundLengthsOutput.Add(messageAndParityWithUncertainLength.Length - i))
						{
							return;
						}
					}
				}

				if (correctedBitCount > 0)
				{
					messageAndParityWithUncertainLength.CopyTo(copy);
				}
			}

			for (var i = 0; i < encoder.ParityLength - 1; ++i)
			{
				var truncatedParityTestLength = encoder.ParityLength - 1 - i;

				encoder.CalculateParity(messageAndParityWithUncertainLength.Slice(0, messageAndParityWithUncertainLength.Length - truncatedParityTestLength), expectedParity);

				var candidateParity = messageAndParityWithUncertainLength.Slice(messageAndParityWithUncertainLength.Length - truncatedParityTestLength, truncatedParityTestLength);

				//if (expectedParity.StartsWith(candidateParity))
				if(expectedParity[0] == candidateParity[0]) // very loose test
				{
					//we take only the largest one
					foundLengthsOutput.Add(messageAndParityWithUncertainLength.Length - truncatedParityTestLength + encoder.ParityLength); 
					return;
				}
			}
		}

		private static bool IsAllZero(in this ReadOnlySpan<byte> span)
		{
			foreach (var t in span)
			{
				if (t != 0)
				{
					return false;
				}
			}

			return true;
		}

		public ref struct FixedList<T>
		{
			private readonly Span<T> m_data;

			public FixedList(in Span<T> data)
			{
				m_data  = data;
				Length = 0;
			}

			public void Clear()
			{
				Length = 0;
			}

			public readonly T this[int index] => m_data[index];

			public bool Add(T i)
			{
				if (Length == m_data.Length)
				{
					return false;
				}

				m_data[Length++] = i;

				return true;
			}

			public int Length { get; private set; }

			public readonly Span<T>.Enumerator GetEnumerator()
			{
				return m_data.Slice(0, Length).GetEnumerator();
			}

			public readonly T[] ToArray()
			{
				return m_data.Slice(0, Length).ToArray();
			}
		}
	}
}
