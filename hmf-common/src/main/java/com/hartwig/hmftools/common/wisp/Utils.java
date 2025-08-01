package com.hartwig.hmftools.common.wisp;

import static java.lang.String.format;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class Utils
{
    public static String generateMutationSequence(
            final RefGenomeInterface refGenome, final int probeLength,
            final String chromosome, final int position, final String ref, final String alt)
    {
        int altLength = alt.length();
        int refLength = ref.length();
        int startLength = probeLength / 2 - altLength / 2;
        int startPos = position - startLength;

        String basesStart = refGenome.getBaseString(chromosome, startPos, position - 1);
        int endBaseLength = probeLength - basesStart.length() - altLength;

        int postPosition = position + refLength;
        String basesEnd = refGenome.getBaseString(chromosome, postPosition, postPosition + endBaseLength - 1);

        String sequence = basesStart + alt + basesEnd;

        if(sequence.length() != probeLength)
        {
            throw new IllegalArgumentException(format("variant(%s:%d %s->%s) invalid sequenceLength(%d): %s",
                    chromosome, position, ref, alt, sequence.length(), sequence));
        }

        return sequence;
    }
}
