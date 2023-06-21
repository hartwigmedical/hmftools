package com.hartwig.hmftools.ctdna.common;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.ctdna.probe.PvConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class CommonUtils
{
    public static final Logger CT_LOGGER = LogManager.getLogger(CommonUtils.class);

    public static final int DEFAULT_PROBE_LENGTH = 120;

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
            CT_LOGGER.error("variant({}:{} {}->{}) invalid sequenceLength({}): {}",
                    chromosome, position, ref, alt, sequence.length(), sequence);
        }

        return sequence;
    }

    public static double calcGcPercent(final String sequence)
    {
        if(sequence.isEmpty())
            return 0;

        int gcCount = 0;
        for(int i = 0; i < sequence.length(); ++i)
        {
            if(sequence.charAt(i) == 'G' || sequence.charAt(i) == 'C')
                ++gcCount;
        }

        return gcCount / (double)sequence.length();
    }
}
