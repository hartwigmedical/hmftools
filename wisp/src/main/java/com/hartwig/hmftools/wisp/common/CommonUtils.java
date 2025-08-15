package com.hartwig.hmftools.wisp.common;

import static java.lang.String.format;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class CommonUtils
{
    public static final Logger CT_LOGGER = LogManager.getLogger(CommonUtils.class);

    public static final String APP_NAME = "Wisp";

    // common fields
    public static final String FLD_TUMOR_ID = "TumorId";
    public static final String FLD_CATEGORY = "Category";
    public static final String FLD_VARIANT = "Variant";

    public static final String BATCH_CONTROL_TAG = "batch_control";

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
            throw new IllegalArgumentException(format("variant(%s:%d %s->%s) invalid sequenceLength(%d): %s",
                    chromosome, position, ref, alt, sequence.length(), sequence));
        }

        return sequence;
    }
}
