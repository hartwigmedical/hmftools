package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_READ_CONTEXT_FLANK_SIZE;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReadContextClassifierTest
{
    private static int READ_LENGTH = 100;
    private static final byte[] QUALITIES = new byte[READ_LENGTH];
    static
    {
        for(int i = 0; i < QUALITIES.length; ++i)
        {
            QUALITIES[i] = 37;
        }
    }

    private static final VariantReadContext SNV_CONTEXT;
    static
    {
        String chromosome = CHR_1;
        int position = 100;
        int readStartPos = 50;
        String refBases = REF_BASES_200;

        String ref = String.valueOf(REF_BASES_200.charAt(position - 1));
        String alt = "A";

        // TODO:
//        assertNotEquals(ref, alt);

        int varReadIndex = position - readStartPos;

        StringBuilder readBases = new StringBuilder(refBases.substring(readStartPos - 1, readStartPos - 1 + READ_LENGTH));
        readBases.setCharAt(varReadIndex, alt.charAt(0));
        SAMRecord read = buildSamRecord(readStartPos, READ_LENGTH + "M", readBases.toString(), QUALITIES);

        SimpleVariant variant = new SimpleVariant(chromosome, position, ref, alt);
        RefSequence refSequence = new RefSequence(1, refBases.getBytes());
        VariantReadContextBuilder variantReadContextBuilder = new VariantReadContextBuilder(DEFAULT_READ_CONTEXT_FLANK_SIZE);
        SNV_CONTEXT = variantReadContextBuilder.createMnvContext(variant, read, varReadIndex, refSequence);
    }

    @Test
    public void testUnmappedRead()
    {
        SAMRecord read = new SAMRecord(null);
        read.setReadUnmappedFlag(true);
        ReadContextClassifier classifier = new ReadContextClassifier(SNV_CONTEXT);
        ReadContextCounter.MatchType matchType = classifier.classifyRead(read);

        assertNull(matchType);
    }

    @Test
    public void testNonoverlappingRead()
    {
        SAMRecord read = buildSamRecord(SNV_CONTEXT.AlignmentEnd + 1, READ_LENGTH + "M", "A".repeat(READ_LENGTH), QUALITIES);
        ReadContextClassifier classifier = new ReadContextClassifier(SNV_CONTEXT);
        ReadContextCounter.MatchType matchType = classifier.classifyRead(read);

        assertNull(matchType);
    }

    @Test
    public void testOverlappingRead()
    {
        SAMRecord read = buildSamRecord(SNV_CONTEXT.AlignmentStart, READ_LENGTH + "M", "A".repeat(READ_LENGTH), QUALITIES);
        ReadContextClassifier classifier = new ReadContextClassifier(SNV_CONTEXT);
        ReadContextCounter.MatchType matchType = classifier.classifyRead(read);

        assertEquals(ReadContextCounter.MatchType.NONE, matchType);
    }
}
