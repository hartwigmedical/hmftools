package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_READ_CONTEXT_FLANK_SIZE;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.region.BaseRegion;
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

    // TODO: Other types of variants.
    private static final int READ_START_POS = 50;
    private static final String REF_READ_STRING;
    private static final String SNV_READ_STRING;
    private static final VariantReadContext SNV_CONTEXT;
    static
    {
        String chromosome = CHR_1;
        int position = 100;
        String refBases = REF_BASES_200;

        String ref = String.valueOf(REF_BASES_200.charAt(position - 1));
        String alt = swapDnaBase(ref);

        REF_READ_STRING = refBases.substring(READ_START_POS - 1, READ_START_POS - 1 + READ_LENGTH);

        int varReadIndex = position - READ_START_POS;
        StringBuilder readString = new StringBuilder(REF_READ_STRING);
        readString.setCharAt(varReadIndex, alt.charAt(0));
        SNV_READ_STRING = readString.toString();
        SAMRecord read = buildSamRecord(READ_START_POS, READ_LENGTH + "M", SNV_READ_STRING, QUALITIES);

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
    public void testNonOverlappingRead()
    {
        SAMRecord read = buildSamRecord(SNV_CONTEXT.AlignmentEnd + 1, READ_LENGTH + "M", "A".repeat(READ_LENGTH), QUALITIES);
        ReadContextClassifier classifier = new ReadContextClassifier(SNV_CONTEXT);
        ReadContextCounter.MatchType matchType = classifier.classifyRead(read);

        assertNull(matchType);
    }

    @Test
    public void testNoMatchingBasesSnv()
    {
        StringBuilder readString = new StringBuilder(SNV_READ_STRING);
        for(int i = 0; i < readString.length(); ++i)
        {
            readString.setCharAt(i, swapDnaBase(readString.charAt(i)));
        }

        SAMRecord read = buildSamRecord(READ_START_POS, READ_LENGTH + "M", readString.toString(), QUALITIES);
        ReadContextClassifier classifier = new ReadContextClassifier(SNV_CONTEXT);
        ReadContextCounter.MatchType matchType = classifier.classifyRead(read);

        assertEquals(ReadContextCounter.MatchType.NONE, matchType);
    }

    @Test
    public void testFullMatchSnvCoreAndRightFlank()
    {
        StringBuilder readString = new StringBuilder(SNV_READ_STRING);
        for(int pos = READ_START_POS; pos < SNV_CONTEXT.AlignmentStart + SNV_CONTEXT.leftFlankLength(); ++pos)
        {
            int idx = pos - READ_START_POS;
            readString.setCharAt(idx, swapDnaBase(readString.charAt(idx)));
        }

        SAMRecord read = buildSamRecord(READ_START_POS, READ_LENGTH + "M", readString.toString(), QUALITIES);
        ReadContextClassifier classifier = new ReadContextClassifier(SNV_CONTEXT);
        ReadContextCounter.MatchType matchType = classifier.classifyRead(read);

        assertEquals(ReadContextCounter.MatchType.FULL, matchType);
    }

    @Test
    public void testPartialMatchSnvPartialCoreAndLeftFlank()
    {
        StringBuilder readString = new StringBuilder(SNV_READ_STRING);
        for(int pos = SNV_CONTEXT.VarReadIndex + SNV_CONTEXT.AlignmentStart + 1; pos < READ_START_POS + READ_LENGTH; ++pos)
        {
            int idx = pos - READ_START_POS;
            readString.setCharAt(idx, swapDnaBase(readString.charAt(idx)));
        }

        SAMRecord read = buildSamRecord(READ_START_POS, READ_LENGTH + "M", readString.toString(), QUALITIES);
        ReadContextClassifier classifier = new ReadContextClassifier(SNV_CONTEXT);
        ReadContextCounter.MatchType matchType = classifier.classifyRead(read);

        assertEquals(ReadContextCounter.MatchType.PARTIAL, matchType);
    }

    @Test
    public void testCoreMatchSnv()
    {
        BaseRegion coreRegion = new BaseRegion(SNV_CONTEXT.AlignmentStart + SNV_CONTEXT.CoreIndexStart, SNV_CONTEXT.AlignmentStart + SNV_CONTEXT.CoreIndexEnd);
        StringBuilder readString = new StringBuilder(SNV_READ_STRING);
        for(int pos = READ_START_POS; pos < READ_START_POS + READ_LENGTH; ++pos)
        {
            if (coreRegion.containsPosition(pos))
            {
                continue;
            }

            int idx = pos - READ_START_POS;
            readString.setCharAt(idx, swapDnaBase(readString.charAt(idx)));
        }

        SAMRecord read = buildSamRecord(READ_START_POS, READ_LENGTH + "M", readString.toString(), QUALITIES);
        ReadContextClassifier classifier = new ReadContextClassifier(SNV_CONTEXT);
        ReadContextCounter.MatchType matchType = classifier.classifyRead(read);

        assertEquals(ReadContextCounter.MatchType.CORE, matchType);
    }

    @Test
    public void testRefMatchSnvPartialCoreAndLeftFlank()
    {
        StringBuilder readString = new StringBuilder(REF_READ_STRING);
        for(int pos = SNV_CONTEXT.VarReadIndex + SNV_CONTEXT.AlignmentStart + 1; pos < READ_START_POS + READ_LENGTH; ++pos)
        {
            int idx = pos - READ_START_POS;
            readString.setCharAt(idx, swapDnaBase(readString.charAt(idx)));
        }

        SAMRecord read = buildSamRecord(READ_START_POS, READ_LENGTH + "M", readString.toString(), QUALITIES);
        ReadContextClassifier classifier = new ReadContextClassifier(SNV_CONTEXT);
        ReadContextCounter.MatchType matchType = classifier.classifyRead(read);

        assertEquals(ReadContextCounter.MatchType.REF, matchType);
    }

    @Test
    public void testRefMatchSnvCore()
    {
        BaseRegion coreRegion = new BaseRegion(SNV_CONTEXT.AlignmentStart + SNV_CONTEXT.CoreIndexStart, SNV_CONTEXT.AlignmentStart + SNV_CONTEXT.CoreIndexEnd);
        StringBuilder readString = new StringBuilder(REF_READ_STRING);
        for(int pos = READ_START_POS; pos < READ_START_POS + READ_LENGTH; ++pos)
        {
            if (coreRegion.containsPosition(pos))
            {
                continue;
            }

            int idx = pos - READ_START_POS;
            readString.setCharAt(idx, swapDnaBase(readString.charAt(idx)));
        }

        SAMRecord read = buildSamRecord(READ_START_POS, READ_LENGTH + "M", readString.toString(), QUALITIES);
        ReadContextClassifier classifier = new ReadContextClassifier(SNV_CONTEXT);
        ReadContextCounter.MatchType matchType = classifier.classifyRead(read);

        assertEquals(ReadContextCounter.MatchType.REF, matchType);
    }
}
