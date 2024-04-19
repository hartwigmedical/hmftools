package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.CORE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.FULL;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.PARTIAL_CORE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.REF;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.buildCigarString;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadContext;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReadContextMatchingTest
{
    private static final String READ_FLANK_BASES = REF_BASES_200.substring(0, 20);

    @Test
    public void testReadCoversCore()
    {
        int position = 100;
        SimpleVariant variant = createSimpleVariant(position, "A", "T");

        VariantReadContext readContext = createReadContext(variant);

        ReadContextMatcher matcher = new ReadContextMatcher(readContext);

        String readBases = REF_BASES_200.substring(0, 40);
        byte[] readQualities = SamRecordTestUtils.buildDefaultBaseQuals(readBases.length());
        String cigar = buildCigarString(readBases.length());

        SAMRecord read = buildSamRecord(1, cigar, readBases, readQualities);

        assertTrue(matcher.coversVariant(read, 20));

        assertTrue(matcher.coversVariant(read, 1));
        assertTrue(matcher.coversVariant(read, 0));
        assertTrue(matcher.coversVariant(read, 39));

        // now test with an indel with a wider alt range
        SimpleVariant var = createSimpleVariant(position, "C", "CCC");
        readContext = createReadContext(var, "ATGTGTG", "CCA");
        matcher = new ReadContextMatcher(readContext);

        assertTrue(matcher.coversVariant(read, 20));
        assertTrue(matcher.coversVariant(read, 7));
        assertFalse(matcher.coversVariant(read, 6));
        assertTrue(matcher.coversVariant(read, 36));
        assertFalse(matcher.coversVariant(read, 37));
    }

    @Test
    public void testRefAndCoreMatches()
    {
        String leftCore = "AC";
        String rightCore = "GT";
        int position = 100;
        String ref = "A";
        String alt = "T";
        SimpleVariant variant = createSimpleVariant(position, ref, alt);
        VariantReadContext readContext = createReadContext(variant, leftCore, rightCore);

        ReadContextMatcher matcher = new ReadContextMatcher(readContext);

        String readBases = readContext.leftFlankStr() + leftCore + ref + rightCore + readContext.rightFlankStr();
        byte[] readQualities = buildDefaultBaseQuals(readBases.length());
        int readVarIndex = readContext.leftFlankLength() + 2;
        String cigar = buildCigarString(readBases.length());

        SAMRecord read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);

        assertEquals(REF, matcher.determineReadMatch(read, readVarIndex));

        // full match
        readBases = readContext.leftFlankStr() + leftCore + alt + rightCore + readContext.rightFlankStr();
        read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);
        assertEquals(FULL, matcher.determineReadMatch(read, readVarIndex));

        // low-qual mismatch but not at variant
        readBases = readContext.leftFlankStr() + "AAAGT" + readContext.rightFlankStr();
        readQualities[readVarIndex - 1] = 10;
        read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);
        assertEquals(REF, matcher.determineReadMatch(read, readVarIndex));

        // too many mismatches
        readBases = readContext.leftFlankStr() + "AAACT" + readContext.rightFlankStr();
        readQualities[readVarIndex - 1] = 10;
        readQualities[readVarIndex + 1] = 10;
        read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);
        assertEquals(NONE, matcher.determineReadMatch(read, readVarIndex));
        assertNotEquals(REF, matcher.determineReadMatch(read, readVarIndex));

        // mismatch at SNV itself cannot be low-qual
        readBases = readContext.leftFlankStr() + "ACGGT" + readContext.rightFlankStr();
        readQualities = buildDefaultBaseQuals(readBases.length());
        readQualities[readVarIndex] = 10;
        read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);
        assertEquals(NONE, matcher.determineReadMatch(read, readVarIndex));

        // can be low-qual at the variant base itself
        readBases = readContext.leftFlankStr() + leftCore + ref + rightCore + readContext.rightFlankStr();
        read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);
        assertEquals(REF, matcher.determineReadMatch(read, readVarIndex));

        // partial core
        readBases = "AAAAA" + readContext.leftFlankStr() + leftCore + alt;
        readVarIndex = readBases.length() - 1;
        readQualities = buildDefaultBaseQuals(readBases.length());
        read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);
        assertEquals(PARTIAL_CORE, matcher.determineReadMatch(read, readVarIndex));

        // same again from the other side
        readBases =  alt + rightCore + readContext.rightFlankStr() + "AAAAA";
        readVarIndex = 0;
        readQualities = buildDefaultBaseQuals(readBases.length());
        read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);
        assertEquals(PARTIAL_CORE, matcher.determineReadMatch(read, readVarIndex));
    }

    @Test
    public void testFlankMatches()
    {
        String leftCore = "AC";
        String rightCore = "GT";
        int position = 100;
        String ref = "A";
        String alt = "T";
        SimpleVariant variant = createSimpleVariant(position, ref, alt);
        VariantReadContext readContext = createReadContext(variant, leftCore, rightCore);

        ReadContextMatcher matcher = new ReadContextMatcher(readContext);

        String flankMismatch = "AAAAATTTTT";
        String readBases = readContext.leftFlankStr() + leftCore + alt + rightCore + flankMismatch;
        byte[] readQualities = buildDefaultBaseQuals(readBases.length());
        int readVarIndex = readContext.leftFlankLength() + 2;
        String cigar = buildCigarString(readBases.length());

        SAMRecord read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);

        assertEquals(CORE, matcher.determineReadMatch(read, readVarIndex));

        // low-qual mismatches in the flanks
        readBases = flankMismatch + leftCore + alt + rightCore + flankMismatch;
        readQualities = buildDefaultBaseQuals(readBases.length());

        for(int i = 0; i < 10; ++i)
        {
            readQualities[i] = 10;
        }

        for(int i = readBases.length() - 11; i < readBases.length(); ++i)
        {
            readQualities[i] = 10;
        }

        read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);

        assertEquals(FULL, matcher.determineReadMatch(read, readVarIndex));

        // missing flank still classifies as FULL if there are no mismatches
        readBases = readContext.leftFlankStr() + leftCore + alt + rightCore;
        readQualities = buildDefaultBaseQuals(readBases.length());
        cigar = buildCigarString(readBases.length());

        read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);

        assertEquals(FULL, matcher.determineReadMatch(read, readVarIndex));
    }

}