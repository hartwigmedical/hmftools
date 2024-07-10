package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.CORE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.FULL;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.PARTIAL_CORE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.SIMPLE_ALT;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.REF;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_SEQUENCE_200;
import static com.hartwig.hmftools.sage.common.TestUtils.buildCigarString;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.VariantUtils.TEST_LEFT_CORE;
import static com.hartwig.hmftools.sage.common.VariantUtils.TEST_RIGHT_CORE;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadContext;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadContextMatcher;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.averageCoreQuality;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReadContextMatchingTest
{
    @Test
    public void testReadCoversCore()
    {
        int position = 100;
        SimpleVariant variant = createSimpleVariant(position, "A", "T");

        VariantReadContext readContext = createReadContext(variant);

        ReadContextMatcher matcher = createReadContextMatcher(readContext);

        String readBases = REF_BASES_200.substring(0, 40);
        byte[] readQualities = SamRecordTestUtils.buildDefaultBaseQuals(readBases.length());
        String cigar = buildCigarString(readBases.length());

        SAMRecord read = buildSamRecord(1, cigar, readBases, readQualities);

        assertTrue(matcher.coversVariant(read, 20));
        assertTrue(matcher.coversVariant(read, 0));
        assertTrue(matcher.coversVariant(read, 39));

        // now test with an indel with a wider alt range
        SimpleVariant var = createSimpleVariant(position, "C", "CCC");
        readContext = createReadContext(var, "ATGTGTG", "CCA");
        matcher = createReadContextMatcher(readContext);

        assertTrue(matcher.coversVariant(read, 20));
        assertTrue(matcher.coversVariant(read, 0));
        assertTrue(matcher.coversVariant(read, 36));
        assertFalse(matcher.coversVariant(read, 37));

        var = createSimpleVariant(position, "ACGT", "A");
        readContext = createReadContext(var, "AA", "GG");
        matcher = createReadContextMatcher(readContext);

        assertTrue(matcher.coversVariant(read, 20));
        assertTrue(matcher.coversVariant(read, 0));
        assertTrue(matcher.coversVariant(read, 38));
        assertFalse(matcher.coversVariant(read, 39));
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

        ReadContextMatcher matcher = createReadContextMatcher(readContext);

        String readBases = readContext.leftFlankStr() + leftCore + ref + rightCore + readContext.rightFlankStr();
        byte[] readQualities = buildDefaultBaseQuals(readBases.length());
        int readVarIndex = readContext.leftFlankLength() + 2;
        String cigar = buildCigarString(readBases.length());

        SAMRecord read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);

        assertEquals(REF, matcher.determineReadMatch(read, readVarIndex));

        // partial ref core assuming satisfies core-covered
        readBases = ref + rightCore + readContext.rightFlankStr();
        read = buildSamRecord(position, cigar, readBases, readQualities);
        assertEquals(REF, matcher.determineReadMatch(read, 0));

        // same again but ending at the ref
        readBases = "ACGTA" + readContext.leftFlankStr() + leftCore + ref;
        read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);
        assertEquals(REF, matcher.determineReadMatch(read, readBases.length() - 1));

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
    public void testSimpleAltMatches()
    {
        String leftCore = "AC";
        String rightCore = "GT";
        int position = 100;
        String ref = "AGT";
        String alt = "TGC";
        SimpleVariant variant = createSimpleVariant(position, ref, alt);
        VariantReadContext readContext = createReadContext(variant, leftCore, rightCore);

        ReadContextMatcher matcher = new ReadContextMatcher(readContext, true, true);

        String readBases = readContext.leftFlankStr() + "TT" + alt + rightCore + readContext.rightFlankStr();
        byte[] readQualities = buildDefaultBaseQuals(readBases.length());
        int readVarIndex = readContext.leftFlankLength() + 2;
        String cigar = buildCigarString(readBases.length());

        SAMRecord read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);

        assertEquals(SIMPLE_ALT, matcher.determineReadMatch(read, readVarIndex));

        readBases = readContext.leftFlankStr() + "TT" + "TCC" + rightCore + readContext.rightFlankStr();
        read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);

        assertEquals(NONE, matcher.determineReadMatch(read, readVarIndex));

        // repeated for an SNV
        ref = "A";
        alt = "T";
        variant = createSimpleVariant(position, ref, alt);
        readContext = createReadContext(variant, leftCore, rightCore);

        matcher = new ReadContextMatcher(readContext, true, true);

        readBases = readContext.leftFlankStr() + leftCore + alt + "CC" + readContext.rightFlankStr();
        readQualities = buildDefaultBaseQuals(readBases.length());
        readVarIndex = readContext.leftFlankLength() + 2;
        cigar = buildCigarString(readBases.length());

        read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);

        assertEquals(SIMPLE_ALT, matcher.determineReadMatch(read, readVarIndex));


        // test 3: delete
        ref = "ATG";
        alt = "A";
        variant = createSimpleVariant(position, ref, alt);
        readContext = createReadContext(variant, leftCore, rightCore);

        matcher = new ReadContextMatcher(readContext, true, true);

        readBases = readContext.leftFlankStr() + leftCore + alt + "CC" + readContext.rightFlankStr();
        readQualities = buildDefaultBaseQuals(readBases.length());
        readVarIndex = readContext.leftFlankLength() + 2;

        cigar = "13M2D12M";
        read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);

        assertEquals(SIMPLE_ALT, matcher.determineReadMatch(read, readVarIndex));

        cigar = "14M2D12M";
        read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);

        assertEquals(NONE, matcher.determineReadMatch(read, readVarIndex));

        // test 4: insert
        ref = "A";
        alt = "ATTT";
        variant = createSimpleVariant(position, ref, alt);
        readContext = createReadContext(variant, leftCore, rightCore);

        matcher = new ReadContextMatcher(readContext, true, true);

        readBases = readContext.leftFlankStr() + leftCore + alt + "CC" + readContext.rightFlankStr();
        readQualities = buildDefaultBaseQuals(readBases.length());
        readVarIndex = readContext.leftFlankLength() + 2;
        cigar = "13M3I12M";

        read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);

        assertEquals(SIMPLE_ALT, matcher.determineReadMatch(read, readVarIndex));

        cigar = "13M4I12M";
        read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);

        assertEquals(NONE, matcher.determineReadMatch(read, readVarIndex));
    }

    @Test
    public void testIndelMatches()
    {
        // basic del ref match
        int position = 50;
        SimpleVariant var = createSimpleVariant(
                position, REF_BASES_200.substring(position, position + 5), REF_BASES_200.substring(position, position + 1));

        String readBases = REF_BASES_200.substring(30, position + 1) + REF_BASES_200.substring(position + 5, 70);
        String readCigar = "21M2D17M";
        SAMRecord read = buildSamRecord(30, readCigar, readBases);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        VariantReadContext readContext = builder.createContext(var, read, 20, REF_SEQUENCE_200);

        ReadContextMatcher matcher = createReadContextMatcher(readContext);

        String refBases = REF_BASES_200.substring(30, 70);
        SAMRecord refRead = buildSamRecord(30, "40M", refBases);

        assertEquals(REF, matcher.determineReadMatch(refRead, 20));
    }

    @Test
    public void testIndelHomologyMatches()
    {
        int position = 64;
        String ref = REF_BASES_200.substring(position, position + 1);
        SimpleVariant var = createSimpleVariant(position, ref, ref + "AAAA");

        String readBases = REF_BASES_200.substring(40, position + 1) + "AAAA" + REF_BASES_200.substring(position + 1, 90);
        String readCigar = "25M4I25M";
        SAMRecord read = buildSamRecord(40, readCigar, readBases);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        VariantReadContext readContext = builder.createContext(var, read, 24, REF_SEQUENCE_200);

        ReadContextMatcher matcher = createReadContextMatcher(readContext);

        assertEquals(FULL, matcher.determineReadMatch(read, 24));

        // values: [20, 23, 24, 11, 12, 15]
        // now with a read with low-qual mismatches in non-permitted locations

        readBases = REF_BASES_200.substring(40, position + 1) + "ACAA" + REF_BASES_200.substring(position + 1, 90);
        readCigar = "25M4I25M";
        read = buildSamRecord(40, readCigar, readBases);
        read.getBaseQualities()[26] = 11;

        assertEquals(FULL, matcher.determineReadMatch(read, 24));

        readBases = REF_BASES_200.substring(40, position + 1) + "AAAC" + REF_BASES_200.substring(position + 1, 90);
        readCigar = "25M4I25M";
        read = buildSamRecord(40, readCigar, readBases);
        read.getBaseQualities()[28] = 11;

        assertEquals(NONE, matcher.determineReadMatch(read, 24));

        // 1:64 G>GAAAA read(TCTGTGACTC-GGAAAAAAAAAAAACT-CCCTGACCCC 12M4I20M) pos(53-84) index(10-11-25) repeat(A-12 index(12-23)) homology(AAAA length=8) alt(11-24) ref(GGAAAAAAAACT)
        readBases = REF_BASES_200.substring(40, position + 1) + "AAAA" + REF_BASES_200.substring(position + 1, position + 8)
                + "T" + REF_BASES_200.substring(position + 9, 90);

        readCigar = "25M4I25M";
        read = buildSamRecord(40, readCigar, readBases);
        read.getBaseQualities()[36] = 11;

        assertEquals(NONE, matcher.determineReadMatch(read, 24));
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

        ReadContextMatcher matcher = createReadContextMatcher(readContext);

        String flankMismatch = "AAAAATTTTT";
        String readBases = readContext.leftFlankStr() + leftCore + alt + rightCore + flankMismatch;
        byte[] readQualities = buildDefaultBaseQuals(readBases.length());
        int readVarIndex = readContext.leftFlankLength() + 2;
        String cigar = buildCigarString(readBases.length());

        SAMRecord read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);

        assertEquals(CORE, matcher.determineReadMatch(read, readVarIndex));

        // low-qual mismatches in the flanks - 3 are permitted
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

        assertEquals(CORE, matcher.determineReadMatch(read, readVarIndex));

        // now limited to 3
        flankMismatch = readContext.rightFlankStr().substring(0, 7) + "TTT";
        readBases = readContext.leftFlankStr() + leftCore + alt + rightCore + flankMismatch;
        readQualities = buildDefaultBaseQuals(readBases.length());

        for(int i = readBases.length() - 4; i < readBases.length(); ++i)
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

    @Test
    public void testAverageCoreQuality()
    {
        String leftCore = TEST_LEFT_CORE;
        String rightCore = TEST_RIGHT_CORE;
        int position = 100;
        String ref = "A";
        String alt = "T";
        SimpleVariant variant = createSimpleVariant(position, ref, alt);

        VariantReadContext readContext = createReadContext(variant, leftCore, rightCore);
        ReadContextMatcher matcher = createReadContextMatcher(readContext);

        String readBases = readContext.leftFlankStr() + leftCore + alt + rightCore + readContext.rightFlankStr();
        byte[] readQualities = buildDefaultBaseQuals(readBases.length());
        int readVarIndex = readContext.leftFlankLength() + leftCore.length();

        readQualities[readVarIndex - 2] = 10;
        readQualities[readVarIndex - 1] = 12;
        readQualities[readVarIndex] = 14;
        readQualities[readVarIndex + 1] = 16;
        readQualities[readVarIndex + 2] = 18;

        String cigar = buildCigarString(readBases.length());

        SAMRecord read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);

        double average = averageCoreQuality(readContext, read, readVarIndex);
        assertEquals(14, average, 0.01);

        // now with partial cores
        readBases = readContext.leftFlankStr() + leftCore + alt;
        readQualities = buildDefaultBaseQuals(readBases.length());

        readQualities[readVarIndex - 2] = 12;
        readQualities[readVarIndex - 1] = 13;
        readQualities[readVarIndex] = 14;

        read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);

        average = averageCoreQuality(readContext, read, readVarIndex);
        assertEquals(13, average, 0.01);

        // partial on the left
        readBases = alt + rightCore + readContext.rightFlankStr();
        readQualities = buildDefaultBaseQuals(readBases.length());

        readQualities[0] = 10;
        readQualities[1] = 11;
        readQualities[2] = 12;

        readVarIndex = 0;
        read = buildSamRecord(position - readVarIndex, cigar, readBases, readQualities);

        average = averageCoreQuality(readContext, read, readVarIndex);
        assertEquals(11, average, 0.01);
    }
}