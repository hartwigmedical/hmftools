package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.FULL;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.buildCigarString;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadContextMatcher;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadCounter;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.common.ReadContextMatcher;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class RealignmentTest
{
    @Test
    public void testInsertRealignment()
    {
        String refBases = "X" + REF_BASES_200.substring(0, 100)
                // 0123456789012345678901234567890123456789
                //             ->
                + "GTTGTTGTTGTCGTTGTTGTTGTTGGTTTTTCTGAGACAGAGTC" + REF_BASES_200.substring(0, 100);

        RefSequence refSequence = new RefSequence(0, refBases.getBytes());

        String insertedBases = "TTGTTGTTGTTGGTTTTTC";

        String varBuildReadBases = refBases.substring(91, 101) + refBases.substring(101, 114) + insertedBases + refBases.substring(114, 150);

        String readCigar = "33M19I36M";

        SAMRecord varBuildRead = buildSamRecord(91, readCigar, varBuildReadBases);
        String ref = refBases.substring(113, 114);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        SimpleVariant var = createSimpleVariant(113, ref, ref + insertedBases);
        VariantReadContext readContext = builder.createContext(var, varBuildRead, 22, refSequence);

        assertEquals(113, readContext.variant().position());
        assertEquals(11, readContext.VarIndex);
        assertEquals("CGTTGTTGTTGTTGGTTTTTCTTGTTGTTGTTGGTTTTTCTGA", readContext.coreStr());

        // variant is at read index 41 when should be 22
        int realignedStartPos = 72;
        readCigar = buildCigarString(varBuildReadBases.length());
        SAMRecord realignedRead = buildSamRecord(realignedStartPos, readCigar, varBuildReadBases);

        ReadContextMatcher matcher = createReadContextMatcher(readContext);

        int realignedReadIndex = Realignment.realignedReadIndexPosition(readContext, realignedRead);
        ReadContextMatch matchType = matcher.determineReadMatch(realignedRead, realignedReadIndex);

        assertEquals(FULL, matchType);
    }

    @Test
    public void testSoftClipRealignment()
    {
        String refBases = REF_BASES_200.substring(0, 100)
                + "GGGGATTGGTCTATCTATCTATCTAGG" + REF_BASES_200.substring(0, 100);
                // 0123456789012345678901234567890123456789

        RefSequence refSequence = new RefSequence(0, refBases.getBytes());

        String insertedBases = "TCTA";

        int varPosition = 108;
        int readPosStart = 89;
        int varIndexInRead = varPosition - readPosStart;
        String varBuildReadBases = refBases.substring(readPosStart, varPosition + 1) + insertedBases + refBases.substring(varPosition + 1, 149);

        String readCigar = "20M4I40M";

        SAMRecord varBuildRead = buildSamRecord(readPosStart, readCigar, varBuildReadBases);
        String ref = refBases.substring(varPosition, varPosition + 1);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        SimpleVariant var = createSimpleVariant(varPosition, ref, ref + insertedBases);
        VariantReadContext readContext = builder.createContext(var, varBuildRead, varIndexInRead, refSequence);

        assertEquals(11, readContext.VarIndex);
        assertEquals("TCTA", readContext.Homology.Bases);
        assertEquals(126, readContext.CorePositionEnd);
        assertEquals("GGTCTATCTATCTATCTATCTAGG", readContext.coreStr());

        ReadContextMatcher matcher = createReadContextMatcher(readContext);

        // confirm the initial read is a full match
        ReadContextMatch matchType = matcher.determineReadMatch(varBuildRead, varIndexInRead);

        assertEquals(FULL, matchType);

        // position:        SC       94
        // position:        SC       4567890123456789012
        // index:           0123456789012345678901234567
        String readBases = "ATTGGTCTATCTATCTATCTATCTAGG" + REF_BASES_200.substring(0, 20);
        int realignedStartPos = readPosStart + 20;
        readCigar = "9S38M";
        SAMRecord realignedRead = buildSamRecord(realignedStartPos, readCigar, readBases);

        RawContext rawContext = RawContext.createFromRead(var, realignedRead);
        assertEquals(8, rawContext.ReadVariantIndex);

        int realignedReadIndex = Realignment.realignedReadIndexPosition(readContext, realignedRead);
        assertEquals(4, realignedReadIndex);
        matchType = matcher.determineReadMatch(realignedRead, realignedReadIndex);

        assertEquals(FULL, matchType);

        ReadContextCounter readContextCounter = createReadCounter(0, readContext);

        readContextCounter.processRead(realignedRead, 0, null);

        assertEquals(1, readContextCounter.readCounts().Realigned);
    }

    private static final String REF_BASE_START = "ACGTTGCAACGTTGCAGGGG"; // 20 bases

    @Test
    public void testInsertRealignment2()
    {
        String refBases = REF_BASE_START
                + "CTTTCTTTTCTTTCTTTTTCTTTA" + REF_BASES_200.substring(0, 50);
        // index   0123456789012345678901234567890123456789
        // pos     20        30        40

        RefSequence refSequence = new RefSequence(0, refBases.getBytes());

        String insertedBases = "TTTCTTTC";

        int varPosition = 35;
        int readPosStart = 20;
        int varIndexInRead = varPosition - readPosStart;
        String varBuildReadBases = refBases.substring(readPosStart, varPosition + 1) + insertedBases + refBases.substring(varPosition + 1);

        String readCigar = "16M8I20M";

        SAMRecord varBuildRead = buildSamRecord(readPosStart, readCigar, varBuildReadBases);
        String ref = refBases.substring(varPosition, varPosition + 1);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        SimpleVariant var = createSimpleVariant(varPosition, ref, ref + insertedBases);
        VariantReadContext readContext = builder.createContext(var, varBuildRead, varIndexInRead, refSequence);

        assertEquals(12, readContext.VarIndex);
        assertEquals("TTTCTTT", readContext.Homology.Bases);
        assertEquals(44, readContext.CorePositionEnd);
        assertEquals("CTTTTTCTTTCTTTCTTTAC", readContext.coreStr());

        ReadContextCounter readContextCounter = createReadCounter(0, readContext);

        readContextCounter.processRead(varBuildRead, 0, null);
        assertEquals(1, readContextCounter.readCounts().Full);

        RawContext rawContext = RawContext.createFromRead(var, varBuildRead);
        assertEquals(15, rawContext.ReadVariantIndex);

        ReadContextMatcher matcher = createReadContextMatcher(readContext);

        // test 1: basic realignment
        String readBases = "CTTTCTTTTTCTTTCTTTCTTTA" + REF_BASES_200.substring(0, 20); // 17 + 20
        int realignedStartPos = readPosStart;
        readCigar = "5M1I9M2D22M";
        SAMRecord realignedRead = buildSamRecord(realignedStartPos, readCigar, readBases);

        int realignedReadIndex = Realignment.realignedReadIndexPosition(readContext, realignedRead);
        ReadContextMatch matchType = matcher.determineReadMatch(realignedRead, realignedReadIndex);

        assertEquals(FULL, matchType);

        readContextCounter.processRead(realignedRead, 0, null);

        assertEquals(1, readContextCounter.readCounts().Realigned);


        // test 2: jitter realignment
        readBases = "CTTTCTTTTTCTTTCTTTCTTTCTTTA" + REF_BASES_200.substring(0, 20);
        readCigar = "9M1I10M2I25M";
        realignedRead = buildSamRecord(realignedStartPos, readCigar, readBases);

        realignedReadIndex = Realignment.realignedReadIndexPosition(readContext, realignedRead);
        matcher = createReadContextMatcher(readContext);
        RealignedType realignedType = Realignment.checkRealignment(
                readContext, matcher, realignedRead,  16, realignedReadIndex, null);

        // works if realignedReadIndex is 6 rather than 10, since the 1xTATT lengthened jitter causes the realignedIndex to be off
        assertEquals(RealignedType.LENGTHENED, realignedType);


        refBases = REF_BASE_START
                + "TTCTTTTCTTTCTTTTTCTTTA" + REF_BASES_200.substring(0, 60);
        // index   0123456789012345678901234567890123456789
        // pos     20        30        40

        refSequence = new RefSequence(0, refBases.getBytes());

        varPosition = 33;
        ref = refBases.substring(varPosition, varPosition + 1);
        var = createSimpleVariant(varPosition, ref, ref + insertedBases);

        readPosStart = 20; // drop back to make the flank long enough
        varIndexInRead = varPosition - readPosStart;
        varBuildReadBases = refBases.substring(readPosStart, varPosition + 1) + insertedBases + refBases.substring(varPosition + 1);
        readCigar = "12M8I30M";
        varBuildRead = buildSamRecord(readPosStart, readCigar, varBuildReadBases);

        readContext = builder.createContext(var, varBuildRead, varIndexInRead, refSequence);

        // 1:33 T>TTTTCTTTC read(TCTTTTCTTTCTTTTTCTTTCTTTCTTTACGCAATATTCG 11M8I20M) pos(21-0) index(10-12-29) repeat(14: TTTC-3) homology(TTTCTTT length(7)) alt(12-20) ref(CTTTTTCTTTAC)
        readBases = "TCTTTTTCTTTCTTTCTTTA" + REF_BASES_200.substring(0, 20); // 20 + 20
        realignedStartPos = readPosStart;
        readCigar = "12M2D28M";
        realignedRead = buildSamRecord(realignedStartPos, readCigar, readBases);

        realignedReadIndex = Realignment.realignedReadIndexPosition(readContext, realignedRead);
        matcher = createReadContextMatcher(readContext);
        matchType = matcher.determineReadMatch(realignedRead, realignedReadIndex);
        assertEquals(FULL, matchType);

        readContextCounter = createReadCounter(0, readContext);
        readContextCounter.processRead(realignedRead, 0, null);
        assertEquals(1, readContextCounter.readCounts().Realigned);


        varPosition = 36;
        varIndexInRead = varPosition - readPosStart;
        var = createSimpleVariant(varPosition, ref, ref + insertedBases);
        readContext = builder.createContext(var, varBuildRead, varIndexInRead, refSequence);

        readBases = "CTTTCTTTTTCTTTCTTTCTTTA" + REF_BASES_200.substring(0, 20); // 23 + 20
        realignedStartPos = 26;
        readCigar = "6S9M2D8M";
        realignedRead = buildSamRecord(realignedStartPos, readCigar, readBases);

        realignedReadIndex = Realignment.realignedReadIndexPosition(readContext, realignedRead);
        matcher = createReadContextMatcher(readContext);
        matchType = matcher.determineReadMatch(realignedRead, realignedReadIndex);
        // assertEquals(FULL, matchType);

        readContextCounter = createReadCounter(0, readContext);
        readContextCounter.processRead(realignedRead, 0, null);
       //  assertEquals(1, readContextCounter.readCounts().Realigned); // TODO
    }

    @Test
    public void testDelBeforeUnclippedStart()
    {
        String refBases = REF_BASES_200.substring(0, 100)
                + "GGGGATTGGTCTATCTATCTATCTAGG" + REF_BASES_200.substring(0, 100);
        //         0123456789012345678901234567890123456789

        RefSequence refSequence = new RefSequence(0, refBases.getBytes());

        String deletedBases = "TCTATCTA";

        int varPosition = 108;
        int readPosStart = 89;
        int varIndexInRead = varPosition - readPosStart;
        String varBuildReadBases = refBases.substring(readPosStart, varPosition + 1) + refBases.substring(varPosition + 1 + deletedBases.length(), 149);

        String readCigar = "20M8D32M";

        SAMRecord varBuildRead = buildSamRecord(readPosStart, readCigar, varBuildReadBases);
        String ref = refBases.substring(varPosition, varPosition + 1);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        SimpleVariant var = createSimpleVariant(varPosition, ref + deletedBases, ref);
        VariantReadContext readContext = builder.createContext(var, varBuildRead, varIndexInRead, refSequence);

        ReadContextCounter readContextCounter = createReadCounter(0, readContext);

        String readBases = "TGGTCTATCTAGG" + REF_BASES_200.substring(0, 20);
        readCigar = "3S30M";
        SAMRecord realignedRead = buildSamRecord(varPosition + deletedBases.length() + 1, readCigar, readBases);

        // core end position for this read is 126
        // index in read for this should be:
        // 3S   index 0-2   pre position 116
        // 30M  index 3-32  positions 116-145

        readContextCounter.processRead(realignedRead, 0, null);
        assertEquals(1, readContextCounter.readCounts().Realigned);
    }
}
