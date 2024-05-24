package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.FULL;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.buildCigarString;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
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

        ReadContextMatcher matcher = new ReadContextMatcher(readContext);

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

        ReadContextMatcher matcher = new ReadContextMatcher(readContext);

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

        RawContext rawContext = RawContext.create(var, realignedRead);
        assertEquals(8, rawContext.ReadVariantIndex);

        int realignedReadIndex = Realignment.realignedReadIndexPosition(readContext, realignedRead);
        assertEquals(4, realignedReadIndex);
        matchType = matcher.determineReadMatch(realignedRead, realignedReadIndex);

        assertEquals(FULL, matchType);

        ReadContextCounter readContextCounter = createReadCounter(0, readContext);

        readContextCounter.processRead(realignedRead, 0, null);

        assertEquals(1, readContextCounter.readCounts().Realigned);
    }
}
