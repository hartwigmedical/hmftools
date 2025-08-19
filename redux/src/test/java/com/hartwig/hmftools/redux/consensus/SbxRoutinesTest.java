package com.hartwig.hmftools.redux.consensus;

import static htsjdk.samtools.SAMUtils.phredToFastq;

public class SbxRoutinesTest
{
    /*
    @Test
    public void testFillQualZeroMismatchesWithRefSingleSNV()
    {
        int alignmentStart = 25;
        int mapq = 10;
        int nm = 1;
        int alignmentScore = 0;
        String readStr = "A".repeat(10) + "C" + "A".repeat(9);
        String cigar = "20M";
        String qualStr = NON_ZERO_QUAL.repeat(10) + ZERO_QUAL.repeat(1) + NON_ZERO_QUAL.repeat(9);

        SAMRecord read = createSamRecordUnpaired("READ_001", CHR_1, alignmentStart, readStr, cigar, false, false, null);
        read.setBaseQualityString(qualStr);
        read.setMappingQuality(mapq);
        read.setAttribute(NUM_MUTATONS_ATTRIBUTE, nm);
        read.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, alignmentScore);

        String refBases = "A".repeat(alignmentStart - 1 + 10) + "G" + "A".repeat(1000);
        RefGenomeInterface refGenome = getRefGenome(refBases);

        SbxRoutines.finaliseRead(refGenome, read);

        int alignmentScoreDiff = BWA_MISMATCH_PENALTY + BWA_MATCH_SCORE;

        int expectedNm = 0;
        int expectedAlignmentScore = alignmentScore + alignmentScoreDiff;
        String expectedReadStr = "A".repeat(10) + "G" + "A".repeat(9);
        String expectedCigar = "20M";
        String expectedQualStr = NON_ZERO_QUAL.repeat(10) + phredToFastq(1) + NON_ZERO_QUAL.repeat(9);

        SAMRecord expectedRead =
                createSamRecordUnpaired(
                        "READ_001", CHR_1, alignmentStart, expectedReadStr, expectedCigar, false, false, null);
        expectedRead.setBaseQualityString(expectedQualStr);
        expectedRead.setMappingQuality(mapq);
        expectedRead.setAttribute(NUM_MUTATONS_ATTRIBUTE, expectedNm);
        expectedRead.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, expectedAlignmentScore);

        assertEquals(expectedRead, read);
    }
    */

}
