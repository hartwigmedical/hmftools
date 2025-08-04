package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_GAP_EXTEND_PENALTY;
import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_GAP_OPEN_PENALTY;
import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_MATCH_SCORE;
import static com.hartwig.hmftools.common.aligner.BwaParameters.BWA_MISMATCH_PENALTY;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.getDuplexIndels;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecordUnpaired;
import static com.hartwig.hmftools.redux.ReduxConstants.INVALID_BASE_QUAL;
import static com.hartwig.hmftools.redux.consensus.SbxRoutines.getAnnotatedBases;
import static com.hartwig.hmftools.redux.consensus.SbxRoutines.processAnnotatedBases;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;
import static htsjdk.samtools.SAMUtils.phredToFastq;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.hartwig.hmftools.common.genome.refgenome.CachedRefGenome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

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
