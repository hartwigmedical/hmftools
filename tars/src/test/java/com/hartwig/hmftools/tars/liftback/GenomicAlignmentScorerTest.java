package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.primaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.selfAlignment;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

public class GenomicAlignmentScorerTest
{
    private static TestGenome genome()
    {
        return new TestGenome().with(CHR_1, 5000, 'A');
    }

    private static LiftedAlignment refAlt(final String chrom, final int pos, final String cigar, final boolean forwardStrand)
    {
        return new LiftedAlignment(chrom, pos, cigar, 0, false, forwardStrand, null);
    }

    @Test
    public void testScoresPlacementsInGenomeSpaceOnTheirOwnStrands()
    {
        TestGenome genome = genome().set(CHR_1, 2000, "T".repeat(50));
        GenomicAlignmentScorer scorer = new GenomicAlignmentScorer(genome.asRefGenome());

        LiftedAlignment self = selfAlignment(CHR_1, 100, "50M");
        LiftedAlignment oppositeStrandAlt = refAlt(CHR_1, 2000, "50M", false);
        List<LiftedAlignment> alignments = new ArrayList<>(List.of(self, oppositeStrandAlt));

        SAMRecord record = primaryRecord(CHR_1, 100, "50M");
        byte[] forward = bases("A".repeat(50));
        record.setReadBases(forward);

        List<LiftedAlignment> scored = scorer.scorePlacements(alignments, record);

        byte[] reverseComplement = Arrays.copyOf(forward, forward.length);
        SequenceUtil.reverseComplement(reverseComplement);

        assertEquals(50, scored.get(0).genomicScore().intValue());
        assertEquals(
                GenomicAlignmentScorer.genomicAlignmentScore(genome.asRefGenome(), CHR_1, 2000, cigarFromStr("50M"), reverseComplement),
                scored.get(1).genomicScore().intValue());
        assertTrue(GenomicAlignmentScorer.genomicAlignmentScore(genome.asRefGenome(), CHR_1, 2000, cigarFromStr("50M"), forward)
                < scored.get(1).genomicScore());
    }

    @Test
    public void testSupplementaryTagDoesNotSuppressPlacementScoring()
    {
        GenomicAlignmentScorer scorer = new GenomicAlignmentScorer(genome().asRefGenome());

        List<LiftedAlignment> alignments = new ArrayList<>(
                List.of(selfAlignment(CHR_1, 100, "50M"), refAlt(CHR_1, 2000, "50M", true)));
        SAMRecord record = primaryRecord(CHR_1, 100, "50M");
        record.setReadBases(bases("A".repeat(50)));
        record.setAttribute("SA", CHR_1 + ",2000,+,50M,0,0;");

        List<LiftedAlignment> scored = scorer.scorePlacements(alignments, record);

        assertEquals(50, scored.get(0).genomicScore().intValue());
        assertEquals(50, scored.get(1).genomicScore().intValue());
    }

    @Test
    public void testEqualsAndMismatchCigarsScoreUpToTheChromosomeEnd()
    {
        GenomicAlignmentScorer scorer = new GenomicAlignmentScorer(genome().asRefGenome());
        SAMRecord record = primaryRecord(CHR_1, 4991, "10M");
        record.setReadBases(bases("A".repeat(10)));

        assertEquals(10, scorer.scoreRecord(refAlt(CHR_1, 4991, "10=", true), record));
        assertEquals(10, scorer.scoreRecord(refAlt(CHR_1, 4991, "10X", true), record));
        assertThrows(IllegalStateException.class, () -> scorer.scoreRecord(refAlt(CHR_1, 4992, "10=", true), record));
        assertThrows(IllegalStateException.class, () -> scorer.scoreRecord(refAlt(CHR_1, 4992, "10X", true), record));
    }

    @Test
    public void testCigarNeedingMoreBasesThanTheReadIsRejected()
    {
        GenomicAlignmentScorer scorer = new GenomicAlignmentScorer(genome().asRefGenome());
        SAMRecord record = primaryRecord(CHR_1, 100, "10M");
        record.setReadBases(bases("A".repeat(8)));

        assertThrows(IllegalStateException.class, () -> scorer.scoreRecord(refAlt(CHR_1, 100, "10M", true), record));
    }

    @Test
    public void testIntronOnlyLiftDoesNotRequireRescore()
    {
        SAMRecord record = primaryRecord(CHR_1, 100, "100M");
        LiftedAlignment lifted = refAlt(CHR_1, 100, "50M100N50M", true);

        assertFalse(GenomicAlignmentScorer.requiresRescore(record, lifted));
    }

    @Test
    public void testChangedScoringOperationsRequireRescore()
    {
        SAMRecord record = primaryRecord(CHR_1, 100, "100M");
        LiftedAlignment lifted = refAlt(CHR_1, 100, "95M5S", true);

        assertTrue(GenomicAlignmentScorer.requiresRescore(record, lifted));
    }

    @Test
    public void testStrandFlipComparesReversedCigarOperations()
    {
        SAMRecord record = primaryRecord(CHR_1, 100, "10S40M");
        record.setReadNegativeStrandFlag(true);
        LiftedAlignment lifted = refAlt(CHR_1, 100, "40M10S", true);

        assertFalse(GenomicAlignmentScorer.requiresRescore(record, lifted));
    }

    @Test
    public void testReferenceWindowsUseOneBasedInclusiveCoordinatesAndFixtureLengths()
    {
        RefGenomeInterface reference = new TestGenome().with(CHR_1, bases("ACGTTGCA")).asRefGenome();

        assertEquals(8, reference.getChromosomeLength(CHR_1));
        assertEquals(Map.of(CHR_1, 8), reference.chromosomeLengths());
        assertArrayEquals(bases("ACG"), GenomicAlignmentScorer.referenceBases(reference, CHR_1, 1, 3));
        assertArrayEquals(bases("GCA"), GenomicAlignmentScorer.referenceBases(reference, CHR_1, 6, 8));
        assertNull(GenomicAlignmentScorer.referenceBases(reference, CHR_1, 0, 3));
        assertNull(GenomicAlignmentScorer.referenceBases(reference, CHR_1, 7, 9));
        assertNull(GenomicAlignmentScorer.referenceBases(reference, CHR_1, 4, 3));
    }

    @Test
    public void testScoresShortOverhangsLikeBwaMem()
    {
        byte[] reference = bases("A".repeat(10));

        // TARS keeps a short anchor only when its recomputed bwa-mem score reaches the overhang threshold.
        assertEquals(5, GenomicAlignmentScorer.ungappedScore(bases("A".repeat(9) + "C"), reference));

        assertEquals(5, GenomicAlignmentScorer.bestScoringPrefixLength(bases("A".repeat(5) + "C".repeat(3)), reference));

        assertEquals(0, GenomicAlignmentScorer.bestScoringPrefixLength(bases("C".repeat(5)), reference));

        // keep extending past an internal mismatch when later matches recover the best score
        assertEquals(10, GenomicAlignmentScorer.bestScoringPrefixLength(bases("A".repeat(5) + "C" + "A".repeat(4)), reference));
    }

    @Test
    public void testScoresGenomicPlacementFromCigar()
    {
        TestGenome genome = new TestGenome().with(CHR_1, 5000, 'A');

        assertEquals(
                5, GenomicAlignmentScorer.genomicAlignmentScore(
                        genome.asRefGenome(), CHR_1, 1, cigarFromStr("10M"), bases("AAAACAAAAA")));
        assertEquals(
                3, GenomicAlignmentScorer.genomicAlignmentScore(
                        genome.asRefGenome(), CHR_1, 1, cigarFromStr("5M2D5M"), bases("AAAAAAAAAA")));
        assertEquals(
                10, GenomicAlignmentScorer.genomicAlignmentScore(
                        genome.asRefGenome(), CHR_1, 1, cigarFromStr("5M100N5M"), bases("AAAAAAAAAA")));
        assertEquals(
                7, GenomicAlignmentScorer.genomicAlignmentScore(
                        genome.asRefGenome(), CHR_1, 1, cigarFromStr("3S7M"), bases("AAAAAAAAAA")));
    }
}
