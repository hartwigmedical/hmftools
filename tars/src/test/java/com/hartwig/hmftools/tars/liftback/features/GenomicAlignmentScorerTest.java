package com.hartwig.hmftools.tars.liftback.features;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.primaryRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.selfAlignment;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.tars.common.BwaScoring;
import com.hartwig.hmftools.tars.liftback.LiftedAlignment;
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
        return new LiftedAlignment(chrom, pos, cigar, 0, false, forwardStrand, 0);
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

        scorer.scorePlacements(alignments, record);

        byte[] reverseComplement = Arrays.copyOf(forward, forward.length);
        SequenceUtil.reverseComplement(reverseComplement);

        assertEquals(50, self.GenomicScore);
        assertEquals(BwaScoring.genomicScore(genome.asRefGenome(), CHR_1, 2000, "50M", reverseComplement),
                oppositeStrandAlt.GenomicScore);
        assertTrue(BwaScoring.genomicScore(genome.asRefGenome(), CHR_1, 2000, "50M", forward)
                < oppositeStrandAlt.GenomicScore);
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

        scorer.scorePlacements(alignments, record);

        assertEquals(50, alignments.get(0).GenomicScore);
        assertEquals(50, alignments.get(1).GenomicScore);
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
}
