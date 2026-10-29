package com.hartwig.hmftools.tars.common;

import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;

import org.junit.Test;

public class BwaScoringTest
{
    @Test
    public void testScoresShortOverhangsLikeBwaMem()
    {
        byte[] reference = bases("A".repeat(10));

        // TARS keeps a short anchor only when its recomputed bwa-mem score reaches the overhang threshold.
        assertEquals(5, BwaScoring.score(bases("A".repeat(9) + "C"), reference));

        assertEquals(5, BwaScoring.maxScoringPrefix(bases("A".repeat(5) + "C".repeat(3)), reference));

        assertEquals(0, BwaScoring.maxScoringPrefix(bases("C".repeat(5)), reference));

        // keep extending past an internal mismatch when later matches recover the best score
        assertEquals(10, BwaScoring.maxScoringPrefix(bases("A".repeat(5) + "C" + "A".repeat(4)), reference));
    }

    @Test
    public void testScoresGenomicPlacementFromCigar()
    {
        TestGenome genome = new TestGenome().with(CHR_1, 5000, 'A');

        assertEquals(5, BwaScoring.genomicScore(genome.asRefGenome(), CHR_1, 1, "10M", bases("AAAACAAAAA")));
        assertEquals(3, BwaScoring.genomicScore(genome.asRefGenome(), CHR_1, 1, "5M2D5M", bases("AAAAAAAAAA")));
        assertEquals(10, BwaScoring.genomicScore(genome.asRefGenome(), CHR_1, 1, "5M100N5M", bases("AAAAAAAAAA")));
        assertEquals(7, BwaScoring.genomicScore(genome.asRefGenome(), CHR_1, 1, "3S7M", bases("AAAAAAAAAA")));
        assertEquals(Integer.MIN_VALUE, BwaScoring.genomicScore(null, CHR_1, 1, "10M", bases("AAAAAAAAAA")));
    }
}
