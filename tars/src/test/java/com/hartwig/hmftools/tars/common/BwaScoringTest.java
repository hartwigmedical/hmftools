package com.hartwig.hmftools.tars.common;

import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertNull;

import java.util.Map;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;

import org.junit.Test;

public class BwaScoringTest
{
    @Test
    public void testReferenceWindowsUseOneBasedInclusiveCoordinatesAndFixtureLengths()
    {
        RefGenomeInterface reference = new TestGenome().with(CHR_1, bases("ACGTTGCA")).asRefGenome();

        assertEquals(8, reference.getChromosomeLength(CHR_1));
        assertEquals(Map.of(CHR_1, 8), reference.chromosomeLengths());
        assertArrayEquals(bases("ACG"), BwaScoring.refWindow(reference, CHR_1, 1, 3));
        assertArrayEquals(bases("GCA"), BwaScoring.refWindow(reference, CHR_1, 6, 8));
        assertNull(BwaScoring.refWindow(reference, CHR_1, 0, 3));
        assertNull(BwaScoring.refWindow(reference, CHR_1, 7, 9));
        assertNull(BwaScoring.refWindow(reference, CHR_1, 4, 3));
    }

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
    }
}
