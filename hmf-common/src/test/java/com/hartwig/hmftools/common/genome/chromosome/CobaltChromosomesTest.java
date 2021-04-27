package com.hartwig.hmftools.common.genome.chromosome;

import static com.hartwig.hmftools.common.cobalt.CobaltTestUtils.female;
import static com.hartwig.hmftools.common.cobalt.CobaltTestUtils.male;
import static com.hartwig.hmftools.common.cobalt.CobaltTestUtils.create;
import static com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes.MIN_Y_COUNT;
import static com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes.MOSIAC_X_CUTOFF;
import static com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes.TRISOMY_CUTOFF;
import static com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes.TWO_X_CUTOFF;
import static com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes.TWO_Y_CUTOFF;
import static com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes.Y_CUTOFF;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.utils.Doubles;

import org.junit.Test;

public class CobaltChromosomesTest {

    @Test
    public void testDefault() {
        final CobaltChromosomes male = male();
        final CobaltChromosomes female = female();

        assertEquals(Gender.MALE, male.gender());
        assertEquals(Gender.FEMALE, female.gender());

        for (int i = 0; i < 22; i++) {
            assertTrue(male.contains(i + ""));
            assertTrue(female.contains(i + ""));

            assertEquals(1, male.get(i + "").typicalRatio(), 0.01);
            assertEquals(1, female.get(i + "").typicalRatio(), 0.01);
        }

        assertTrue(male.contains("X"));
        assertTrue(female.contains("X"));
        assertEquals(0.5, male.get("X").typicalRatio(), 0.01);
        assertEquals(1, female.get("X").typicalRatio(), 0.01);

        assertTrue(male.contains("Y"));
        assertFalse(female.contains("Y"));
        assertEquals(0.5, male.get("Y").typicalRatio(), 0.01);
    }

    @Test
    public void testFemale() {
        MedianRatio chrX = create("X", 1, MIN_Y_COUNT);
        MedianRatio chrY = create("Y", 0.05, 100);
        List<MedianRatio> chromosomes = ImmutableList.of(chrX, chrY);
        CobaltChromosomes victim = new CobaltChromosomes(chromosomes);
        assertEquals(Gender.FEMALE, victim.gender());
        assertFalse(victim.hasGermlineAberrations());

        CobaltChromosome victimX = victim.get("X");
        assertTrue(victimX.isNormal());
        assertTrue(victimX.isDiploid());

        assertEquals(1, victim.get("X").actualRatio(), 0.01);
        assertFalse(victim.contains("Y"));
    }

    @Test
    public void testMale() {
        MedianRatio chrX = create("X", 0.5, MIN_Y_COUNT);
        MedianRatio chrY = create("Y", Y_CUTOFF, MIN_Y_COUNT);
        List<MedianRatio> chromosomes = ImmutableList.of(chrX, chrY);
        CobaltChromosomes victim = new CobaltChromosomes(chromosomes);
        assertEquals(Gender.MALE, victim.gender());
        assertFalse(victim.hasGermlineAberrations());

        CobaltChromosome victimX = victim.get("X");
        assertTrue(victimX.isNormal());
        assertFalse(victimX.isDiploid());

        assertEquals(0.5, victim.get("X").actualRatio(), 0.01);
        assertEquals(0.5, victim.get("Y").actualRatio(), 0.01);
    }

    @Test
    public void testKlinefelter() {
        MedianRatio chrX = create("X", TWO_X_CUTOFF, MIN_Y_COUNT);
        MedianRatio chrY = create("Y", Y_CUTOFF, MIN_Y_COUNT);
        List<MedianRatio> chromosomes = ImmutableList.of(chrX, chrY);
        CobaltChromosomes victim = new CobaltChromosomes(chromosomes);
        assertEquals(Gender.MALE, victim.gender());
        assertEquals(1, victim.germlineAberrations().size());
        assertTrue(victim.germlineAberrations().contains(GermlineAberration.KLINEFELTER));

        assertEquals(1, victim.get("X").actualRatio(), 0.01);
        assertEquals(0.5, victim.get("Y").actualRatio(), 0.01);
    }

    @Test
    public void testFemaleEvenWithSmallYCount() {
        MedianRatio chrX = create("X", 1, MIN_Y_COUNT);
        MedianRatio chrY = create("Y", Y_CUTOFF, MIN_Y_COUNT - 1);
        List<MedianRatio> chromosomes = ImmutableList.of(chrX, chrY);
        CobaltChromosomes victim = new CobaltChromosomes(chromosomes);
        assertEquals(Gender.FEMALE, victim.gender());
        assertFalse(victim.hasGermlineAberrations());
    }

    @Test
    public void testFemaleEvenWithSmallYMedian() {
        MedianRatio chrX = create("X", 1, MIN_Y_COUNT);
        MedianRatio chrY = create("Y", Y_CUTOFF - Y_CUTOFF / 10d, MIN_Y_COUNT);
        List<MedianRatio> chromosomes = ImmutableList.of(chrX, chrY);
        CobaltChromosomes victim = new CobaltChromosomes(chromosomes);
        assertEquals(Gender.FEMALE, victim.gender());
        assertFalse(victim.hasGermlineAberrations());
    }

    @Test
    public void testMosiacX() {
        assertTrue(Doubles.lessThan(TWO_X_CUTOFF, MOSIAC_X_CUTOFF));
        MedianRatio chr1 = create("1", 1, MIN_Y_COUNT);
        MedianRatio chrX = create("X", TWO_X_CUTOFF, MIN_Y_COUNT);
        List<MedianRatio> chromosomes = ImmutableList.of(chr1, chrX);
        CobaltChromosomes victim = new CobaltChromosomes(chromosomes);
        assertEquals(Gender.FEMALE, victim.gender());
        assertEquals(1, victim.germlineAberrations().size());
        assertTrue(victim.germlineAberrations().contains(GermlineAberration.MOSAIC_X));

        assertEquals(1, victim.get("1").actualRatio(), 0.01);
        assertEquals(TWO_X_CUTOFF, victim.get("X").actualRatio(), 0.01);
    }

    @Test
    public void testGCBiasNotMosiacX() {
        assertTrue(Doubles.lessThan(TWO_X_CUTOFF, MOSIAC_X_CUTOFF));
        MedianRatio chr1 = create("1", TWO_X_CUTOFF, MIN_Y_COUNT);
        MedianRatio chrX = create("X", TWO_X_CUTOFF, MIN_Y_COUNT);
        List<MedianRatio> chromosomes = ImmutableList.of(chr1, chrX);
        CobaltChromosomes victim = new CobaltChromosomes(chromosomes);
        assertEquals(Gender.FEMALE, victim.gender());
        assertFalse(victim.hasGermlineAberrations());
    }

    @Test
    public void testXYY() {
        MedianRatio chr1 = create("1", 1);
        MedianRatio chrX = create("X", 0.5);
        MedianRatio chrY = create("Y", TWO_Y_CUTOFF, MIN_Y_COUNT);
        List<MedianRatio> chromosomes = ImmutableList.of(chr1, chrX, chrY);
        CobaltChromosomes victim = new CobaltChromosomes(chromosomes);
        assertEquals(Gender.MALE, victim.gender());
        assertEquals(1, victim.germlineAberrations().size());
        assertTrue(victim.germlineAberrations().contains(GermlineAberration.XYY));

        assertEquals(1, victim.get("1").actualRatio(), 0.01);
        assertEquals(0.5, victim.get("X").actualRatio(), 0.01);
        assertEquals(1, victim.get("Y").actualRatio(), 0.01);
    }

    @Test
    public void testTrisomy() {
        MedianRatio chr1 = create("1", TRISOMY_CUTOFF);
        MedianRatio chrX = create("X", TRISOMY_CUTOFF);
        MedianRatio chr13 = create("13", TRISOMY_CUTOFF);
        MedianRatio chr15 = create("15", TRISOMY_CUTOFF);
        MedianRatio chr18 = create("18", TRISOMY_CUTOFF);
        MedianRatio chr21 = create("21", TRISOMY_CUTOFF);
        MedianRatio chr22 = create("22", TRISOMY_CUTOFF);
        List<MedianRatio> chromosomes = ImmutableList.of(chr1, chr13, chr15, chr18, chr21, chr22, chrX);
        CobaltChromosomes victim = new CobaltChromosomes(chromosomes);
        assertEquals(Gender.FEMALE, victim.gender());
        assertEquals(5, victim.germlineAberrations().size());
        assertTrue(victim.germlineAberrations().contains(GermlineAberration.TRISOMY_13));
        assertTrue(victim.germlineAberrations().contains(GermlineAberration.TRISOMY_15));
        assertTrue(victim.germlineAberrations().contains(GermlineAberration.TRISOMY_18));
        assertTrue(victim.germlineAberrations().contains(GermlineAberration.TRISOMY_21));
        assertTrue(victim.germlineAberrations().contains(GermlineAberration.TRISOMY_X));

        assertEquals(1.5, victim.get("13").actualRatio(), 0.01);
        assertEquals(1.5, victim.get("15").actualRatio(), 0.01);
        assertEquals(1.5, victim.get("18").actualRatio(), 0.01);
        assertEquals(1.5, victim.get("21").actualRatio(), 0.01);
        assertEquals(1, victim.get("22").actualRatio(), 0.01);
        assertEquals(1.5, victim.get("X").actualRatio(), 0.01);
    }
}
