package com.hartwig.hmftools.common.genome.chromosome;

import static com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes.MIN_RATIO_COUNT;
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
import com.hartwig.hmftools.common.cobalt.ImmutableMedianRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.compress.utils.Lists;
import org.junit.Test;

public class CobaltChromosomesTest {

    public static CobaltChromosomes female() {
        List<MedianRatio> ratios = Lists.newArrayList();
        for (int i = 0; i < 22; i++) {
            ratios.add(create(String.valueOf(i), 1));
        }

        ratios.add(create("X", 1));
        return new CobaltChromosomes(ratios);
    }

    public static CobaltChromosomes male() {
        List<MedianRatio> ratios = Lists.newArrayList();
        for (int i = 0; i < 22; i++) {
            ratios.add(create(String.valueOf(i), 1));
        }

        ratios.add(create("X", 0.5));
        ratios.add(create("Y", 0.5));
        return new CobaltChromosomes(ratios);
    }

    @Test
    public void testFemale() {
        MedianRatio chrX = create("X", 1, MIN_RATIO_COUNT);
        MedianRatio chrY = create("Y", 0.05, 100);
        List<MedianRatio> chromosomes = ImmutableList.of(chrX, chrY);
        CobaltChromosomes victim = new CobaltChromosomes(chromosomes);
        assertEquals(Gender.FEMALE, victim.gender());
        assertEquals(0, victim.aberrations().size());

        assertEquals(1, victim.get("X").actualRatio(), 0.01);
        assertFalse(victim.contains("Y"));
    }

    @Test
    public void testMale() {
        MedianRatio chrX = create("X", 0.5, MIN_RATIO_COUNT);
        MedianRatio chrY = create("Y", Y_CUTOFF, MIN_RATIO_COUNT);
        List<MedianRatio> chromosomes = ImmutableList.of(chrX, chrY);
        CobaltChromosomes victim = new CobaltChromosomes(chromosomes);
        assertEquals(Gender.MALE, victim.gender());
        assertEquals(0, victim.aberrations().size());

        assertEquals(0.5, victim.get("X").actualRatio(), 0.01);
        assertEquals(0.5, victim.get("Y").actualRatio(), 0.01);
    }

    @Test
    public void testKleinfelter() {
        MedianRatio chrX = create("X", TWO_X_CUTOFF, MIN_RATIO_COUNT);
        MedianRatio chrY = create("Y", Y_CUTOFF, MIN_RATIO_COUNT);
        List<MedianRatio> chromosomes = ImmutableList.of(chrX, chrY);
        CobaltChromosomes victim = new CobaltChromosomes(chromosomes);
        assertEquals(Gender.MALE, victim.gender());
        assertEquals(1, victim.aberrations().size());
        assertTrue(victim.aberrations().contains(ChromosomalAberration.KLINEFELTER));

        assertEquals(1, victim.get("X").actualRatio(), 0.01);
        assertEquals(0.5, victim.get("Y").actualRatio(), 0.01);
    }

    @Test
    public void testFemaleEvenWithSmallYCount() {
        MedianRatio chrX = create("X", 1, MIN_RATIO_COUNT);
        MedianRatio chrY = create("Y", Y_CUTOFF, MIN_RATIO_COUNT - 1);
        List<MedianRatio> chromosomes = ImmutableList.of(chrX, chrY);
        CobaltChromosomes victim = new CobaltChromosomes(chromosomes);
        assertEquals(Gender.FEMALE, victim.gender());
        assertEquals(0, victim.aberrations().size());
    }

    @Test
    public void testFemaleEvenWithSmallYMedian() {
        MedianRatio chrX = create("X", 1, MIN_RATIO_COUNT);
        MedianRatio chrY = create("Y", Y_CUTOFF - Y_CUTOFF/10d, MIN_RATIO_COUNT);
        List<MedianRatio> chromosomes = ImmutableList.of(chrX, chrY);
        CobaltChromosomes victim = new CobaltChromosomes(chromosomes);
        assertEquals(Gender.FEMALE, victim.gender());
        assertEquals(0, victim.aberrations().size());
    }

    @Test
    public void testMosiacX() {
        assertTrue(Doubles.lessThan(TWO_X_CUTOFF, MOSIAC_X_CUTOFF));
        MedianRatio chr1 = create("1", 1, MIN_RATIO_COUNT);
        MedianRatio chrX = create("X", TWO_X_CUTOFF, MIN_RATIO_COUNT);
        List<MedianRatio> chromosomes = ImmutableList.of(chr1, chrX);
        CobaltChromosomes victim = new CobaltChromosomes(chromosomes);
        assertEquals(Gender.FEMALE, victim.gender());
        assertEquals(1, victim.aberrations().size());
        assertTrue(victim.aberrations().contains(ChromosomalAberration.MOSAIC_X));

        assertEquals(1, victim.get("1").actualRatio(), 0.01);
        assertEquals(TWO_X_CUTOFF, victim.get("X").actualRatio(), 0.01);
    }

    @Test
    public void testGCBiasNotMosiacX() {
        assertTrue(Doubles.lessThan(TWO_X_CUTOFF, MOSIAC_X_CUTOFF));
        MedianRatio chr1 = create("1", TWO_X_CUTOFF, MIN_RATIO_COUNT);
        MedianRatio chrX = create("X", TWO_X_CUTOFF, MIN_RATIO_COUNT);
        List<MedianRatio> chromosomes = ImmutableList.of(chr1, chrX);
        CobaltChromosomes victim = new CobaltChromosomes(chromosomes);
        assertEquals(Gender.FEMALE, victim.gender());
        assertEquals(0, victim.aberrations().size());
    }

    @Test
    public void testXYY() {
        MedianRatio chr1 = create("1", 1);
        MedianRatio chrX = create("X", 0.5);
        MedianRatio chrY = create("Y", TWO_Y_CUTOFF);
        List<MedianRatio> chromosomes = ImmutableList.of(chr1, chrX, chrY);
        CobaltChromosomes victim = new CobaltChromosomes(chromosomes);
        assertEquals(Gender.MALE, victim.gender());
        assertEquals(1, victim.aberrations().size());
        assertTrue(victim.aberrations().contains(ChromosomalAberration.XYY));

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
        assertEquals(5, victim.aberrations().size());
        assertTrue(victim.aberrations().contains(ChromosomalAberration.TRISOMY_13));
        assertTrue(victim.aberrations().contains(ChromosomalAberration.TRISOMY_15));
        assertTrue(victim.aberrations().contains(ChromosomalAberration.TRISOMY_18));
        assertTrue(victim.aberrations().contains(ChromosomalAberration.TRISOMY_21));
        assertTrue(victim.aberrations().contains(ChromosomalAberration.TRISOMY_X));

        assertEquals(1.5, victim.get("13").actualRatio(), 0.01);
        assertEquals(1.5, victim.get("15").actualRatio(), 0.01);
        assertEquals(1.5, victim.get("18").actualRatio(), 0.01);
        assertEquals(1.5, victim.get("21").actualRatio(), 0.01);
        assertEquals(1, victim.get("22").actualRatio(), 0.01);
        assertEquals(1.5, victim.get("X").actualRatio(), 0.01);
    }


    private static MedianRatio create(String contig, double ratio) {
        return create(contig, ratio, MIN_RATIO_COUNT);
    }

    private static MedianRatio create(String contig, double ratio, int count) {
        return ImmutableMedianRatio.builder().count(count).chromosome(contig).medianRatio(ratio).build();
    }

}
