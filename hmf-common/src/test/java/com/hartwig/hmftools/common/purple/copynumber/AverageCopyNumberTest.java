package com.hartwig.hmftools.common.purple.copynumber;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class AverageCopyNumberTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void excludeSexChromosomes() {
        PurpleCopyNumber cn1 = create("1", 2, 100);
        PurpleCopyNumber cn2 = create("2", 3, 100);
        PurpleCopyNumber cnX = create("X", 4, 100);
        PurpleCopyNumber cnY = create("Y", 5, 100);
        final List<PurpleCopyNumber> copyNumbers = Lists.newArrayList(cn1, cn2, cnX, cnY);

        assertEquals(2.5, AverageCopyNumber.averageCopyNumber(copyNumbers), EPSILON);
    }

    @Test
    public void excludeBafZero() {
        PurpleCopyNumber cn1 = create("1", 2, 100);
        PurpleCopyNumber cn2 = create("2", 3, 0);
        final List<PurpleCopyNumber> copyNumbers = Lists.newArrayList(cn1, cn2);

        assertEquals(2, AverageCopyNumber.averageCopyNumber(copyNumbers), EPSILON);
    }

    @Test
    public void bafWeighted() {
        PurpleCopyNumber cn1 = create("1", 2, 100);
        PurpleCopyNumber cn2 = create("2", 4, 300);
        final List<PurpleCopyNumber> copyNumbers = Lists.newArrayList(cn1, cn2);

        assertEquals(3.5, AverageCopyNumber.averageCopyNumber(copyNumbers), EPSILON);
    }

    @Test
    public void chromosomalAverage() {
        PurpleCopyNumber cn1a = create("1", 2, 100);
        PurpleCopyNumber cn1b = create("1", 3, 100);
        PurpleCopyNumber cn1c = create("1", 4, 100);
        PurpleCopyNumber cn2 = create("2", 3, 100);
        PurpleCopyNumber cnX = create("X", 4, 100);
        PurpleCopyNumber cnY = create("Y", 5, 100);
        final List<PurpleCopyNumber> copyNumbers = Lists.newArrayList(cn1a, cn1b, cn1c, cn2, cnX, cnY);
        final Map<String, Double> averagePerChromosome = AverageCopyNumber.averageChromosomalCopyNumber(copyNumbers);
        assertEquals(3, averagePerChromosome.get("1"), EPSILON);
        assertEquals(3, averagePerChromosome.get("2"), EPSILON);
        assertEquals(4, averagePerChromosome.get("X"), EPSILON);
        assertEquals(5, averagePerChromosome.get("Y"), EPSILON);
    }

    @NotNull
    private PurpleCopyNumber create(@NotNull final String chromosome, double copyNumber, int bafCount) {
        return PurpleDatamodelTest.createCopyNumber(chromosome, 1, 100, copyNumber).bafCount(bafCount).build();
    }

}
