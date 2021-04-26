package com.hartwig.hmftools.common.utils.pcf;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PCFFileTest {

    private static final int WINDOW = 1000;
    private static final String BASE_PATH = Resources.getResource("pcf").getPath();

    @Test
    public void testBafFile() throws IOException {
        final ListMultimap<Chromosome, PCFPosition> resultMap =
                PCFFile.readPositions(WINDOW, PCFSource.TUMOR_BAF, BASE_PATH + File.separator + "baf.pcf");

        final List<PCFPosition> chromOneResults = resultMap.get(HumanChromosome._1);

        assertEquals(4, chromOneResults.size());
        assertPosition(93548001, chromOneResults.get(0));
        assertPosition(193800001, chromOneResults.get(1));
        assertPosition(193803001, chromOneResults.get(2));
        assertPosition(193804001, chromOneResults.get(3));

        final List<PCFPosition> chromTwoResults = resultMap.get(HumanChromosome._2);

        assertEquals(4, chromTwoResults.size());
        assertPosition(93548001, chromTwoResults.get(0));
        assertPosition(193800001, chromTwoResults.get(1));
        assertPosition(193803001, chromTwoResults.get(2));
        assertPosition(193804001, chromTwoResults.get(3));
    }

    @Test
    public void testRatioFile() throws IOException {
        final ListMultimap<Chromosome, PCFPosition> resultMap =
                PCFFile.readPositions(WINDOW, PCFSource.TUMOR_BAF, BASE_PATH + File.separator + "ratio.pcf");

        final List<PCFPosition> chromosomeOneResults = resultMap.get(HumanChromosome._1);

        assertEquals(5, chromosomeOneResults.size());
        assertPosition(835001, 1, 835001, chromosomeOneResults.get(0));
        assertPosition(2583001, 2583001, 2583001, chromosomeOneResults.get(1));
        assertPosition(2584001, 2584001, 2695001, chromosomeOneResults.get(2));
        assertPosition(2695001, 2584001, 2695001, chromosomeOneResults.get(3));
        assertPosition(4363001, 4363001, 4363001, chromosomeOneResults.get(4));

        final List<PCFPosition> chromosomeThreeResults = resultMap.get(HumanChromosome._3);
        assertEquals(5, chromosomeThreeResults.size());
        assertPosition(90449001, 1, 90449001, chromosomeThreeResults.get(0));
        assertPosition(90452001, 90452001, 90452001, chromosomeThreeResults.get(1));
        assertPosition(90455001, 90455001, 90457001, chromosomeThreeResults.get(2));
        assertPosition(90457001, 90455001, 90457001, chromosomeThreeResults.get(3));
        assertPosition(90458001, 90458001, 90458001, chromosomeThreeResults.get(4));
    }

    private static void assertPosition(long position, long min, long max, @NotNull final PCFPosition victim) {
        assertEquals(position, victim.position());
        assertEquals(min, victim.minPosition());
        assertEquals(max, victim.maxPosition());
    }

    private static void assertPosition(long position, @NotNull final PCFPosition victim) {
        assertEquals(position, victim.position());
    }
}
