package com.hartwig.hmftools.common.utils.pcf;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;
import static com.hartwig.hmftools.common.utils.pcf.PCFSource.TUMOR_RATIO;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.commons.io.FileUtils;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PCFFileTest extends PCFTestBase
{
    private static final int WINDOW = 1000;
    private static final String BASE_PATH = Resources.getResource("pcf").getPath();

    @Test
    public void writeTest() throws Exception
    {
        List<ChrBaseRegion> regions1 = new ArrayList<>();
        regions1.add(region(_1, 101, 200));
        regions1.add(region(_1, 201, 300));
        regions1.add(region(_1, 401, 500));
        PCFIntervals intervals1 = new PCFIntervals(_1, regions1);

        List<ChrBaseRegion> regions2 = new ArrayList<>();
        regions2.add(region(_2, 1101, 1200));
        regions2.add(region(_2, 1201, 1300));
        regions2.add(region(_2, 1401, 1500));
        PCFIntervals intervals2 = new PCFIntervals(_2, regions2);

        List<ChrBaseRegion> regions3 = new ArrayList<>();
        PCFIntervals intervals3 = new PCFIntervals(_3, regions3);

        GenomeIntervals data = new GenomeIntervals(List.of(intervals1, intervals2, intervals3));

        File tempDir = FileUtils.getTempDirectory();
        String tempFile = new File(tempDir, "whatever.pcf").getAbsolutePath();
        PCFFile.write(tempFile, data);

        ListMultimap<Chromosome, PCFPosition> read = PCFFile.readPositions(100, TUMOR_RATIO, tempFile);
        assertEquals(2, read.asMap().size());
        List<PCFPosition> positions1 = read.get(HumanChromosome._1);
        assertEquals(5, positions1.size());
        assertPosition(101, 1, 101, positions1.get(0));
        assertPosition(201, 201, 201, positions1.get(1));
        assertPosition(301, 301, 401, positions1.get(2));
        assertPosition(401, 301, 401, positions1.get(3));
        assertPosition(501, 501, 501, positions1.get(4));
    }

    @Test
    public void testBafFile() throws IOException
    {
        final ListMultimap<Chromosome, PCFPosition> resultMap =
                PCFFile.readPositions(WINDOW, PCFSource.TUMOR_BAF, BASE_PATH + File.separator + "baf.pcf");

        final List<PCFPosition> chromOneResults = resultMap.get(_1);

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
    public void testRatioFile() throws IOException
    {
        final ListMultimap<Chromosome, PCFPosition> resultMap =
                PCFFile.readPositions(WINDOW, PCFSource.TUMOR_BAF, BASE_PATH + File.separator + "ratio.pcf");

        final List<PCFPosition> chromosomeOneResults = resultMap.get(_1);

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

    private static void assertPosition(long position, long min, long max, @NotNull final PCFPosition victim)
    {
        assertEquals(position, victim.position());
        assertEquals(min, victim.minPosition());
        assertEquals(max, victim.maxPosition());
    }

    private static void assertPosition(long position, @NotNull final PCFPosition victim)
    {
        assertEquals(position, victim.position());
    }
}
