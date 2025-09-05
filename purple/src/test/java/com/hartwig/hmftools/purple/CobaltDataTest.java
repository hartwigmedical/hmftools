package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._X;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.pcf.GenomeIntervals;
import com.hartwig.hmftools.common.utils.pcf.PCFFile;
import com.hartwig.hmftools.common.utils.pcf.PCFIntervals;
import com.hartwig.hmftools.purple.exclusions.ExcludedRegionsFileTest;

import org.apache.commons.io.FileUtils;
import org.junit.Before;
import org.junit.Test;

public class CobaltDataTest
{
    private static final double DELTA = 0.0001;
    File tempDir = FileUtils.getTempDirectory();
    private final String sample = "blah";

    @Before
    public void setup() throws Exception
    {
        // Write the ratios file.
        List<CobaltRatio> toStore = new ArrayList<>();
        addRatios(toStore, "1", 10);
        addRatios(toStore, "2", 20);
        addRatios(toStore, "X", 5);
        String regionsFile = CobaltRatioFile.generateFilename(tempDir.getAbsolutePath(), sample);
        CobaltRatioFile.write(regionsFile, toStore);

        // We need a pcf file too.
        List<ChrBaseRegion> regions1 = new ArrayList<>();
        regions1.add(region(_1, 1001, 2000));
        regions1.add(region(_1, 2001, 3000));
        PCFIntervals intervals1 = new PCFIntervals(_1, regions1);

        List<ChrBaseRegion> regions2 = new ArrayList<>();
        regions2.add(region(_2, 1101, 1200));
        regions2.add(region(_2, 1201, 1300));
        PCFIntervals intervals2 = new PCFIntervals(_2, regions2);

        GenomeIntervals data = new GenomeIntervals(List.of(intervals1, intervals2));
        String filename = PCFFile.generateRatioFilename(tempDir.getAbsolutePath(), sample);
        PCFFile.write(filename, data);
    }

    @Test
    public void noExcludedRegions() throws Exception
    {
        CobaltData cobaltData = new CobaltData(null, sample, tempDir.getAbsolutePath(), Gender.FEMALE, true, false, null);
        assertEquals(3, cobaltData.Ratios.size());
        // The ratio data for chromosomes 1 and 2 should be precisely what was written to the ratios file.
        List<CobaltRatio> ratios1 = cobaltData.Ratios.get(_1);
        assertEquals(10, ratios1.size());
        for(int i = 0; i < ratios1.size(); i++)
        {
            assertEquals(cr("1", i), ratios1.get(i));
        }
        List<CobaltRatio> ratios2 = cobaltData.Ratios.get(_2);
        assertEquals(20, ratios2.size());
        for(int i = 0; i < ratios2.size(); i++)
        {
            assertEquals(cr("2", i), ratios2.get(i));
        }
        // For chromosome X the ratios are adjusted for the supplied gender (female).
        List<CobaltRatio> ratiosX = cobaltData.Ratios.get(_X);
        assertEquals(5, ratiosX.size());
        for(int i = 0; i < ratiosX.size(); i++)
        {
            assertEquals("X", ratiosX.get(i).chromosome());
            assertEquals(i * 1000L + 1, ratiosX.get(i).position());
            assertEquals(i * 0.01, ratiosX.get(i).referenceReadDepth(), DELTA);
            assertEquals(1.0, ratiosX.get(i).referenceGCRatio(), DELTA);
        }
    }

    @Test
    public void ratiosAreFilteredByExcludedRegions() throws Exception
    {
        File regionsFile = new File(tempDir, "regions.tsv");
        List<String> lines = new ArrayList<>();
        lines.add(ExcludedRegionsFileTest.line("1", 1200, 1320, "BLAH1"));
        lines.add(ExcludedRegionsFileTest.line("1", 3900, 4020, "BLAH2"));
        lines.add(ExcludedRegionsFileTest.line("1", 9900, 10_020, "BLAH3"));
        lines.add(ExcludedRegionsFileTest.line("1", 39900, 40_020, "BLAH4"));
        lines.add(ExcludedRegionsFileTest.line("2", 5400, 5520, "BLAH5"));

        ExcludedRegionsFileTest.writeExcludedRegionsFile(regionsFile, lines);
        CobaltData cobaltData =
                new CobaltData(null, sample, tempDir.getAbsolutePath(), Gender.FEMALE, true, false, regionsFile.getAbsolutePath());
        assertEquals(3, cobaltData.Ratios.size());

        List<CobaltRatio> ratios1 = cobaltData.Ratios.get(_1);
        assertEquals(10, ratios1.size());
        assertEquals(cr("1", 0), ratios1.get(0));
        assertEquals(cr("1", 1).mask(), ratios1.get(1)); // 1001-2000, intersects BLAH1
        assertEquals(cr("1", 2), ratios1.get(2));
        assertEquals(cr("1", 3).mask(), ratios1.get(3)); // 3001-4000, intersects BLAH2
        assertEquals(cr("1", 4).mask(), ratios1.get(4)); // 4001-5000, intersects BLAH2
        assertEquals(cr("1", 5), ratios1.get(5));
        assertEquals(cr("1", 6), ratios1.get(6));
        assertEquals(cr("1", 7), ratios1.get(7));
        assertEquals(cr("1", 8), ratios1.get(8));
        assertEquals(cr("1", 9).mask(), ratios1.get(9)); // 9001-10000, intersects BLAH3

        List<CobaltRatio> ratios2 = cobaltData.Ratios.get(_2);
        assertEquals(20, ratios2.size());

        // Only one of the cobalt ratios gets masked.
        for(int i = 0; i < 5; i++)
        {
            assertEquals(cr("2", i), ratios2.get(i));
        }
        assertEquals(cr("2", 5).mask(), ratios2.get(5));
        for(int i = 6; i < ratios2.size(); i++)
        {
            assertEquals(cr("2", i), ratios2.get(i));
        }
    }

    private void addRatios(List<CobaltRatio> list, String chromosome, int count)
    {
        for(int i = 0; i < count; i++)
        {
            list.add(cr(chromosome, i));
        }
    }

    private CobaltRatio cr(String chromosome, int index)
    {
        int position = index * 1000 + 1;
        double d = index * 0.01;
        return new CobaltRatio(chromosome, position, d, d, d, d, d, d, d);
    }

    ChrBaseRegion region(HumanChromosome chromosome, int start, int end)
    {
        return new ChrBaseRegion(V38.versionedChromosome(chromosome), start, end);
    }
}
