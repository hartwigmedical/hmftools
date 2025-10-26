package com.hartwig.hmftools.cobalt.diploid;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;

import static junit.framework.TestCase.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.google.common.io.Files;
import com.hartwig.hmftools.cobalt.e2e.DiploidFileSection;
import com.hartwig.hmftools.cobalt.testutils.DiploidRegionsFileWriter;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

import org.junit.Test;

public class DiploidRegionLoaderTest
{
    @Test
    public void regionsTest() throws IOException
    {
        File tempDir = Files.createTempDir();
        File diploidBedFile = new File(tempDir, "diploid_bed.tsv.gz");
        DiploidRegionsFileWriter diploidRegionsFileWriter = new DiploidRegionsFileWriter();
        diploidRegionsFileWriter.addSection(new DiploidFileSection(_1, 1000, 3000));
        diploidRegionsFileWriter.addSection(new DiploidFileSection(_1, 5000, 6000));
        diploidRegionsFileWriter.addSection(new DiploidFileSection(_2, 1000, 3000));
        diploidRegionsFileWriter.addSection(new DiploidFileSection(_2, 7000, 13000));
        diploidRegionsFileWriter.write(diploidBedFile);

        DiploidRegionLoader diploidRegionLoader = new DiploidRegionLoader( diploidBedFile.getAbsolutePath());

        ListMultimap<Chromosome, DiploidStatus> regions = diploidRegionLoader.regions();
        assertEquals(19, regions.size());
        List<DiploidStatus> regions1 = regions.get(_1);
        assertEquals(6, regions1.size());
        int i=0;
        check(regions1.get(i++), _1, 1, false);
        check(regions1.get(i++), _1, 1001, true);
        check(regions1.get(i++), _1, 2001, true);
        check(regions1.get(i++), _1, 3001, false);
        check(regions1.get(i++), _1, 4001, false);
        check(regions1.get(i), _1, 5001, true);

        List<DiploidStatus> regions2 = regions.get(_2);
        assertEquals(13, regions2.size());
        i=0;
        check(regions2.get(i++), _2, 1, false);
        check(regions2.get(i++), _2, 1001, true);
        check(regions2.get(i++), _2, 2001, true);
        check(regions2.get(i++), _2, 3001, false);
        check(regions2.get(i++), _2, 4001, false);
        check(regions2.get(i++), _2, 5001, false);
        check(regions2.get(i++), _2, 6001, false);
        check(regions2.get(i++), _2, 7001, true);
        check(regions2.get(i), _2, 8001, true);
        check(regions2.get(12), _2, 12001, true);
    }

    private void check(DiploidStatus status, Chromosome chromosome, int position, boolean isDiploid)
    {
        assertEquals(chromosome, status.humanChromosome());
        assertEquals(position, status.start());
        assertEquals(isDiploid, status.isDiploid);
    }
}