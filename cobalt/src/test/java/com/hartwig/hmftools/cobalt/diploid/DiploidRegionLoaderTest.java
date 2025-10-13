package com.hartwig.hmftools.cobalt.diploid;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;

import static junit.framework.TestCase.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.io.Files;
import com.hartwig.hmftools.cobalt.ChromosomePositionCodec;
import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.cobalt.e2e.DiploidFileSection;
import com.hartwig.hmftools.cobalt.testutils.DiploidRegionsFileWriter;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

import htsjdk.samtools.util.Locatable;

import tech.tablesaw.api.*;

public class DiploidRegionLoaderTest
{
    @Test
    public void testBuildRatios() throws IOException
    {
        File tempDir = Files.createTempDir();
        File diploidBedFile = new File(tempDir, "diploid_bed.tsv.gz");
        DiploidRegionsFileWriter diploidRegionsFileWriter = new DiploidRegionsFileWriter();
        diploidRegionsFileWriter.addSection(new DiploidFileSection(_1, 1000, 3000));
        diploidRegionsFileWriter.addSection(new DiploidFileSection(_1, 5000, 6000));
        diploidRegionsFileWriter.addSection(new DiploidFileSection(_2, 1000, 3000));
        diploidRegionsFileWriter.addSection(new DiploidFileSection(_2, 7000, 13000));
        diploidRegionsFileWriter.write(diploidBedFile);
        String chr1 = "1";
        String chr2 = "2";
        ChromosomePositionCodec chromosomePosCodec = new ChromosomePositionCodec();
        DiploidRegionLoader diploidRegionLoader = null;

        try
        {
            diploidRegionLoader = new DiploidRegionLoader(chromosomePosCodec, diploidBedFile.getAbsolutePath());
        }
        catch(Exception e)
        {

        }

        Table result = diploidRegionLoader.build();
        assertEquals(11, result.rowCount());
        StringColumn chrColumn = result.stringColumn(CobaltColumns.CHROMOSOME);
        assertReadRatio("1", 1001, result.where(chrColumn.isEqualTo(chr1)).row(0));
        assertReadRatio("1", 2001, result.where(chrColumn.isEqualTo(chr1)).row(1));
        assertReadRatio("1", 5001, result.where(chrColumn.isEqualTo(chr1)).row(2));
        assertReadRatio("2", 1001, result.where(chrColumn.isEqualTo(chr2)).row(0));
        assertReadRatio("2", 2001, result.where(chrColumn.isEqualTo(chr2)).row(1));
    }

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

        ChromosomePositionCodec chromosomePosCodec = new ChromosomePositionCodec();
        DiploidRegionLoader diploidRegionLoader = null;
        diploidRegionLoader = new DiploidRegionLoader(chromosomePosCodec, diploidBedFile.getAbsolutePath());

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

    private void assertReadRatio(String contig, long position, Row victim)
    {
        assertEquals(contig, victim.getString(CobaltColumns.CHROMOSOME));
        assertEquals(position, victim.getInt(CobaltColumns.POSITION));
    }

    private static Locatable locatable(String contig, int start, int end)
    {
        return new Locatable()
        {
            @Override
            public String getContig()
            {
                return contig;
            }

            @Override
            public int getStart()
            {
                return start;
            }

            @Override
            public int getEnd()
            {
                return end;
            }
        };
    }
}
