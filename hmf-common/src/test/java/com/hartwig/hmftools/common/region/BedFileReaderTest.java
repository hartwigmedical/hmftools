package com.hartwig.hmftools.common.region;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.junit.Test;

public class BedFileReaderTest
{
    @Test
    public void readBedFileTest() throws Exception
    {
        String bedFilePath = Resources.getResource("bed/regions1.bed").getPath();
        List<ChrBaseRegion> regions = BedFileReader.loadBedFile(bedFilePath);

        assertEquals(3, regions.size());
        // Note that the start is 1+ what is in the file.
        assertEquals(new ChrBaseRegion("16", 105265, 105385), regions.get(0));
        assertEquals(new ChrBaseRegion("16", 1849806, 1849926), regions.get(1));
    }

    @Test
    public void readChrBedFileTest()
    {
        String bedFilePath = Resources.getResource("bed/regions2.bed").getPath();
        Map<Chromosome, List<BaseRegion>> map = BedFileReader.loadBedFileChrMap(bedFilePath);

        assertEquals(3, map.size());
        List<BaseRegion> chr1Regions = map.get(HumanChromosome._1);
        assertEquals(2, chr1Regions.size());
        assertEquals(new BaseRegion(11, 12), chr1Regions.get(0));
        assertEquals(new BaseRegion(19, 20), chr1Regions.get(1));

        List<BaseRegion> chrXRegions = map.get(HumanChromosome._1);
        assertEquals(2, chrXRegions.size());
        assertEquals(new BaseRegion(11, 12), chrXRegions.get(0));
        assertEquals(new BaseRegion(19, 20), chrXRegions.get(1));
    }
}
