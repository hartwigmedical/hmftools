package com.hartwig.hmftools.common.region;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.junit.Test;

public class TaggedRegionTest
{
    @Test
    public void loadFromBedFileTest()
    {
        String bedFilePath = Resources.getResource("bed/regions2.bed").getPath();
        Map<Chromosome, List<TaggedRegion>> map = TaggedRegion.loadRegionsFromBedFile(bedFilePath);

        assert map != null;
        assertEquals(3, map.size());
        List<TaggedRegion> chr1Regions = map.get(HumanChromosome._1);
        assertEquals(2, chr1Regions.size());
        assertEquals(new TaggedRegion("1", 11, 12, "AAA"), chr1Regions.get(0));
        assertEquals(new TaggedRegion("1", 19, 20, "BBB"), chr1Regions.get(1));
    }

    @Test
    public void loadFromBedFileWithMissingLabelsColumnTest()
    {
        String bedFilePath = Resources.getResource("bed/regions3.bed").getPath();
        Map<Chromosome, List<TaggedRegion>> map = TaggedRegion.loadRegionsFromBedFile(bedFilePath);

        assert map != null;
        assertEquals(3, map.size());
        List<TaggedRegion> chr1Regions = map.get(HumanChromosome._1);
        assertEquals(2, chr1Regions.size());
        assertEquals(new TaggedRegion("1", 11, 12, ""), chr1Regions.get(0));
        assertEquals(new TaggedRegion("1", 19, 20, ""), chr1Regions.get(1));
    }

    @Test
    public void formattedTest()
    {
        TaggedRegion taggedRegion = new TaggedRegion("1", 19, 100, "AAA");
        assertEquals("AAA:19-100", taggedRegion.formatted());
    }
}
