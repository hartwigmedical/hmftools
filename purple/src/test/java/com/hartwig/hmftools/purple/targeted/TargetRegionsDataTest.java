package com.hartwig.hmftools.purple.targeted;

import static org.junit.Assert.assertEquals;

import java.nio.file.Paths;
import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.TaggedRegion;

import org.junit.Before;
import org.junit.Test;

public class TargetRegionsDataTest extends TargetRegionsTestBase
{
    @Before
    public void setup()
    {
        super.setup();
        String panelFilePath =
                Paths.get("src", "test", "resources", "panel", "panel_1.bed").toAbsolutePath().toString();
        targetRegionsData.loadTargetRegionsBed(panelFilePath, ensemblDataCache);

    }

    @Test
    public void targetRegionsTest()
    {
        List<TaggedRegion> chr1Regions = targetRegionsData.targetRegions(chr1);
        assertEquals(3, chr1Regions.size());
        assertEquals(new TaggedRegion(HumanChromosome._1.toString(), 26694021, 26694080, "ARID1A_0"), chr1Regions.get(0));

        List<TaggedRegion> chr7Regions = targetRegionsData.targetRegions(chr7);
        assertEquals(4, chr7Regions.size());
        assertEquals(new TaggedRegion(HumanChromosome._7.toString(), 55030021, 55035000, "EGFR_1"), chr7Regions.get(3));
    }
}
