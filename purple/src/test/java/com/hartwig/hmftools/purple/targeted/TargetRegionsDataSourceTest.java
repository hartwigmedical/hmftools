package com.hartwig.hmftools.purple.targeted;

import static com.hartwig.hmftools.common.purple.GermlineStatus.AMPLIFICATION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.DIPLOID;
import static com.hartwig.hmftools.common.purple.GermlineStatus.EXCLUDED;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.NOISE;
import static com.hartwig.hmftools.common.purple.GermlineStatus.UNKNOWN;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.hartwig.hmftools.common.genome.position.GP;
import com.hartwig.hmftools.common.purple.PurpleSegment;
import com.hartwig.hmftools.common.region.TaggedRegion;

import org.junit.Before;
import org.junit.Test;

public class TargetRegionsDataSourceTest extends TargetRegionsTestBase
{
    @Before
    public void setup()
    {
        super.setup();
    }

    @Test
    public void targetRegions()
    {
        targetRegionsData.loadTargetRegionsBed(panelFilePath("panel_1.bed"), ensemblDataCache);
        TargetRegionsDataSource dataSource = new TargetRegionsDataSource(targetRegionsData, refGenomeVersion, List.of());
        List<TaggedRegion> chr1Regions = dataSource.targetRegions(chromosome1);
        assertEquals(3, chr1Regions.size());
        assertEquals(26694020 + 1, chr1Regions.get(0).start());
        assertEquals(26694080, chr1Regions.get(0).end());
        assertEquals("ARID1A_0", chr1Regions.get(0).mTag);

        assertEquals(4, dataSource.targetRegions(chromosome7).size());
    }

    @Test
    public void germlineStatusTest()
    {
        PurpleSegment ps0 = ps("chr1", 1001, 2000, DIPLOID);
        PurpleSegment ps1 = ps("chr1", 2001, 3000, HOM_DELETION);
        PurpleSegment ps2 = ps("chr1", 3001, 4000, HET_DELETION);
        PurpleSegment ps3 = ps("chr1", 4001, 5000, NOISE);
        PurpleSegment ps4 = ps("chr2", 1001, 2000, UNKNOWN);
        PurpleSegment ps5 = ps("chr2", 2001, 3000, AMPLIFICATION);
        PurpleSegment ps6 = ps("chr2", 3001, 4000, EXCLUDED);

        List<PurpleSegment> segments = List.of(ps0, ps1, ps2, ps3, ps4, ps5, ps6);
        TargetRegionsDataSource dataSource = new TargetRegionsDataSource(targetRegionsData, refGenomeVersion, segments);

        assertEquals(DIPLOID, dataSource.germlineStatus(new GP(1500, "chr1")));
        assertEquals(DIPLOID, dataSource.germlineStatus(new GP(1600, "chr1")));
        assertNull(dataSource.germlineStatus(new GP(2500, "chr3")));
        assertEquals(AMPLIFICATION, dataSource.germlineStatus(new GP(2500, "chr2")));
        assertEquals(AMPLIFICATION, dataSource.germlineStatus(new GP(2800, "chr2")));
        assertNull(dataSource.germlineStatus(new GP(2500, "chr3")));
    }
}
