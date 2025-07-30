package com.hartwig.hmftools.purple;

import static org.immutables.value.internal.$guava$.collect.$ImmutableList.of;
import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.region.TaggedRegion;

import org.junit.Before;
import org.junit.Test;

public class TargetRegionsCopyNumberTest extends TargetRegionsTestBase
{
    @Before
    public void setup()
    {
        super.setup();
    }

    @Test
    public void tsvTest()
    {
        CobaltRatio cobaltRatio = new CobaltRatio("chr1", 55_000_000, -1, -1, -1, 1.0, 183.845, 0.55, 1.01);
        TaggedRegion region = new TaggedRegion(55_000_080, 55_000_130, "BLAH");
        TargetRegionsCopyNumber trc = new TargetRegionsCopyNumber(cobaltRatio, of(region));
        assertEquals("chr1\t55000000\t55001000\tBLAH:55000080-55000130\tfalse\t183.845\t1.01\t0.55", trc.toTSV());
    }
}
