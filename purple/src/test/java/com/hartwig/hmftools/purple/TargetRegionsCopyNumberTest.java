package com.hartwig.hmftools.purple;

import static org.immutables.value.internal.$guava$.collect.$ImmutableList.of;
import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.SegmentSupport;
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
        CobaltRatio cobaltRatio = new CobaltRatio("chr1", 55_000_001, -1, -1, -1, 1.0, 183.845, 0.55, 1.01);
        TaggedRegion region = new TaggedRegion(55_000_080, 55_000_130, "BLAH");
        PurpleCopyNumber purpleCopyNumber= ImmutablePurpleCopyNumber.builder()
                .chromosome(cobaltRatio.chromosome())
                .start(1000)
                .end(100_000_000)
                .averageTumorCopyNumber(12.1)
                .segmentStartSupport(SegmentSupport.TELOMERE)
                .segmentEndSupport(SegmentSupport.NONE)
                .method(CopyNumberMethod.UNKNOWN)
                .bafCount(123)
                .depthWindowCount(2)
                .gcContent(0.55)
                .minStart(0)
                .maxStart(150_000_000)
                .averageObservedBAF(0.56)
                .averageActualBAF(0.57).build();
        double minorAlleleCopyNumber = purpleCopyNumber.minorAlleleCopyNumber();
        TargetRegionsCopyNumber trc = new TargetRegionsCopyNumber(cobaltRatio, of(region), purpleCopyNumber);
        final String cobaltPart = "chr1\t55000001\t55001000\tBLAH:55000080-55000130\tfalse\t183.845\t1.01\t0.55";
        final String purplePart = "\t1000\t100000000\t12.1\t" + minorAlleleCopyNumber + "\t2\t123\t0.55";
        assertEquals(cobaltPart+purplePart, trc.toTSV());
    }
}
