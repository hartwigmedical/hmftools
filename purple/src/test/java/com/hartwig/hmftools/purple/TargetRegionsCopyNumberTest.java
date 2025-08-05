package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.purple.CopyNumberMethod.GERMLINE_HET2HOM_DELETION;

import static org.immutables.value.internal.$guava$.collect.$ImmutableList.of;
import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;
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
    public void headingsTest()
    {
        String expected = "chromosome\twindowStart\twindowEnd\tbedRegions\tmasked\taverageDepth\twindowGCContent\twindowTumorRatio\tregionStart"
                + "\tregionEnd\tcopyNumber\tminorAlleleCopyNumber\tdepthWindowCount\tbafCount\tGCContent\tCNMethod";
        assertEquals(expected, TargetRegionsCopyNumber.tsvFileHeader());
    }

    @Test
    public void tsvTest()
    {
        CobaltRatio cobaltRatio = new CobaltRatio("chr1", 55_000_001, -1, -1, -1, 1.011111111122222, 183.845111111111, 0.551111111, 1.0122222222222);
        TaggedRegion region = new TaggedRegion(cobaltRatio.chromosome(), 55_000_080, 55_000_130, "BLAH");
        PurpleCopyNumber purpleCopyNumber= ImmutablePurpleCopyNumber.builder()
                .chromosome(cobaltRatio.chromosome())
                .start(1000)
                .end(100_000_000)
                .averageTumorCopyNumber(12.11111111111111)
                .segmentStartSupport(SegmentSupport.TELOMERE)
                .segmentEndSupport(SegmentSupport.NONE)
                .method(GERMLINE_HET2HOM_DELETION)
                .bafCount(123)
                .depthWindowCount(2)
                .gcContent(0.550202020200202)
                .minStart(0)
                .maxStart(150_000_000)
                .averageObservedBAF(0.5611111111111111111)
                .averageActualBAF(0.57333333333333).build();
        double minorAlleleCopyNumber = purpleCopyNumber.minorAlleleCopyNumber();
        assertEquals(5.1674, minorAlleleCopyNumber, 0.00001); // see below
        String method = GERMLINE_HET2HOM_DELETION.name();
        TargetRegionsCopyNumber trc = new TargetRegionsCopyNumber(cobaltRatio, of(region), purpleCopyNumber);
        final String cobaltPart = "chr1\t55000001\t55001000\tBLAH:55000080-55000130\tfalse\t183.8451\t1.0122\t0.5511";
        final String purplePart = "\t1000\t100000000\t12.1111\t5.1674\t2\t123\t0.5502\t" + method;
        assertEquals(cobaltPart+purplePart, trc.toTSV());
    }
}
