package com.hartwig.hmftools.common.purple.copynumber.tolerance;

import static com.hartwig.hmftools.common.purple.copynumber.tolerance.CopyNumberDeviation.HC_MAX_COPY_NUMBER_TOLERANCE;
import static com.hartwig.hmftools.common.purple.copynumber.tolerance.CopyNumberDeviation.LC_MAX_COPY_NUMBER_TOLERANCE;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.copynumber.tolerance.CopyNumberDeviation;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CopyNumberDeviationTest {
    private static final double EPSILON = 1e-10;

    @Test
    public void testMaxCopyNumberDeviation() {
        final FittedRegion firstWithSV = create(1, 1000, SegmentSupport.MULTIPLE);
        final FittedRegion firstWithoutSV = create(1, 1000, SegmentSupport.NONE);
        final FittedRegion secondWithSV = create(1001, 2000, SegmentSupport.MULTIPLE);
        final FittedRegion secondWithoutSV = create(1001, 2000, SegmentSupport.NONE);
        final FittedRegion thirdWithSV = create(2001, 3000, SegmentSupport.MULTIPLE);
        final FittedRegion thirdWithoutSV = create(2001, 3000, SegmentSupport.NONE);

        assertEquals(HC_MAX_COPY_NUMBER_TOLERANCE, CopyNumberDeviation.maxCopyNumberDeviation(secondWithSV, firstWithSV), EPSILON);
        assertEquals(HC_MAX_COPY_NUMBER_TOLERANCE, CopyNumberDeviation.maxCopyNumberDeviation(secondWithSV, firstWithoutSV), EPSILON);
        assertEquals(HC_MAX_COPY_NUMBER_TOLERANCE, CopyNumberDeviation.maxCopyNumberDeviation(secondWithSV, thirdWithSV), EPSILON);
        assertEquals(LC_MAX_COPY_NUMBER_TOLERANCE, CopyNumberDeviation.maxCopyNumberDeviation(secondWithSV, thirdWithoutSV), EPSILON);

        assertEquals(LC_MAX_COPY_NUMBER_TOLERANCE, CopyNumberDeviation.maxCopyNumberDeviation(secondWithoutSV, firstWithSV), EPSILON);
        assertEquals(LC_MAX_COPY_NUMBER_TOLERANCE, CopyNumberDeviation.maxCopyNumberDeviation(secondWithoutSV, firstWithoutSV), EPSILON);
        assertEquals(HC_MAX_COPY_NUMBER_TOLERANCE, CopyNumberDeviation.maxCopyNumberDeviation(secondWithoutSV, thirdWithSV), EPSILON);
        assertEquals(LC_MAX_COPY_NUMBER_TOLERANCE, CopyNumberDeviation.maxCopyNumberDeviation(secondWithoutSV, thirdWithoutSV), EPSILON);
    }

    @NotNull
    private static FittedRegion create(long start, long end, SegmentSupport support) {
        return PurpleDatamodelTest.createDefaultFittedRegion("1", start, end).bafCount(0).support(support).build();
    }
}
