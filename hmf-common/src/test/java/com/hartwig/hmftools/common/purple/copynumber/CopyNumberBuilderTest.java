package com.hartwig.hmftools.common.purple.copynumber;

import static com.hartwig.hmftools.common.purple.copynumber.CopyNumberBuilder.HC_MAX_COPY_NUMBER_TOLERANCE;
import static com.hartwig.hmftools.common.purple.copynumber.CopyNumberBuilder.LC_MAX_COPY_NUMBER_TOLERANCE;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.segment.StructuralVariantSupport;

import org.junit.Test;

public class CopyNumberBuilderTest {
    private static final double EPSILON = 1e-10;

    @Test
    public void testMaxCopyNumberDeviation() {
        final FittedRegion firstWithSV = create(1, 1000, StructuralVariantSupport.MULTIPLE);
        final FittedRegion firstWithoutSV = create(1, 1000, StructuralVariantSupport.NONE);
        final FittedRegion secondWithSV = create(1001, 2000, StructuralVariantSupport.MULTIPLE);
        final FittedRegion secondWithoutSV = create(1001, 2000, StructuralVariantSupport.NONE);
        final FittedRegion thirdWithSV = create(2001, 3000, StructuralVariantSupport.MULTIPLE);
        final FittedRegion thirdWithoutSV = create(2001, 3000, StructuralVariantSupport.NONE);

        final PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.MALE, 1, 1);
        final CopyNumberBuilder builderWithSVSupport = new CopyNumberBuilder(true, purityAdjuster, secondWithSV);
        final CopyNumberBuilder builderWithoutSVSupport = new CopyNumberBuilder(true, purityAdjuster, secondWithoutSV);

        assertEquals(HC_MAX_COPY_NUMBER_TOLERANCE, builderWithSVSupport.maxCopyNumberDeviation(firstWithSV), EPSILON);
        assertEquals(HC_MAX_COPY_NUMBER_TOLERANCE, builderWithSVSupport.maxCopyNumberDeviation(firstWithoutSV), EPSILON);
        assertEquals(HC_MAX_COPY_NUMBER_TOLERANCE, builderWithSVSupport.maxCopyNumberDeviation(thirdWithSV), EPSILON);
        assertEquals(LC_MAX_COPY_NUMBER_TOLERANCE, builderWithSVSupport.maxCopyNumberDeviation(thirdWithoutSV), EPSILON);

        assertEquals(LC_MAX_COPY_NUMBER_TOLERANCE, builderWithoutSVSupport.maxCopyNumberDeviation(firstWithSV), EPSILON);
        assertEquals(LC_MAX_COPY_NUMBER_TOLERANCE, builderWithoutSVSupport.maxCopyNumberDeviation(firstWithoutSV), EPSILON);
        assertEquals(HC_MAX_COPY_NUMBER_TOLERANCE, builderWithoutSVSupport.maxCopyNumberDeviation(thirdWithSV), EPSILON);
        assertEquals(LC_MAX_COPY_NUMBER_TOLERANCE, builderWithoutSVSupport.maxCopyNumberDeviation(thirdWithoutSV), EPSILON);
    }

    private static FittedRegion create(long start, long end, StructuralVariantSupport support) {
        return PurpleDatamodelTest.createDefaultFittedRegion("1", start, end).bafCount(0).structuralVariantSupport(support).build();
    }

}
