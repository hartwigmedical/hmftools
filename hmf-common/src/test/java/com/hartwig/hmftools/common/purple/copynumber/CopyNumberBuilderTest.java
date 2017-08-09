package com.hartwig.hmftools.common.purple.copynumber;

import static com.hartwig.hmftools.common.purple.copynumber.CopyNumberBuilder.HC_MAX_COPY_NUMBER_TOLERANCE;
import static com.hartwig.hmftools.common.purple.copynumber.CopyNumberBuilder.LC_MAX_COPY_NUMBER_TOLERANCE;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.segment.StructuralVariantSupport;

import org.jetbrains.annotations.NotNull;
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

    @Test
    public void testPurityAdjustedBaf() {
        testPurityAdjustedBaf(0.1);
        testPurityAdjustedBaf(0.2);
        testPurityAdjustedBaf(0.3);
        testPurityAdjustedBaf(0.4);
        testPurityAdjustedBaf(0.5);
        testPurityAdjustedBaf(0.6);
        testPurityAdjustedBaf(0.7);
        testPurityAdjustedBaf(0.8);
        testPurityAdjustedBaf(0.9);
    }

    private void testPurityAdjustedBaf(double purity) {
        testPurityAdjustedBaf(purity, 2, 1);
        testPurityAdjustedBaf(purity, 2, 2);
        testPurityAdjustedBaf(purity, 3, 2);
        testPurityAdjustedBaf(purity, 3, 3);
        testPurityAdjustedBaf(purity, 4, 2);
        testPurityAdjustedBaf(purity, 4, 3);
        testPurityAdjustedBaf(purity, 4, 4);
        testPurityAdjustedBaf(purity, 5, 3);
        testPurityAdjustedBaf(purity, 5, 4);
    }

    private void testPurityAdjustedBaf(final double purity, final int ploidy, final int alleleCount) {
        final CopyNumberBuilder builder = createBuilder(purity, 1, 100, ploidy);

        double expectedPurityAdjustedBAF = 1d * alleleCount / ploidy;
        double observedBAF = modelBAF(purity, ploidy, alleleCount);
        assertEquals(expectedPurityAdjustedBAF, builder.purityAdjustedBAF(observedBAF), EPSILON);
    }

    private CopyNumberBuilder createBuilder(double purity, long start, long end, double copyNumber) {
        return new CopyNumberBuilder(true, new PurityAdjuster(Gender.MALE, purity, 1), create(start, end, copyNumber));
    }

    private static FittedRegion create(long start, long end, double copyNumber) {
        return create("1", start, end, 0, 0, copyNumber);
    }

    private static FittedRegion create(long start, long end, int bafCount, double baf, double copyNumber) {
        return create("1", start, end, bafCount, baf, copyNumber);
    }

    private static FittedRegion create(long start, long end, StructuralVariantSupport support) {
        return PurpleDatamodelTest.createDefaultFittedRegion("1", start, end).bafCount(0).structuralVariantSupport(support).build();
    }

    @NotNull
    private static FittedRegion create(@NotNull String chromosome, long start, long end, int bafCount, double baf, double tumorCopyNumber) {
        return PurpleDatamodelTest.createDefaultFittedRegion(chromosome, start, end)
                .bafCount(bafCount)
                .observedBAF(baf)
                .tumorCopyNumber(tumorCopyNumber)
                .refNormalisedCopyNumber(tumorCopyNumber)
                .build();
    }

    private static double modelBAF(final double purity, final int ploidy, final int alleleCount) {
        assert (alleleCount >= ploidy / 2d);
        return (1 + purity * (alleleCount - 1)) / (2 + purity * (ploidy - 2));
    }
}
