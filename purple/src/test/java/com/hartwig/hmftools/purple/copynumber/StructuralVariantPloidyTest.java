package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.purple.FittingTestUtils.buildPurityAdjuster;

import static org.junit.Assert.assertEquals;

import java.util.Optional;

import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegPloidy;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

@SuppressWarnings("OptionalUsedAsFieldOrParameterType")
public class StructuralVariantPloidyTest
{
    static final String CHROMOSOME = "X";
    static final PurityAdjuster PURE = buildPurityAdjuster(Gender.FEMALE, 1d, 1d);

    private static final double EPSILON = 1e-10;
    private static final double PLOIDY = 1;

    @Test
    public void testCompleteImpliedCopyNumbers()
    {
        final StructuralVariantLegPloidy positive = create(1, Optional.of(4d), Optional.of(3d));
        assertLeft(4d, 1d, positive);
        assertRight(3d, 1d, positive);

        final StructuralVariantLegPloidy negative = create(-1, Optional.of(3d), Optional.of(4d));
        assertLeft(3d, 1d, negative);
        assertRight(4d, 1d, negative);
    }

    @Test
    public void testImpliedCopyNumberMissingLeft()
    {
        final StructuralVariantLegPloidy positive = create(1, Optional.empty(), Optional.of(3d));
        assertLeft(4d, 1d, positive);
        assertRight(0d, 0d, positive);

        final StructuralVariantLegPloidy negative = create(-1, Optional.empty(), Optional.of(4d));
        assertLeft(3d, 1d, negative);
        assertRight(0d, 0d, negative);
    }

    @Test
    public void testImpliedCopyNumberMissingRight()
    {
        final StructuralVariantLegPloidy positive = create(1, Optional.of(4d), Optional.empty());
        assertLeft(0d, 0d, positive);
        assertRight(3d, 1d, positive);

        final StructuralVariantLegPloidy negative = create(-1, Optional.of(3d), Optional.empty());
        assertLeft(0d, 0d, negative);
        assertRight(4d, 1d, negative);
    }

    private static void assertLeft(double expectedCopyNumber, double expectedWeight, @NotNull final StructuralVariantLegPloidy ploidy)
    {
        assertEquals(expectedCopyNumber, ploidy.impliedLeftCopyNumber(), EPSILON);
        assertEquals(expectedWeight, ploidy.impliedLeftCopyNumberWeight(), EPSILON);
    }

    private static void assertRight(double expectedCopyNumber, double expectedWeight, @NotNull final StructuralVariantLegPloidy ploidy)
    {
        assertEquals(expectedCopyNumber, ploidy.impliedRightCopyNumber(), EPSILON);
        assertEquals(expectedWeight, ploidy.impliedRightCopyNumberWeight(), EPSILON);
    }

    private static StructuralVariantLegPloidy create(
            int orientation, final Optional<Double> leftCopyNumber, final Optional<Double> rightCopyNumber)
    {
        return StructuralVariantPloidyFactoryTest.svLegPloidy(CHROMOSOME, orientation, leftCopyNumber, rightCopyNumber, PLOIDY);
    }
}
