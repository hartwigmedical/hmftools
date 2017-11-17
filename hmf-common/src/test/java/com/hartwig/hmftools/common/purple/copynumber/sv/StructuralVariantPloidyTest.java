package com.hartwig.hmftools.common.purple.copynumber.sv;

import static org.junit.Assert.assertEquals;

import java.util.Optional;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

@SuppressWarnings("OptionalUsedAsFieldOrParameterType")
public class StructuralVariantPloidyTest {

    private static final String CHROMOSOME = "1";
    private static final double PLOIDY = 1;
    private static final double EPSILON = 1e-10;

    @Test
    public void testCompleteImpliedCopyNumbers() {
        final StructuralVariantPloidy positive = create(1, Optional.of(4d), Optional.of(3d));
        assertLeft(4d, 1d, positive);
        assertRight(3d, 1d, positive);

        final StructuralVariantPloidy negative = create(-1, Optional.of(3d), Optional.of(4d));
        assertLeft(3d, 1d, negative);
        assertRight(4d, 1d, negative);
    }

    @Test
    public void testImpliedCopyNumberMissingLeft() {
        final StructuralVariantPloidy positive = create(1, Optional.empty(), Optional.of(3d));
        assertLeft(4d, 1d, positive);
        assertRight(0d, 0d, positive);

        final StructuralVariantPloidy negative = create(-1, Optional.empty(), Optional.of(4d));
        assertLeft(3d, 1d, negative);
        assertRight(0d, 0d, negative);
    }

    @Test
    public void testImpliedCopyNumberMissingRight() {
        final StructuralVariantPloidy positive = create(1, Optional.of(4d), Optional.empty());
        assertLeft(0d, 0d, positive);
        assertRight(3d, 1d, positive);

        final StructuralVariantPloidy negative = create(-1, Optional.of(3d), Optional.empty());
        assertLeft(0d, 0d, negative);
        assertRight(4d, 1d, negative);
    }


    private void assertLeft(double expectedCopyNumber, double expectedWeight, @NotNull final StructuralVariantPloidy ploidy) {
        assertEquals(expectedCopyNumber, ploidy.impliedLeftCopyNumber(), EPSILON);
        assertEquals(expectedWeight, ploidy.impliedLeftCopyNumberWeight(), EPSILON);
    }

    private void assertRight(double expectedCopyNumber, double expectedWeight, @NotNull final StructuralVariantPloidy ploidy) {
        assertEquals(expectedCopyNumber, ploidy.impliedRightCopyNumber(), EPSILON);
        assertEquals(expectedWeight, ploidy.impliedRightCopyNumberWeight(), EPSILON);
    }

    @NotNull
    static StructuralVariantPloidy create(int orientation, @NotNull final Optional<Double> leftCopyNumber,
            @NotNull final Optional<Double> rightCopyNumber) {
        return ImmutableStructuralVariantPloidy.builder()
                .chromosome(CHROMOSOME)
                .position(1)
                .orientation(orientation)
                .vaf(0.5)
                .weight(1)
                .averageImpliedPloidy(PLOIDY)
                .unweightedImpliedPloidy(PLOIDY)
                .leftCopyNumber(leftCopyNumber)
                .rightCopyNumber(rightCopyNumber)
                .build();
    }
}
