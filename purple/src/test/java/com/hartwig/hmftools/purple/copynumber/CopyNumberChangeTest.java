package com.hartwig.hmftools.purple.copynumber;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;
import com.hartwig.hmftools.purple.copynumber.sv.CopyNumberChange;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegPloidy;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CopyNumberChangeTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testTwoPositiveOrientationsHaveNoCopyNumberChange() {
        final StructuralVariantLegPloidy leg1 = svPloidy(1, 1, 2, 2);
        final StructuralVariantLegPloidy leg2 = svPloidy(1, 1.5, 2, 2);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2);
        assertPloidy(0, legList);
        assertCopyNumberChange(0, legList);

        final CopyNumberChange victim = new CopyNumberChange(legList);
        assertEquals(0, victim.copyNumberChange(leg1), EPSILON);
        assertEquals(0, victim.copyNumberChange(leg2), EPSILON);
    }

    @Test
    public void testTwoPositiveOrientationsHaveSmallCopyNumberChange() {
        final StructuralVariantLegPloidy leg1 = svPloidy(1, 1, 3, 2);
        final StructuralVariantLegPloidy leg2 = svPloidy(1, 1.5, 3, 2);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2);
        assertPloidy(1, legList);
        assertCopyNumberChange(-1, legList);

        final CopyNumberChange victim = new CopyNumberChange(legList);
        assertEquals(0.4, victim.copyNumberChange(leg1), EPSILON);
        assertEquals(0.6, victim.copyNumberChange(leg2), EPSILON);
    }

    @Test
    public void testTwoPositiveOrientationsHaveLargeCopyNumberChange() {
        final StructuralVariantLegPloidy leg1 = svPloidy(1, 1, 4, 1);
        final StructuralVariantLegPloidy leg2 = svPloidy(1, 1.5, 4, 1);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2);
        assertPloidy(3, legList);
        assertCopyNumberChange(-3, legList);

        final CopyNumberChange victim = new CopyNumberChange(legList);
        assertEquals(1.2, victim.copyNumberChange(leg1), EPSILON);
        assertEquals(1.8, victim.copyNumberChange(leg2), EPSILON);
    }

    @Test
    public void testTwoPositiveOrientationsIncreaseCopyNumber() {
        final StructuralVariantLegPloidy leg1 = svPloidy(1, 1, 3, 4.5);
        final StructuralVariantLegPloidy leg2 = svPloidy(1, 1.5, 3, 4.5);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2);
        assertPloidy(-1.5, legList);
        assertCopyNumberChange(1.5, legList);

        final CopyNumberChange victim = new CopyNumberChange(legList);
        assertEquals(-0.6, victim.copyNumberChange(leg1), EPSILON);
        assertEquals(-0.9, victim.copyNumberChange(leg2), EPSILON);
    }

    @Test
    public void testTwoNegativeOrientationsGivePositiveCopyNumberChange() {
        final StructuralVariantLegPloidy leg1 = svPloidy(-1, 1, 1, 4);
        final StructuralVariantLegPloidy leg2 = svPloidy(-1, 1.5, 1, 4);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2);
        assertPloidy(3, legList);
        assertCopyNumberChange(3, legList);

        final CopyNumberChange victim = new CopyNumberChange(legList);
        assertEquals(1.2, victim.copyNumberChange(leg1), EPSILON);
        assertEquals(1.8, victim.copyNumberChange(leg2), EPSILON);
    }

    @Test
    public void testTwoNegativeOrientationsGiveNegativeCopyNumberChange() {
        final StructuralVariantLegPloidy leg1 = svPloidy(-1, 1, 4, 1);
        final StructuralVariantLegPloidy leg2 = svPloidy(-1, 1.5, 4, 1);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2);
        assertPloidy(-3, legList);
        assertCopyNumberChange(-3, legList);

        final CopyNumberChange victim = new CopyNumberChange(legList);
        assertEquals(-1.2, victim.copyNumberChange(leg1), EPSILON);
        assertEquals(-1.8, victim.copyNumberChange(leg2), EPSILON);
    }

    @Test
    public void testTwoDifferentDirections() {
        final StructuralVariantLegPloidy leg1 = svPloidy(-1, 1, 2, 1.6);
        final StructuralVariantLegPloidy leg2 = svPloidy(1, 1.5, 2, 1.6);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2);
        assertPloidy(2.5, legList);
        assertCopyNumberChange(-0.4, legList);

        final CopyNumberChange victim = new CopyNumberChange(legList);
        assertEquals(1.05, victim.copyNumberChange(leg1), EPSILON);
        assertEquals(1.45, victim.copyNumberChange(leg2), EPSILON);
    }

    @Test
    public void testTwoDifferentDirectionsThroughZero() {
        final StructuralVariantLegPloidy leg1 = svPloidy(-1, 1, 1, 0.6);
        final StructuralVariantLegPloidy leg2 = svPloidy(1, 1.5, 1, 0.6);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2);
        assertPloidy(1.6, legList);
        assertCopyNumberChange(-0.4, legList);

        final CopyNumberChange victim = new CopyNumberChange(legList);
        assertEquals(0.6, victim.copyNumberChange(leg1), EPSILON);
        assertEquals(1, victim.copyNumberChange(leg2), EPSILON);
    }

    @Test
    public void testMultipleInDifferentDirection() {
        final StructuralVariantLegPloidy leg1 = svPloidy(-1, 1, 2, 1.4);
        final StructuralVariantLegPloidy leg2 = svPloidy(1, 0.7, 2, 1.4);
        final StructuralVariantLegPloidy leg3 = svPloidy(1, 0.8, 2, 1.4);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2, leg3);
        assertPloidy(2.5, legList);
        assertCopyNumberChange(-0.6, legList);

        final CopyNumberChange victim = new CopyNumberChange(legList);
        assertEquals(0.95, victim.copyNumberChange(leg1), EPSILON);
        assertEquals(0.743, victim.copyNumberChange(leg2), EPSILON);
        assertEquals(0.807, victim.copyNumberChange(leg3), EPSILON);
    }

    @Test
    public void testMultipleInDifferentDirectionPositive() {
        final StructuralVariantLegPloidy leg1 = svPloidy(-1, 0.6, 2, 2.5);
        final StructuralVariantLegPloidy leg2 = svPloidy(1, 1, 2, 2.5);
        final StructuralVariantLegPloidy leg3 = svPloidy(1, 1, 2, 2.5);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2, leg3);

        assertPloidy(2.6, legList);
        assertCopyNumberChange(0.5, legList);

        final CopyNumberChange victim = new CopyNumberChange(legList);
        assertEquals(1.55, victim.copyNumberChange(leg1), EPSILON);
        assertEquals(0.525, victim.copyNumberChange(leg2), EPSILON);
        assertEquals(0.525, victim.copyNumberChange(leg3), EPSILON);
    }

    @Test
    public void testNoCopyNumberChange() {
        final StructuralVariantLegPloidy leg1 = svPloidy(-1, 1, 2, 2);
        final StructuralVariantLegPloidy leg2 = svPloidy(1, 0.75, 2, 2);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2);
        final CopyNumberChange victim = new CopyNumberChange(legList);
        assertPloidy(1.75, legList);
        assertCopyNumberChange(0, legList);

        assertEquals(0.875, victim.copyNumberChange(leg1), EPSILON);
        assertEquals(0.875, victim.copyNumberChange(leg2), EPSILON);
    }

    @Test
    public void testNoCopyNumberChangeThroughZero() {
        final StructuralVariantLegPloidy leg1 = svPloidy(-1, 1, 0.8, 0.8);
        final StructuralVariantLegPloidy leg2 = svPloidy(1, 0.75, 0.8, 0.8);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2);
        final CopyNumberChange victim = new CopyNumberChange(legList);
        assertPloidy(1.55, legList);
        assertCopyNumberChange(0, legList);

        assertEquals(0.775, victim.copyNumberChange(leg1), EPSILON);
        assertEquals(0.775, victim.copyNumberChange(leg2), EPSILON);
    }

    private void assertCopyNumberChange(double expected, List<StructuralVariantLegPloidy> legList) {
        final CopyNumberChange victim = new CopyNumberChange(legList);

        double sumCopyNumberChange = legList.stream().mapToDouble(x -> -x.orientation() * victim.copyNumberChange(x)).sum();
        assertEquals(expected, sumCopyNumberChange, EPSILON);
    }

    private void assertPloidy(double expected, @NotNull final List<StructuralVariantLegPloidy> legList) {
        final CopyNumberChange victim = new CopyNumberChange(legList);

        double sumCopyNumberChange = legList.stream().mapToDouble(victim::copyNumberChange).sum();
        assertEquals(expected, sumCopyNumberChange, EPSILON);
    }

    @NotNull
    private static StructuralVariantLegPloidy svPloidy(int orientation, double ploidy, double leftCopyNumber, double rightCopyNumber) {
        return StructuralVariantPloidyFactoryTest.svLegPloidy(orientation, Optional.of(leftCopyNumber), Optional.of(rightCopyNumber), ploidy).build();
    }
}
