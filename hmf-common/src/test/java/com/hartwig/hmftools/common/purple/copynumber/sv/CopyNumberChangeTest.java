package com.hartwig.hmftools.common.purple.copynumber.sv;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CopyNumberChangeTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testTwoPositiveOrientationsGivePositiveCopyNumberChange() {
        double copyNumberChange = 1;

        final StructuralVariantLegPloidy leg1 = svPloidy(1, 1, copyNumberChange);
        final StructuralVariantLegPloidy leg2 = svPloidy(1, 1.5, copyNumberChange);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2);
        assertCopyNumberChange(copyNumberChange, legList);

        final CopyNumberChange victim = new CopyNumberChange(legList);
        assertEquals(-0.4, victim.copyNumberChange(leg1), EPSILON);
        assertEquals(-0.6, victim.copyNumberChange(leg2), EPSILON);
    }

    @Test
    public void testTwoPositiveOrientationsGiveNegativeCopyNumberChange() {
        double copyNumberChange = -2;

        final StructuralVariantLegPloidy leg1 = svPloidy(1, 1, copyNumberChange);
        final StructuralVariantLegPloidy leg2 = svPloidy(1, 1.5, copyNumberChange);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2);
        assertCopyNumberChange(copyNumberChange, legList);

        final CopyNumberChange victim = new CopyNumberChange(legList);
        assertEquals(0.8, victim.copyNumberChange(leg1), EPSILON);
        assertEquals(1.2, victim.copyNumberChange(leg2), EPSILON);
    }

    @Test
    public void testTwoNegativeOrientationsGivePositiveCopyNumberChange() {
        double copyNumberChange = 3;

        final StructuralVariantLegPloidy leg1 = svPloidy(-1, 1, copyNumberChange);
        final StructuralVariantLegPloidy leg2 = svPloidy(-1, 1.5, copyNumberChange);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2);
        assertCopyNumberChange(copyNumberChange, legList);

        final CopyNumberChange victim = new CopyNumberChange(legList);
        assertEquals(1.2, victim.copyNumberChange(leg1), EPSILON);
        assertEquals(1.8, victim.copyNumberChange(leg2), EPSILON);
    }

    @Test
    public void testTwoNegativeOrientationsGiveNegativeCopyNumberChange() {
        double copyNumberChange = -3;

        final StructuralVariantLegPloidy leg1 = svPloidy(-1, 1, copyNumberChange);
        final StructuralVariantLegPloidy leg2 = svPloidy(-1, 1.5, copyNumberChange);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2);
        assertCopyNumberChange(copyNumberChange, legList);

        final CopyNumberChange victim = new CopyNumberChange(legList);
        assertEquals(-1.2, victim.copyNumberChange(leg1), EPSILON);
        assertEquals(-1.8, victim.copyNumberChange(leg2), EPSILON);
    }

    @Test
    public void testTwoDifferentDirections() {

        double copyNumberChange = -0.4;

        final StructuralVariantLegPloidy leg1 = svPloidy(-1, 1, copyNumberChange);
        final StructuralVariantLegPloidy leg2 = svPloidy(1, 1.5, copyNumberChange);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2);
        //        assertCopyNumberChange(copyNumberChange, legList);

        final CopyNumberChange victim = new CopyNumberChange(legList);
        assertEquals(1.05, victim.copyNumberChange(leg1), EPSILON);
        assertEquals(1.45, victim.copyNumberChange(leg2), EPSILON);

        System.out.println(victim.copyNumberChange(leg1));
        System.out.println(victim.copyNumberChange(leg2));

    }

    @Test
    public void testMultipleInDifferentDirection() {

        double copyNumberChange = -0.6;

        final StructuralVariantLegPloidy leg1 = svPloidy(-1, 1, copyNumberChange);
        final StructuralVariantLegPloidy leg2 = svPloidy(1, 0.7, copyNumberChange);
        final StructuralVariantLegPloidy leg3 = svPloidy(1, 0.8, copyNumberChange);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2, leg3);
        assertAbsCopyNumberChange(legList);
        assertCopyNumberChange(copyNumberChange, legList);

        final CopyNumberChange victim = new CopyNumberChange(legList);
        assertEquals(0.95, victim.copyNumberChange(leg1), EPSILON);
        assertEquals(0.743, victim.copyNumberChange(leg2), EPSILON);
        assertEquals(0.807, victim.copyNumberChange(leg3), EPSILON);
    }

    @Test
    public void testMultipleInDifferentDirectionPositive() {

        double copyNumberChange = -0.5;

        final StructuralVariantLegPloidy leg1 = svPloidy(-1, 0.6, copyNumberChange);
        final StructuralVariantLegPloidy leg2 = svPloidy(1, 1, copyNumberChange);
        final StructuralVariantLegPloidy leg3 = svPloidy(1, 1, copyNumberChange);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2, leg3);
        final CopyNumberChange victim = new CopyNumberChange(legList);

        assertCopyNumberChange(copyNumberChange, legList);
        assertAbsCopyNumberChange(legList);

        double c1 = victim.copyNumberChange(leg1);
        double c2 = victim.copyNumberChange(leg2);
        double c3 = victim.copyNumberChange(leg3);

        System.out.println(c1 - c2 - c3);
        System.out.println(c1 + c2 + c3);

        //
        //        assertEquals(0.6, victim.copyNumberChange(leg1), EPSILON);
        //        assertEquals(0.55, victim.copyNumberChange(leg2), EPSILON);
        //        assertEquals(0.55, victim.copyNumberChange(leg3), EPSILON);
    }

    @Test
    public void testNoCopyNumberChange() {

        double copyNumberChange = 0;

        final StructuralVariantLegPloidy leg1 = svPloidy(-1, 1, copyNumberChange);
        final StructuralVariantLegPloidy leg2 = svPloidy(1, 0.75, copyNumberChange);
        final List<StructuralVariantLegPloidy> legList = Lists.newArrayList(leg1, leg2);
        final CopyNumberChange victim = new CopyNumberChange(legList);
        assertAbsCopyNumberChange(legList);
        assertCopyNumberChange(copyNumberChange, legList);

        assertEquals(0.875, victim.copyNumberChange(leg1), EPSILON);
        assertEquals(0.875, victim.copyNumberChange(leg2), EPSILON);
    }

    private void assertCopyNumberChange(double expected, List<StructuralVariantLegPloidy> legList) {
        final CopyNumberChange victim = new CopyNumberChange(legList);

        double sumCopyNumberChange = legList.stream().mapToDouble(x -> -x.orientation() * victim.copyNumberChange(x)).sum();
        assertEquals(expected, sumCopyNumberChange, EPSILON);
    }

    private void assertAbsCopyNumberChange(@NotNull final List<StructuralVariantLegPloidy> legList) {
        final CopyNumberChange victim = new CopyNumberChange(legList);

        double sumPloidy = legList.stream().mapToDouble(StructuralVariantLegPloidy::averageImpliedPloidy).sum();
        double sumCopyNumberChange = legList.stream().mapToDouble(victim::copyNumberChange).sum();
        assertEquals(sumPloidy, sumCopyNumberChange, EPSILON);
    }

    @NotNull
    private static StructuralVariantLegPloidy svPloidy(int orientation, double ploidy, double copyNumberChange) {
        return PurpleDatamodelTest.svLegPloidy(orientation, Optional.of(0d), Optional.of(copyNumberChange), ploidy).build();
    }

}
