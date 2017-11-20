package com.hartwig.hmftools.common.purple.copynumber.sv;

import static com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegsFactory.reduce;
import static com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantPloidyTest.CHROMOSOME;

import static org.apache.commons.math3.util.Precision.EPSILON;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class StructuralVariantLegsFactoryTest {


    @Test
    public void testReduce() {
        final StructuralVariantLeg largePositive = createLeg(1, 0.9);
        final StructuralVariantLeg smallPositive = createLeg(1, 0.1);
        final StructuralVariantLeg largeNegative = createLeg(-1, 0.8);
        final StructuralVariantLeg smallNegative = createLeg(-1, 0.2);

        assertReduce(largePositive, Lists.newArrayList(largePositive, smallPositive));
        assertReduce(largePositive, Lists.newArrayList(smallPositive, largePositive));
        assertReduce(largeNegative, Lists.newArrayList(largeNegative, smallNegative));
        assertReduce(largeNegative, Lists.newArrayList(smallNegative, largeNegative));

        assertReduce(createLeg(1, 0.7), Lists.newArrayList(smallNegative, largePositive));
        assertReduce(createLeg(1, 0.1), Lists.newArrayList(largeNegative, largePositive));
        assertReduce(createLeg(-1, 0.7), Lists.newArrayList(largeNegative, smallPositive));
    }

    @Test
    public void testIgnoreInserts() {
        final StructuralVariant first = sv(100, 200, StructuralVariantType.INS, 0.25, 0.25);
        final StructuralVariant second = sv(300, 400, StructuralVariantType.INS, 0.25, 0.25);
        final List<StructuralVariant> variants = Lists.newArrayList(first, second);

        final List<StructuralVariantLegs> legs = StructuralVariantLegsFactory.create(variants);
        assertEquals(0, legs.size());
    }

    @Test
    public void testSameVAF() {
        final StructuralVariant first = sv(100, 200, StructuralVariantType.INV, 0.25, 0.25);
        final StructuralVariant second = sv(200, 300, StructuralVariantType.BND, 0.25, 0.25);
        final List<StructuralVariant> variants = Lists.newArrayList(first, second);

        final List<StructuralVariantLegs> legs = StructuralVariantLegsFactory.create(variants);
        assertEquals(2, legs.size());
        assertLeg(1, 0.25, legs.get(0).start().get());
        assertLeg(1, 0.25, legs.get(0).end().get());
        assertLeg(1, 0.25, legs.get(1).start().get());
        assertLeg(1, 0.25, legs.get(1).end().get());
    }

    @Test
    public void testMaxPositive() {
        final StructuralVariant first = sv(100, 200, StructuralVariantType.INV, 0.25, 0.25);
        final StructuralVariant second = sv(200, 300, StructuralVariantType.BND, 0.4, 0.25);
        final List<StructuralVariant> variants = Lists.newArrayList(first, second);

        final List<StructuralVariantLegs> legs = StructuralVariantLegsFactory.create(variants);
        assertEquals(2, legs.size());
        assertLeg(1, 0.25, legs.get(0).start().get());
        assertFalse(legs.get(0).end().isPresent());
        assertLeg(1, 0.4, legs.get(1).start().get());
        assertLeg(1, 0.25, legs.get(1).end().get());
    }

    @Test
    public void testOpposingSign() {
        final StructuralVariant first = sv(100, 200, StructuralVariantType.DEL, 0.25, 0.4);
        final StructuralVariant second = sv(200, 300, StructuralVariantType.BND, 0.25, 0.25);
        final List<StructuralVariant> variants = Lists.newArrayList(first, second);

        final List<StructuralVariantLegs> legs = StructuralVariantLegsFactory.create(variants);
        assertEquals(3, legs.size());
        assertLeg(1, 0.25, legs.get(0).start().get());
        assertFalse(legs.get(0).end().isPresent());
        assertFalse(legs.get(1).start().isPresent());
        assertLeg(1, 0.25, legs.get(1).end().get());
        assertLeg(-1, 0.15, legs.get(2).start().get());
        assertFalse(legs.get(2).end().isPresent());

    }
    private void assertLeg(int orientation, double vaf, @NotNull final StructuralVariantLeg victim) {
        assertEquals(vaf, victim.vaf(), EPSILON);
        assertEquals(orientation, victim.orientation());
    }

    private void assertReduce(@NotNull final StructuralVariantLeg expected, @NotNull final List<StructuralVariantLeg> legs) {
        final StructuralVariantLeg result = reduce(legs);
        assertEquals(expected.vaf(), result.vaf(), EPSILON);
        assertEquals(expected.orientation(), result.orientation());
    }

    @NotNull
    private StructuralVariantLeg createLeg(int orientation, double vaf) {
        return createLeg(1, orientation, vaf);
    }

    @NotNull
    static StructuralVariant sv(long start, long end, StructuralVariantType type, double startAF, double endAF) {
        return PurpleDatamodelTest.createStructuralVariant(CHROMOSOME, start, CHROMOSOME, end, type).startAF(startAF).endAF(endAF).build();
    }

    @NotNull
    private StructuralVariantLeg createLeg(long position, int orientation, double vaf) {
        return ImmutableStructuralVariantLeg.builder().chromosome(CHROMOSOME).position(position).orientation(orientation).vaf(vaf).build();
    }
}
