package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.common.purple.PurpleTestUtils.createStartLeg;
import static com.hartwig.hmftools.purple.copynumber.StructuralVariantPloidyTest.CHROMOSOME;
import static com.hartwig.hmftools.purple.copynumber.sv.SvLegsFactory.reduce;

import static org.apache.commons.math3.util.Precision.EPSILON;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegs;
import com.hartwig.hmftools.purple.copynumber.sv.SvLegsFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class StructuralVariantLegsFactoryTest
{
    @Test
    public void testReduce()
    {
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
    public void testIgnoreInserts()
    {
        final StructuralVariant first = sv(100, 200, StructuralVariantType.INS, 0.25, 0.25);
        final StructuralVariant second = sv(300, 400, StructuralVariantType.INS, 0.25, 0.25);
        final List<StructuralVariant> variants = Lists.newArrayList(first, second);

        final List<StructuralVariantLegs> legs = SvLegsFactory.create(variants);
        assertEquals(0, legs.size());
    }

    @Test
    public void testDuplicateLeg()
    {
        final StructuralVariant first = sv(100, 100, StructuralVariantType.INV, 0.25, 0.25);

        final List<StructuralVariantLegs> legs = SvLegsFactory.create(Lists.newArrayList(first));
        assertEquals(1, legs.size());
        assertNotNull(legs.get(0).start());
        assertNull(legs.get(0).end());
    }

    @Test
    public void testSameVAF()
    {
        final StructuralVariant first = sv(100, 200, StructuralVariantType.INV, 0.25, 0.25);
        final StructuralVariant second = sv(200, 300, StructuralVariantType.BND, 0.25, 0.25);
        final List<StructuralVariant> variants = Lists.newArrayList(first, second);

        final List<StructuralVariantLegs> legs = SvLegsFactory.create(variants);
        assertEquals(2, legs.size());

        assertNotNull(legs.get(0).start());
        assertNotNull(legs.get(0).end());
        assertNull(legs.get(1).start());
        assertNotNull(legs.get(1).end());

        assertLeg(1, 0.25, legs.get(0).start());
        assertLeg(1, 0.25, legs.get(0).end());
        assertLeg(1, 0.25, legs.get(1).end());
    }

    @Test
    public void testMaxPositive()
    {
        final StructuralVariant first = sv(100, 200, StructuralVariantType.INV, 0.25, 0.25);
        final StructuralVariant second = sv(200, 300, StructuralVariantType.BND, 0.4, 0.25);
        final List<StructuralVariant> variants = Lists.newArrayList(first, second);

        final List<StructuralVariantLegs> legs = SvLegsFactory.create(variants);
        assertEquals(2, legs.size());
        assertNotNull(legs.get(0).start());
        assertLeg(1, 0.25, legs.get(0).start());
        assertNull(legs.get(0).end());
        assertNotNull(legs.get(1).start());
        assertLeg(1, 0.4, legs.get(1).start());
        assertNotNull(legs.get(1).end());
        assertLeg(1, 0.25, legs.get(1).end());
    }

    @Test
    public void testOpposingSign()
    {
        final StructuralVariant first = sv(100, 200, StructuralVariantType.DEL, 0.25, 0.4);
        final StructuralVariant second = sv(199, 300, StructuralVariantType.BND, 0.25, 0.25);
        final List<StructuralVariant> variants = Lists.newArrayList(first, second);

        final List<StructuralVariantLegs> legs = SvLegsFactory.create(variants);
        assertEquals(3, legs.size());
        assertNotNull(legs.get(0).start());
        assertLeg(1, 0.25, legs.get(0).start());
        assertNull(legs.get(0).end());
        assertNull(legs.get(1).start());
        assertNotNull(legs.get(1).end());
        assertLeg(1, 0.25, legs.get(1).end());
        assertNotNull(legs.get(2).start());
        assertLeg(-1, 0.15, legs.get(2).start());
        assertNull(legs.get(2).end());
    }

    @Test
    public void testSingleBreakend()
    {
        final StructuralVariant first = breakend(100, 0.25);

        final List<StructuralVariantLegs> legs = SvLegsFactory.create(Lists.newArrayList(first));
        assertEquals(1, legs.size());
        assertNotNull(legs.get(0).start());
        assertNull(legs.get(0).end());
    }

    private void assertLeg(int orientation, double vaf, @NotNull final StructuralVariantLeg victim)
    {
        assertEquals(vaf, victim.alleleFrequency(), EPSILON);
        assertEquals(orientation, victim.orientation());
    }

    private void assertReduce(@NotNull final StructuralVariantLeg expected, @NotNull final List<StructuralVariantLeg> legs)
    {
        final StructuralVariantLeg result = reduce(expected, legs);
        assertEquals(expected.alleleFrequency(), result.alleleFrequency(), EPSILON);
        assertEquals(expected.orientation(), result.orientation());
    }

    @NotNull
    private static StructuralVariantLeg createLeg(int orientation, double vaf)
    {
        return createLeg(1, orientation, vaf);
    }

    @NotNull
    private static StructuralVariant sv(int start, int end, StructuralVariantType type, double startAF, double endAF)
    {
        return PurpleTestUtils.createStructuralVariant(CHROMOSOME, start, CHROMOSOME, end, type, startAF, endAF).build();
    }

    @NotNull
    private static StructuralVariant breakend(int start, double startAF)
    {
        return PurpleTestUtils.createStructuralVariantSingleBreakend(CHROMOSOME, start, startAF).build();
    }

    @NotNull
    static StructuralVariantLeg createLeg(int position, int orientation, double vaf)
    {
        return createStartLeg(CHROMOSOME, position, StructuralVariantType.DEL)
                .orientation((byte) orientation)
                .alleleFrequency(vaf)
                .build();
    }

    @NotNull
    static StructuralVariantLeg createLeg(int position, int orientation, double vaf, int tumorVariantFragmentCount)
    {
        return createStartLeg(CHROMOSOME, position, StructuralVariantType.DEL)
                .orientation((byte) orientation)
                .alleleFrequency(vaf)
                .tumorVariantFragmentCount(tumorVariantFragmentCount)
                .build();
    }
}
