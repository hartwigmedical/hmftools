package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.common.drivercatalog.DriverImpact.isFrameshift;
import static com.hartwig.hmftools.common.drivercatalog.DriverImpact.isInframe;
import static com.hartwig.hmftools.common.drivercatalog.DriverImpact.isMissense;
import static com.hartwig.hmftools.common.drivercatalog.DriverImpact.isNonsense;
import static com.hartwig.hmftools.common.drivercatalog.DriverImpact.isSplice;
import static com.hartwig.hmftools.common.variant.Hotspot.NON_HOTSPOT;
import static com.hartwig.hmftools.purple.MiscTestUtils.createVariant;

import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.test.SomaticVariantTestFactory;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

import org.junit.Before;
import org.junit.Test;

public class TsgImpactComparatorTest {

    private SomaticVariant missense;
    private SomaticVariant nonsense;
    private SomaticVariant spliceIndel;
    private SomaticVariant spliceSNP;
    private SomaticVariant frameshift;
    private SomaticVariant inframe;

    @Before
    public void setup()
    {
        nonsense = createVariant(VariantType.MNP, CodingEffect.NONSENSE_OR_FRAMESHIFT, 0, NON_HOTSPOT, 0.5);
        missense = createVariant(VariantType.SNP, CodingEffect.MISSENSE, 0, NON_HOTSPOT, 0.5);
        spliceIndel = createVariant(VariantType.INDEL, CodingEffect.SPLICE, 0, NON_HOTSPOT, 0.5);
        spliceSNP = createVariant(VariantType.SNP, CodingEffect.SPLICE, 0, NON_HOTSPOT, 0.5);
        frameshift = createVariant(VariantType.INDEL, CodingEffect.NONSENSE_OR_FRAMESHIFT, 0, NON_HOTSPOT, 0.5);
        inframe = createVariant(VariantType.INDEL, CodingEffect.MISSENSE, 0, NON_HOTSPOT, 0.5);
    }

    @Test
    public void testSetup()
    {
        assertTrue(isNonsense(nonsense.type(), nonsense.variantImpact().CanonicalCodingEffect));
        assertTrue(isMissense(missense.type(), missense.variantImpact().CanonicalCodingEffect));
        assertTrue(isSplice(spliceSNP.variantImpact().CanonicalCodingEffect));
        assertTrue(isSplice(spliceIndel.variantImpact().CanonicalCodingEffect));
        assertTrue(isFrameshift(frameshift.type(), frameshift.variantImpact().CanonicalCodingEffect));
        assertTrue(isInframe(inframe.type(), inframe.variantImpact().CanonicalCodingEffect));
    }

    @Test
    public void testSorting()
    {
        List<SomaticVariant> list = Lists.newArrayList(missense, nonsense, spliceIndel, spliceSNP, inframe, frameshift);
        Collections.shuffle(list);

        list.sort(new TsgImpactComparator());
        assertTrue(isFrameshift(list.get(0).type(), list.get(0).variantImpact().CanonicalCodingEffect));
        assertTrue(isNonsense(list.get(1).type(), list.get(1).variantImpact().CanonicalCodingEffect));
        assertTrue(isSplice(list.get(2).variantImpact().CanonicalCodingEffect));
        assertTrue(isSplice(list.get(3).variantImpact().CanonicalCodingEffect));
        assertTrue(isMissense(list.get(4).type(), list.get(5).variantImpact().CanonicalCodingEffect));
        assertTrue(isInframe(list.get(5).type(), list.get(5).variantImpact().CanonicalCodingEffect));
    }

    private static com.hartwig.hmftools.common.variant.SomaticVariant create(final VariantType type, final CodingEffect codingEffect)
    {
        return SomaticVariantTestFactory.builder().type(type).canonicalCodingEffect(codingEffect).build();
    }
}
