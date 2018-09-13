package com.hartwig.hmftools.common.drivercatalog;

import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

public class TsgImpactComparatorTest {

    private EnrichedSomaticVariant missense;
    private EnrichedSomaticVariant nonsense;
    private EnrichedSomaticVariant spliceIndel;
    private EnrichedSomaticVariant spliceSNP;
    private EnrichedSomaticVariant frameshift;
    private EnrichedSomaticVariant inframe;

    @Before
    public void setup() {
        nonsense = create(VariantType.MNP, CodingEffect.NONSENSE_OR_FRAMESHIFT);
        missense = create(VariantType.SNP, CodingEffect.MISSENSE);
        spliceIndel = create(VariantType.INDEL, CodingEffect.SPLICE);
        spliceSNP = create(VariantType.SNP, CodingEffect.SPLICE);
        frameshift = create(VariantType.INDEL, CodingEffect.NONSENSE_OR_FRAMESHIFT);
        inframe = create(VariantType.INDEL, CodingEffect.MISSENSE);
    }

    @Test
    public void testSetup() {
        assertTrue(TsgImpactComparator.isNonsense(nonsense));
        assertTrue(TsgImpactComparator.isMissense(missense));
        assertTrue(TsgImpactComparator.isSplice(spliceSNP));
        assertTrue(TsgImpactComparator.isSplice(spliceIndel));
        assertTrue(TsgImpactComparator.isFrameshift(frameshift));
        assertTrue(TsgImpactComparator.isInframe(inframe));
    }

    @Test
    public void testSorting() {
        final List<EnrichedSomaticVariant> list = Lists.newArrayList(missense, nonsense, spliceIndel, spliceSNP, inframe, frameshift);
        Collections.shuffle(list);

        list.sort(new TsgImpactComparator());
        assertTrue(TsgImpactComparator.isFrameshift(list.get(0)));
        assertTrue(TsgImpactComparator.isNonsense(list.get(1)));
        assertTrue(TsgImpactComparator.isSplice(list.get(2)));
        assertTrue(TsgImpactComparator.isSplice(list.get(3)));
        assertTrue(TsgImpactComparator.isMissense(list.get(4)));
        assertTrue(TsgImpactComparator.isInframe(list.get(5)));
    }

    @NotNull
    private EnrichedSomaticVariant create(@NotNull final VariantType type, @NotNull final CodingEffect codingEffect) {
        return SomaticVariantTestBuilderFactory.createEnriched().type(type).canonicalCodingEffect(codingEffect).build();
    }

}
