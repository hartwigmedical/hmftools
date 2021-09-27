package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.common.drivercatalog.DriverImpact.isFrameshift;
import static com.hartwig.hmftools.common.drivercatalog.DriverImpact.isInframe;
import static com.hartwig.hmftools.common.drivercatalog.DriverImpact.isMissense;
import static com.hartwig.hmftools.common.drivercatalog.DriverImpact.isNonsense;
import static com.hartwig.hmftools.common.drivercatalog.DriverImpact.isSplice;

import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.test.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;
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
        assertTrue(isNonsense(nonsense));
        assertTrue(isMissense(missense));
        assertTrue(isSplice(spliceSNP));
        assertTrue(isSplice(spliceIndel));
        assertTrue(isFrameshift(frameshift));
        assertTrue(isInframe(inframe));
    }

    @Test
    public void testSorting() {
        List<SomaticVariant> list = Lists.newArrayList(missense, nonsense, spliceIndel, spliceSNP, inframe, frameshift);
        Collections.shuffle(list);

        list.sort(new TsgImpactComparator());
        assertTrue(isFrameshift(list.get(0)));
        assertTrue(isNonsense(list.get(1)));
        assertTrue(isSplice(list.get(2)));
        assertTrue(isSplice(list.get(3)));
        assertTrue(isMissense(list.get(4)));
        assertTrue(isInframe(list.get(5)));
    }

    @NotNull
    private static SomaticVariant create(@NotNull VariantType type, @NotNull CodingEffect codingEffect) {
        return SomaticVariantTestBuilderFactory.create().type(type).canonicalCodingEffect(codingEffect).build();
    }
}
