package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.orange.algo.purple.PurpleVariant;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleVariantFactory;

import org.junit.Test;

public class VariantDedupTest {

    @Test
    public void canDedupVariants() {
        PurpleVariant variant1 = TestPurpleVariantFactory.builder()
                .gene("EGFR")
                .canonicalImpact(TestPurpleVariantFactory.impactBuilder()
                        .hgvsCodingImpact("c.1")
                        .hgvsProteinImpact("p.Glu746_Pro753delinsMetSer")
                        .addEffects(VariantEffect.PHASED_INFRAME_DELETION)
                        .build())
                .variantCopyNumber(0.9)
                .build();

        PurpleVariant variant2 = TestPurpleVariantFactory.builder()
                .from(variant1)
                .canonicalImpact(TestPurpleVariantFactory.impactBuilder().from(variant1.canonicalImpact()).hgvsCodingImpact("c.2").build())
                .variantCopyNumber(1.2)
                .build();

        PurpleVariant variant3 = TestPurpleVariantFactory.builder()
                .from(variant1)
                .canonicalImpact(TestPurpleVariantFactory.impactBuilder().from(variant1.canonicalImpact()).hgvsCodingImpact("c.3").build())
                .variantCopyNumber(0.9)
                .build();

        PurpleVariant variant4 = TestPurpleVariantFactory.builder()
                .gene("APC")
                .canonicalImpact(TestPurpleVariantFactory.impactBuilder()
                        .hgvsProteinImpact("p.Met1fs")
                        .addEffects(VariantEffect.FRAMESHIFT)
                        .build())
                .variantCopyNumber(0.8)
                .build();

        List<PurpleVariant> variants = Lists.newArrayList(variant1, variant2, variant3, variant4);
        List<PurpleVariant> dedup = VariantDedup.apply(variants);

        assertEquals(2, dedup.size());
        assertTrue(dedup.contains(variant3));
        assertTrue(dedup.contains(variant4));
    }
}