package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.orange.algo.purple.PurpleVariant;
import com.hartwig.hmftools.orange.algo.purple.PurpleVariantTestFactory;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.junit.Test;

public class VariantsTest {

    @Test
    public void canDedupVariants() {
        PurpleVariant variant1 = PurpleVariantTestFactory.builder()
                .gene("EGFR")
                .canonicalImpact(PurpleVariantTestFactory.impactBuilder()
                        .hgvsCodingImpact("c.1")
                        .hgvsProteinImpact("p.Glu746_Pro753delinsMetSer")
                        .addEffects(VariantEffect.PHASED_INFRAME_DELETION)
                        .build())
                .variantCopyNumber(0.9)
                .build();

        PurpleVariant variant2 = PurpleVariantTestFactory.builder()
                .from(variant1)
                .canonicalImpact(PurpleVariantTestFactory.impactBuilder().from(variant1.canonicalImpact()).hgvsCodingImpact("c.2").build())
                .variantCopyNumber(1.2)
                .build();

        PurpleVariant variant3 = PurpleVariantTestFactory.builder()
                .from(variant1)
                .canonicalImpact(PurpleVariantTestFactory.impactBuilder().from(variant1.canonicalImpact()).hgvsCodingImpact("c.3").build())
                .variantCopyNumber(0.9)
                .build();

        PurpleVariant variant4 = PurpleVariantTestFactory.builder()
                .gene("APC")
                .canonicalImpact(PurpleVariantTestFactory.impactBuilder()
                        .hgvsProteinImpact("p.Met1fs")
                        .addEffects(VariantEffect.FRAMESHIFT)
                        .build())
                .variantCopyNumber(0.8)
                .build();

        List<PurpleVariant> variants = Lists.newArrayList(variant1, variant2, variant3, variant4);
        List<PurpleVariant> dedup = Variants.dedup(variants);

        assertEquals(2, dedup.size());
        assertTrue(dedup.contains(variant3));
        assertTrue(dedup.contains(variant4));
    }

    @Test
    public void canRenderRNADepthField() {
        PurpleVariant missingRNA = PurpleVariantTestFactory.builder().rnaDepth(null).build();
        assertEquals(ReportResources.NOT_AVAILABLE, Variants.rnaDepthField(missingRNA));

        PurpleVariant proper = PurpleVariantTestFactory.builder()
                .rnaDepth(PurpleVariantTestFactory.depthBuilder().alleleReadCount(10).totalReadCount(20).build())
                .build();
        assertEquals("10/20 (50%)", Variants.rnaDepthField(proper));

        PurpleVariant noDepth = PurpleVariantTestFactory.builder()
                .rnaDepth(PurpleVariantTestFactory.depthBuilder().alleleReadCount(0).totalReadCount(0).build())
                .build();
        assertEquals("0/0", Variants.rnaDepthField(noDepth));
    }

    @Test
    public void canGenerateVariantEvents() {
        assertEquals("p.Gly12Cys",
                Variants.toVariantEvent("p.Gly12Cys", "c.123A>C", Sets.newHashSet(VariantEffect.MISSENSE), CodingEffect.MISSENSE));
        assertEquals("c.123A>C splice",
                Variants.toVariantEvent("p.?", "c.123A>C", Sets.newHashSet(VariantEffect.MISSENSE), CodingEffect.SPLICE));
        assertEquals("c.123A>C", Variants.toVariantEvent("", "c.123A>C", Sets.newHashSet(VariantEffect.MISSENSE), CodingEffect.MISSENSE));
        assertEquals("missense_variant", Variants.toVariantEvent("", "", Sets.newHashSet(VariantEffect.MISSENSE), CodingEffect.MISSENSE));
    }
}