package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantTestFactory;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.junit.Test;

public class VariantsTest {

    @Test
    public void canDedupVariants() {
        ReportableVariant variant1 = ReportableVariantTestFactory.builder()
                .canonicalEffect(VariantEffect.PHASED_INFRAME_DELETION.effect())
                .gene("EGFR")
                .canonicalHgvsCodingImpact("c.1")
                .canonicalHgvsProteinImpact("p.Glu746_Pro753delinsMetSer")
                .alleleCopyNumber(0.9)
                .build();

        ReportableVariant variant2 =
                ReportableVariantTestFactory.builder().from(variant1).canonicalHgvsCodingImpact("c.2").alleleCopyNumber(1.2).build();

        ReportableVariant variant3 =
                ReportableVariantTestFactory.builder().from(variant1).canonicalHgvsCodingImpact("c.3").alleleCopyNumber(0.9).build();

        ReportableVariant variant4 = ReportableVariantTestFactory.builder()
                .canonicalEffect(VariantEffect.FRAMESHIFT.effect())
                .gene("APC")
                .canonicalHgvsCodingImpact("p.Met1fs")
                .alleleCopyNumber(0.8)
                .build();

        List<ReportableVariant> variants = Lists.newArrayList(variant1, variant2, variant3, variant4);
        List<ReportableVariant> dedup = Variants.dedup(variants);

        assertEquals(2, dedup.size());
        assertTrue(dedup.contains(variant3));
        assertTrue(dedup.contains(variant4));
    }

    @Test
    public void canRenderRNADepthField() {
        ReportableVariant missingRNA = ReportableVariantTestFactory.builder().rnaAlleleReadCount(null).rnaTotalReadCount(null).build();
        assertEquals(ReportResources.NOT_AVAILABLE, Variants.rnaDepthField(missingRNA));

        ReportableVariant proper = ReportableVariantTestFactory.builder().rnaAlleleReadCount(10).rnaTotalReadCount(20).build();
        assertEquals("10/20 (50%)", Variants.rnaDepthField(proper));

        ReportableVariant noDepth = ReportableVariantTestFactory.builder().rnaAlleleReadCount(0).rnaTotalReadCount(0).build();
        assertEquals("0/0", Variants.rnaDepthField(noDepth));
    }
}