package com.hartwig.hmftools.common.protect;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.protect.variant.OtherEffectsTestFactory;
import com.hartwig.hmftools.common.purple.interpretation.GainLoss;
import com.hartwig.hmftools.common.purple.interpretation.GainLossTestFactory;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.test.SomaticVariantTestFactory;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantTestFactory;
import com.hartwig.hmftools.common.variant.Variant;

import org.junit.Test;

public class ProtectEventGeneratorTest {

    @Test
    public void canTestToVariantEvent() {
        assertEquals("p.Gly12Cys", ProtectEventGenerator.toVariantEvent("p.Gly12Cys", "c.123A>C", "missense_variant"));
        assertEquals("c.123A>C", ProtectEventGenerator.toVariantEvent("p.?", "c.123A>C", "missense_variant"));
        assertEquals("c.123A>C", ProtectEventGenerator.toVariantEvent("", "c.123A>C", "missense_variant"));
        assertEquals("upstream", ProtectEventGenerator.toVariantEvent("", "", "upstream_gene_variant"));
        assertEquals("missense_variant", ProtectEventGenerator.toVariantEvent("", "", "missense_variant"));
    }

    @Test
    public void canGenerateEventForReportableVariant() {
        ReportableVariant base = ReportableVariantTestFactory.builder().isCanonical(true).canonicalHgvsCodingImpact("coding").build();
        assertEquals("coding", ProtectEventGenerator.variantEvent(base));

        ReportableVariant nonCanonical = ReportableVariantTestFactory.builder()
                .from(base)
                .isCanonical(false)
                .otherReportedEffects(OtherEffectsTestFactory.create())
                .build();
        assertNotNull(ProtectEventGenerator.variantEvent(nonCanonical));
    }

    @Test
    public void canGenerateEventForVariant() {
        Variant base = SomaticVariantTestFactory.builder().canonicalEffect("some effect").build();
        assertEquals("some effect", ProtectEventGenerator.variantEvent(base));

        String upstreamEffect = ProtectEventGenerator.UPSTREAM_GENE_VARIANT;
        Variant upstream = ImmutableSomaticVariantImpl.builder().from(base).canonicalEffect(upstreamEffect).build();
        assertEquals("upstream", ProtectEventGenerator.variantEvent(upstream));

        Variant coding = ImmutableSomaticVariantImpl.builder().from(upstream).canonicalHgvsCodingImpact("coding impact").build();
        assertEquals("coding impact", ProtectEventGenerator.variantEvent(coding));

        Variant protein = ImmutableSomaticVariantImpl.builder().from(coding).canonicalHgvsProteinImpact("protein impact").build();
        assertEquals("protein impact", ProtectEventGenerator.variantEvent(protein));
    }

    @Test
    public void canGenerateEventForCopyNumber() {
        GainLoss gainLoss = GainLossTestFactory.createTestGainLoss();
        assertEquals(gainLoss.interpretation().display(), ProtectEventGenerator.copyNumberEvent(gainLoss));
    }

    @Test
    public void canGenerateEventForFusion() {
        LinxFusion fusion = LinxTestFactory.builder().geneStart("start").geneEnd("end").build();
        assertEquals("start - end fusion", ProtectEventGenerator.fusionEvent(fusion));
    }
}