package com.hartwig.hmftools.common.protect;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.purple.loader.GainLoss;
import com.hartwig.hmftools.common.purple.loader.GainLossTestFactory;
import com.hartwig.hmftools.common.test.SomaticVariantTestFactory;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantTestFactory;
import com.hartwig.hmftools.common.variant.Variant;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class EventGeneratorTest {

    @Test
    public void canTestToVariantEvent() {
        assertEquals("p.Gly12Cys", EventGenerator.toVariantEvent("p.Gly12Cys", "c.123A>C", "missense_variant", CodingEffect.MISSENSE));
        assertEquals("c.123A>C splice", EventGenerator.toVariantEvent("p.?", "c.123A>C", "missense_variant", CodingEffect.SPLICE));
        assertEquals("c.123A>C", EventGenerator.toVariantEvent("", "c.123A>C", "missense_variant", CodingEffect.MISSENSE));
        assertEquals("splice", EventGenerator.toVariantEvent("", "", "splice", CodingEffect.SPLICE));
        assertEquals("missense_variant", EventGenerator.toVariantEvent("", "", "missense_variant", CodingEffect.MISSENSE));
    }

    @Test
    public void canGenerateEventForReportableVariant() {
        ReportableVariant base = ReportableVariantTestFactory.builder().isCanonical(true).canonicalHgvsCodingImpact("coding").build();
        assertEquals("coding", EventGenerator.variantEvent(base));

        ReportableVariant nonCanonical = ReportableVariantTestFactory.builder()
                .from(base)
                .isCanonical(false)
                .otherReportedEffects(createAltTranscriptInfo())
                .build();
        assertNotNull(EventGenerator.variantEvent(nonCanonical));
    }

    @NotNull
    private static String createAltTranscriptInfo() {
        return "trans|coding|impact|effect|MISSENSE";
    }

    @Test
    public void canGenerateEventForVariant() {
        Variant base = SomaticVariantTestFactory.builder().canonicalEffect("some effect").build();
        assertEquals("some effect", EventGenerator.variantEvent(base));

        Variant coding = ImmutableSomaticVariantImpl.builder().from(base).canonicalHgvsCodingImpact("coding impact").build();
        assertEquals("coding impact", EventGenerator.variantEvent(coding));

        Variant protein = ImmutableSomaticVariantImpl.builder().from(base).canonicalHgvsProteinImpact("protein impact").build();
        assertEquals("protein impact", EventGenerator.variantEvent(protein));
    }

    @Test
    public void canGenerateEventForCopyNumber() {
        GainLoss gainLoss = GainLossTestFactory.createTestGainLoss();
        assertEquals(gainLoss.interpretation().display(), EventGenerator.copyNumberEvent(gainLoss));
    }

    @Test
    public void canGenerateEventForFusion() {
        LinxFusion fusion = LinxTestFactory.fusionBuilder().geneStart("start").geneEnd("end").build();
        assertEquals("start - end fusion", EventGenerator.fusionEvent(fusion));
    }
}