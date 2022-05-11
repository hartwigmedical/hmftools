package com.hartwig.hmftools.common.protect;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.purple.PurpleTestFactory;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.test.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.Variant;

import org.junit.Test;

public class ProtectEventGeneratorTest {

    @Test
    public void canGenerateEventForVariant() {
        Variant base = SomaticVariantTestBuilderFactory.create().canonicalEffect("some effect").build();
        assertEquals("some effect", ProtectEventGenerator.variantEvent(base));

        Variant upstream =
                ImmutableSomaticVariantImpl.builder().from(base).canonicalEffect(ProtectEventGenerator.UPSTREAM_GENE_VARIANT).build();
        assertEquals("upstream", ProtectEventGenerator.variantEvent(upstream));

        Variant coding = ImmutableSomaticVariantImpl.builder().from(upstream).canonicalHgvsCodingImpact("coding impact").build();
        assertEquals("coding impact", ProtectEventGenerator.variantEvent(coding));

        Variant protein = ImmutableSomaticVariantImpl.builder().from(coding).canonicalHgvsProteinImpact("protein impact").build();
        assertEquals("protein impact", ProtectEventGenerator.variantEvent(protein));
    }

    @Test
    public void canGenerateEventForCopyNumber() {
        ReportableGainLoss gainLoss = PurpleTestFactory.createTestReportableGainLoss();
        assertEquals(gainLoss.interpretation().display(), ProtectEventGenerator.copyNumberEvent(gainLoss));
    }

    @Test
    public void canGenerateEventForFusion() {
        LinxFusion fusion = LinxTestFactory.testBuilder().geneStart("start").geneEnd("end").build();
        assertEquals("start - end fusion", ProtectEventGenerator.fusionEvent(fusion));
    }
}