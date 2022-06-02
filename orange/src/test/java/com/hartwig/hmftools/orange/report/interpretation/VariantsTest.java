package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.VariantTestFactory;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.junit.Test;

public class VariantsTest {

    @Test
    public void canRenderRNADepthField() {
        ReportableVariant missingRNA = VariantTestFactory.builder().rnaAlleleReadCount(null).rnaTotalReadCount(null).build();
        assertEquals(ReportResources.NOT_AVAILABLE, Variants.rnaDepthField(missingRNA));

        ReportableVariant proper = VariantTestFactory.builder().rnaAlleleReadCount(10).rnaTotalReadCount(20).build();
        assertEquals("10/20 (50%)", Variants.rnaDepthField(proper));

        ReportableVariant noDepth = VariantTestFactory.builder().rnaAlleleReadCount(0).rnaTotalReadCount(0).build();
        assertEquals("0/0", Variants.rnaDepthField(noDepth));
    }
}