package com.hartwig.hmftools.common.drivercatalog;

import static com.hartwig.hmftools.common.drivercatalog.ReportablePredicate.MAX_ONCO_REPEAT_COUNT;

import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelFactoryTest;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.common.variant.VariantType;

import org.junit.Test;

public class ReportablePredicateTest {

    private final DriverGenePanel genePanel = DriverGenePanelFactoryTest.testGenePanel();

    @Test
    public void testIgnoreIndelsWithLargeRepeatCount() {
        final SomaticVariant variant = SomaticVariantTestBuilderFactory.create()
                .gene("AR")
                .repeatCount(MAX_ONCO_REPEAT_COUNT)
                .type(VariantType.INDEL)
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .build();

        final SomaticVariant variantLargeRepeatCount =
                ImmutableSomaticVariantImpl.builder().from(variant).repeatCount(MAX_ONCO_REPEAT_COUNT + 1).build();

        ReportablePredicate oncoPredicate = new ReportablePredicate(DriverCategory.ONCO, genePanel);

        assertTrue(oncoPredicate.test(variant));
        assertFalse(oncoPredicate.test(variantLargeRepeatCount));
    }
}
