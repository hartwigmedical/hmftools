package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.drivercatalog.panel.ReportablePredicate.MAX_ONCO_REPEAT_COUNT;

import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelFactory;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.ReportablePredicate;
import com.hartwig.hmftools.common.test.SomaticVariantTestFactory;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.junit.Test;

public class ReportablePredicateTest {

    private final DriverGenePanel genePanel = loadTestPanel();

    @Test
    public void testIgnoreIndelsWithLargeRepeatCount() {
        final SomaticVariant variant = SomaticVariantTestFactory.builder()
                .gene("AR")
                .repeatCount(MAX_ONCO_REPEAT_COUNT)
                .type(VariantType.INDEL)
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .build();

        final SomaticVariant variantLargeRepeatCount =
                ImmutableSomaticVariantImpl.builder().from(variant).repeatCount(MAX_ONCO_REPEAT_COUNT + 1).build();

        ReportablePredicate oncoPredicate = new ReportablePredicate(ONCO, genePanel.driverGenes());

        assertTrue(oncoPredicate.test(
                variant.gene(), variant.type(), variant.repeatCount(), variant.isHotspot(),
                variant.canonicalCodingEffect(), variant.canonicalEffect()));

        assertFalse(oncoPredicate.test(
                variantLargeRepeatCount.gene(), variantLargeRepeatCount.type(), variantLargeRepeatCount.repeatCount(),
                variantLargeRepeatCount.isHotspot(), variantLargeRepeatCount.canonicalCodingEffect(), variantLargeRepeatCount.canonicalEffect()));
    }

    private DriverGenePanel loadTestPanel()
    {
        List<DriverGene> driverGenes = Lists.newArrayList();

        driverGenes.add(ImmutableDriverGene.builder()
                .gene("AR")
                .reportMissenseAndInframe(true)
                .reportNonsenseAndFrameshift(false)
                .reportSplice(false)
                .reportDeletion(false)
                .reportDisruption(false)
                .reportAmplification(true)
                .reportSomaticHotspot(true)
                .likelihoodType(ONCO)
                .reportGermlineDisruption(false)
                .reportGermlineHotspot(DriverGeneGermlineReporting.NONE)
                .reportGermlineVariant(DriverGeneGermlineReporting.NONE)
                .reportGermlineDisruption(false)
                .build());

        return DriverGenePanelFactory.create(driverGenes);
    }
}
