package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.drivercatalog.panel.ReportablePredicate.MAX_ONCO_REPEAT_COUNT;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo.VAR_IMPACT_OTHER_REPORT_DELIM;

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
import com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

import org.junit.Test;

public class ReportablePredicateTest
{
    private final String GENE_AR = "AR";
    private final DriverGenePanel genePanel = loadTestPanel();

    @Test
    public void testAlternateTranscriptImpact()
    {
        AltTranscriptReportableInfo altTransInfo1 = new AltTranscriptReportableInfo(
                "GENE_01", "TRANS_03", "", "", VariantEffect.INTRONIC.effect(), CodingEffect.NONE);

        AltTranscriptReportableInfo altTransInfo2 = new AltTranscriptReportableInfo(
                "GENE_01","TRANS_02", "", "", VariantEffect.MISSENSE.effect(), CodingEffect.MISSENSE);

        String altTransInfo = altTransInfo1.serialise() + VAR_IMPACT_OTHER_REPORT_DELIM + altTransInfo2.serialise();

        VariantImpact impact = new VariantImpact(
                GENE_AR, "TRANS_01", VariantEffect.INTRONIC.effect(), CodingEffect.NONE, "",
                "", false, altTransInfo, CodingEffect.MISSENSE, 1);

        ReportablePredicate predicate = new ReportablePredicate(ONCO, genePanel.driverGenes());

        assertTrue(predicate.isReportable(impact, VariantType.SNP, 0, false));
    }

    @Test
    public void testIgnoreIndelsWithLargeRepeatCount()
    {
        final SomaticVariant variant = SomaticVariantTestFactory.builder()
                .gene(GENE_AR)
                .repeatCount(MAX_ONCO_REPEAT_COUNT)
                .type(VariantType.INDEL)
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .build();

        final SomaticVariant variantLargeRepeatCount =
                ImmutableSomaticVariantImpl.builder().from(variant).repeatCount(MAX_ONCO_REPEAT_COUNT + 1).build();

        ReportablePredicate oncoPredicate = new ReportablePredicate(ONCO, genePanel.driverGenes());

        assertTrue(oncoPredicate.isReportable(
                variant.gene(), variant.type(), variant.repeatCount(), variant.isHotspot(),
                variant.canonicalCodingEffect(), variant.canonicalEffect()));

        assertFalse(oncoPredicate.isReportable(
                variantLargeRepeatCount.gene(), variantLargeRepeatCount.type(), variantLargeRepeatCount.repeatCount(),
                variantLargeRepeatCount.isHotspot(), variantLargeRepeatCount.canonicalCodingEffect(), variantLargeRepeatCount.canonicalEffect()));
    }

    private DriverGenePanel loadTestPanel()
    {
        List<DriverGene> driverGenes = Lists.newArrayList();

        driverGenes.add(ImmutableDriverGene.builder()
                .gene(GENE_AR)
                .reportMissenseAndInframe(true)
                .reportNonsenseAndFrameshift(false)
                .reportSplice(false)
                .reportDeletion(false)
                .reportDisruption(false)
                .reportAmplification(true)
                .reportSomaticHotspot(true)
                .likelihoodType(ONCO)
                .reportGermlineDisruption(DriverGeneGermlineReporting.NONE)
                .reportGermlineDeletion(DriverGeneGermlineReporting.NONE)
                .reportGermlineHotspot(DriverGeneGermlineReporting.NONE)
                .reportGermlineVariant(DriverGeneGermlineReporting.NONE)
                .reportPGX(false)
                .build());

        return DriverGenePanelFactory.create(driverGenes);
    }
}
