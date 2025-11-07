package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.common.driver.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.purple.PurpleCommon.DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO;
import static com.hartwig.hmftools.common.purple.PurpleCommon.DEFAULT_DRIVER_HET_DELETION_THRESHOLD;
import static com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo.VAR_IMPACT_OTHER_REPORT_DELIM;

import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.driver.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.driver.panel.ReportablePredicate;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

import org.junit.Test;

public class ReportablePredicateTest
{
    private final String GENE_AR = "AR";

    @Test
    public void testAlternateTranscriptImpact()
    {
        AltTranscriptReportableInfo altTransInfo1 = new AltTranscriptReportableInfo(
                "GENE_01", "TRANS_03", "", "", VariantEffect.INTRONIC.effect(), CodingEffect.NONE);

        AltTranscriptReportableInfo altTransInfo2 = new AltTranscriptReportableInfo(
                "GENE_01", "TRANS_02", "", "", VariantEffect.MISSENSE.effect(), CodingEffect.MISSENSE);

        String altTransInfo = altTransInfo1.serialise() + VAR_IMPACT_OTHER_REPORT_DELIM + altTransInfo2.serialise();

        VariantImpact impact = new VariantImpact(
                GENE_AR, "TRANS_01", VariantEffect.INTRONIC.effect(), CodingEffect.NONE, "",
                "", false, altTransInfo, CodingEffect.MISSENSE, 1);

        List<DriverGene> driverGenes = Lists.newArrayList();

        driverGenes.add(ImmutableDriverGene.builder()
                .gene(GENE_AR)
                .reportMissenseAndInframe(true)
                .reportNonsenseAndFrameshift(false)
                .reportSplice(false)
                .reportDeletion(false)
                .reportDisruption(false)
                .reportAmplification(true)
                .amplificationRatio(DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO)
                .reportHetDeletion(true)
                .hetDeletionThreshold(DEFAULT_DRIVER_HET_DELETION_THRESHOLD)
                .reportSomaticHotspot(true)
                .likelihoodType(ONCO)
                .reportGermlineDisruption(DriverGeneGermlineReporting.NONE)
                .reportGermlineDeletion(DriverGeneGermlineReporting.NONE)
                .reportGermlineHotspot(DriverGeneGermlineReporting.NONE)
                .reportGermlineVariant(DriverGeneGermlineReporting.NONE)
                .reportPGX(false)
                .build());

        ReportablePredicate predicate = new ReportablePredicate(ONCO, driverGenes);

        assertTrue(predicate.isReportable(impact, VariantType.SNP, false));
    }
}
