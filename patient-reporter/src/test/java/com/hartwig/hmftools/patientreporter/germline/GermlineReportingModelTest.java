package com.hartwig.hmftools.patientreporter.germline;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.variant.ImmutableReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantSource;
import com.hartwig.hmftools.common.variant.VariantTestFactory;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.patientreporter.PatientReporterTestFactory;
import com.hartwig.hmftools.patientreporter.algo.AnalysedReportData;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GermlineReportingModelTest {

    @Test
    public void canDetermineNotifyForGermlineVariants() {
        LimsGermlineReportingLevel germlineReportingLevel = LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION;
        AnalysedReportData testReportData = PatientReporterTestFactory.loadTestAnalysedReportData();

        ReportableVariant reportableVariant1 = testReportableVariant("MUTYH", GenotypeStatus.HOM_ALT);
        Set<String> germlineGenesWithIndependentHits1 = Sets.newHashSet();

        assertTrue(testReportData.germlineReportingModel()
                .notifyGermlineVariant(reportableVariant1, germlineReportingLevel, germlineGenesWithIndependentHits1));

        ReportableVariant reportableVariant = testReportableVariant("MUTYH", GenotypeStatus.UNKNOWN);
        Set<String> germlineGenesWithIndependentHits = Sets.newHashSet();

        assertFalse(testReportData.germlineReportingModel()
                .notifyGermlineVariant(reportableVariant, germlineReportingLevel, germlineGenesWithIndependentHits));

        assertFalse(testReportData.germlineReportingModel()
                .notifyGermlineVariant(reportableVariant,
                        LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                        germlineGenesWithIndependentHits));
    }

    @NotNull
    private static ReportableVariant testReportableVariant(@NotNull String gene, @NotNull GenotypeStatus genotypeStatus) {
        return ImmutableReportableVariant.builder()
                .from(VariantTestFactory.createTestReportableVariant())
                .type(VariantType.SNP)
                .source(ReportableVariantSource.GERMLINE)
                .gene(gene)
                .genotypeStatus(genotypeStatus)
                .build();
    }
}