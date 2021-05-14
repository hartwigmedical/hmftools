package com.hartwig.hmftools.patientreporter.algo;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.patientreporter.PatientReporterConfig;
import com.hartwig.hmftools.patientreporter.PatientReporterTestFactory;
import com.hartwig.hmftools.protect.purple.ImmutableReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariantSource;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class GenomicAnalyzerTest {

    @Test
    public void canRunOnTestRun() throws IOException {
        AnalysedReportData testReportData = PatientReporterTestFactory.loadTestAnalysedReportData();

        GenomicAnalyzer analyzer = new GenomicAnalyzer(testReportData.germlineReportingModel(),
                testReportData.virusDbModel(),
                testReportData.virusSummaryModel(),
                testReportData.virusBlackListModel());

        PatientReporterConfig config = PatientReporterTestFactory.createTestReporterConfig();

        assertNotNull(analyzer.run("sample", config, LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION));
    }

    @Test
    public void canTestHasOtherGermlineVariantWithDifferentPhaseSet() {
        List<ReportableVariant> reportableVariants1 =
                testReportableVariants("MUTYH", GenotypeStatus.HET, null, "MUTYH", GenotypeStatus.HET, null);
        ReportableVariant reportableVariantToCompare1 = testReportableVariant("MUTYH", GenotypeStatus.HET, null);
        assertFalse(GenomicAnalyzer.hasOtherGermlineVariantWithDifferentPhaseSet(reportableVariants1, reportableVariantToCompare1));

        List<ReportableVariant> reportableVariants2 =
                testReportableVariants("MUTYH", GenotypeStatus.HET, null, "MUTYH", GenotypeStatus.HET, 123);
        ReportableVariant reportableVariantToCompare2 = testReportableVariant("MUTYH", GenotypeStatus.HET, 123);
        assertTrue(GenomicAnalyzer.hasOtherGermlineVariantWithDifferentPhaseSet(reportableVariants2, reportableVariantToCompare2));

    }

    @NotNull
    public List<ReportableVariant> testReportableVariants(@NotNull String gene1, @NotNull GenotypeStatus genotypeStatus1,
            @Nullable Integer localPhaseSet1, @NotNull String gene2, @NotNull GenotypeStatus genotypeStatus2,
            @Nullable Integer localPhaseSet2) {

        ReportableVariant variant1 = ImmutableReportableVariant.builder()
                .type(VariantType.SNP)
                .source(ReportableVariantSource.GERMLINE)
                .gene(gene1)
                .genotypeStatus(genotypeStatus1)
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .canonicalTranscript(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.UNDEFINED)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .totalReadCount(0)
                .alleleReadCount(0)
                .totalCopyNumber(0)
                .alleleCopyNumber(0D)
                .hotspot(Hotspot.HOTSPOT)
                .clonalLikelihood(1D)
                .driverLikelihood(0D)
                .biallelic(false)
                .localPhaseSet(localPhaseSet1)
                .build();

        ReportableVariant variant2 = ImmutableReportableVariant.builder()
                .type(VariantType.SNP)
                .source(ReportableVariantSource.GERMLINE)
                .gene(gene2)
                .genotypeStatus(genotypeStatus2)
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .canonicalTranscript(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.UNDEFINED)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .totalReadCount(0)
                .alleleReadCount(0)
                .totalCopyNumber(0)
                .alleleCopyNumber(0D)
                .hotspot(Hotspot.HOTSPOT)
                .clonalLikelihood(1D)
                .driverLikelihood(0D)
                .biallelic(false)
                .localPhaseSet(localPhaseSet2)
                .build();

        return Lists.newArrayList(variant1, variant2);
    }

    @NotNull
    public ReportableVariant testReportableVariant(@NotNull String gene, @NotNull GenotypeStatus genotypeStatus,
            @Nullable Integer localPhaseSet) {

        return ImmutableReportableVariant.builder()
                .type(VariantType.SNP)
                .source(ReportableVariantSource.GERMLINE)
                .gene(gene)
                .genotypeStatus(genotypeStatus)
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .canonicalTranscript(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.UNDEFINED)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .totalReadCount(0)
                .alleleReadCount(0)
                .totalCopyNumber(0)
                .alleleCopyNumber(0D)
                .hotspot(Hotspot.HOTSPOT)
                .clonalLikelihood(1D)
                .driverLikelihood(0D)
                .biallelic(false)
                .localPhaseSet(localPhaseSet)
                .build();
    }
}
