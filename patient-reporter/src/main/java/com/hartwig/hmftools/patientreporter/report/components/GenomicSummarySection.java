package com.hartwig.hmftools.patientreporter.report.components;

import static com.hartwig.hmftools.patientreporter.report.Commons.dataTableStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.fontStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.summaryHeaderStyle;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;

import java.text.DecimalFormat;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.report.data.GeneCopyNumberDataSource;
import com.hartwig.hmftools.patientreporter.report.data.GeneFusionDataSource;
import com.hartwig.hmftools.patientreporter.report.data.SomaticVariantDataSource;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneFusion;
import com.hartwig.hmftools.patientreporter.variants.ReportableSomaticVariant;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;

public final class GenomicSummarySection {

    private GenomicSummarySection() {
    }

    @NotNull
    public static ComponentBuilder<?, ?> build(@NotNull AnalysedPatientReport report) {
        String avgTumorPloidyString = new DecimalFormat("#.#").format(report.averageTumorPloidy());
        int genomicAlterationHeaderWidth = 105;

        ComponentBuilder<?, ?> genomicAlterationSummary = cmp.verticalList(cmp.horizontalList(cmp.text("Genomic Alterations")
                        .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)),
                cmp.horizontalList(cmp.text(
                        "Summary on genomic alterations (somatic variants, copy number changes, gene disruptions and gene fusions). "
                                + "Details can be found further down the report.")
                        .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                        .setStyle(fontStyle().setFontSize(8))),
                cmp.verticalGap(10),
                cmp.horizontalList(cmp.text("Genes with driver variant")
                                .setStyle(summaryHeaderStyle())
                                .setFixedWidth(genomicAlterationHeaderWidth),
                        cmp.text(somaticVariantsWithDriver(report.somaticVariants())).setStyle(dataTableStyle())),
                cmp.horizontalList(cmp.text("Number of reported variants")
                                .setStyle(summaryHeaderStyle())
                                .setFixedWidth(genomicAlterationHeaderWidth),
                        cmp.text(countSomaticVariants(report.somaticVariants())).setStyle(dataTableStyle())),
                cmp.horizontalList(cmp.text("Genes with copy-gain")
                                .setStyle(summaryHeaderStyle())
                                .setFixedWidth(genomicAlterationHeaderWidth),
                        cmp.text(amplificationGenes(report.geneCopyNumbers())).setStyle(dataTableStyle())),
                cmp.horizontalList(cmp.text("Genes with copy-loss")
                                .setStyle(summaryHeaderStyle())
                                .setFixedWidth(genomicAlterationHeaderWidth),
                        cmp.text(lossGenes(report.geneCopyNumbers())).setStyle(dataTableStyle())),
                cmp.horizontalList(cmp.text("Gene fusions").setStyle(summaryHeaderStyle()).setFixedWidth(genomicAlterationHeaderWidth),
                        cmp.text(geneFusions(report.geneFusions())).setStyle(dataTableStyle())));

        ComponentBuilder<?, ?> tumorCharacteristicsPart = cmp.verticalList(cmp.horizontalList(cmp.text("Tumor Characteristics")
                        .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)),
                cmp.horizontalList(cmp.text(
                        "Whole genome sequencing based tumor characteristics. Details can be found down further the report.")
                        .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                        .setStyle(fontStyle().setFontSize(8))),
                cmp.verticalGap(10),
                cmp.horizontalList(cmp.text("Tumor purity").setStyle(summaryHeaderStyle()),
                        cmp.text(impliedPurityString(report)).setStyle(dataTableStyle())),
                cmp.horizontalList(cmp.text("Average tumor ploidy").setStyle(summaryHeaderStyle()),
                        cmp.text(PatientReportFormat.correctValueForFitReliability(avgTumorPloidyString, report.hasReliablePurityFit()))
                                .setStyle(dataTableStyle())),
                cmp.horizontalList(cmp.text("Tumor mutational load").setStyle(summaryHeaderStyle()),
                        cmp.text(PatientReportFormat.correctValueForFitReliability(interpretMutationalLoad(report.tumorMutationalLoad()),
                                report.hasReliablePurityFit())).setStyle(dataTableStyle())),
                cmp.horizontalList(cmp.text("Microsatellite (in)stability").setStyle(summaryHeaderStyle()),
                        cmp.text(PatientReportFormat.correctValueForFitReliability(interpretMSI(report.microsatelliteIndelsPerMb()),
                                report.hasReliablePurityFit())).setStyle(dataTableStyle())));

        return cmp.horizontalList(genomicAlterationSummary, cmp.horizontalGap(40), tumorCharacteristicsPart);
    }

    @NotNull
    private static String impliedPurityString(@NotNull AnalysedPatientReport report) {
        return report.hasReliablePurityFit() ? PatientReportFormat.formatPercent(report.impliedPurity()) : "[below detection threshold]";
    }

    @NotNull
    private static String interpretMutationalLoad(int mutationalLoad) {
        if (mutationalLoad > MutationalLoadSection.ML_THRESHOLD) {
            return "High";
        } else {
            return "Low";
        }
    }

    @NotNull
    private static String interpretMSI(double microsatelliteIndicator) {
        if (microsatelliteIndicator > MicrosatelliteSection.MSI_THRESHOLD) {
            return "Unstable";
        } else {
            return "Stable";
        }
    }

    @NotNull
    private static String countSomaticVariants(@NotNull List<ReportableSomaticVariant> variants) {
        return Integer.toString(variants.size());
    }

    @NotNull
    private static String somaticVariantsWithDriver(@NotNull List<ReportableSomaticVariant> variants) {
        List<String> somaticVariants = Lists.newArrayList();
        for (final ReportableSomaticVariant variant : SomaticVariantDataSource.sort(variants)) {
            if (SomaticVariantDataSource.driverField(variant).equals("High")) {
                somaticVariants.add(variant.gene());
            }
        }
        return String.join(", ", somaticVariants);
    }

    @NotNull
    private static String amplificationGenes(@NotNull List<GeneCopyNumber> copyNumbers) {
        List<String> geneCopyNumbersAmplification = Lists.newArrayList();
        for (GeneCopyNumber copyNumber : copyNumbers) {
            if (GeneCopyNumberDataSource.type(copyNumber).equals("gain")) {
                geneCopyNumbersAmplification.add(copyNumber.gene());
            }
        }
        return String.join(", ", geneCopyNumbersAmplification);
    }

    @NotNull
    private static String lossGenes(@NotNull List<GeneCopyNumber> copyNumbers) {
        List<String> geneCopyNumbersLoss = Lists.newArrayList();
        for (GeneCopyNumber copyNumber : copyNumbers) {
            if (GeneCopyNumberDataSource.type(copyNumber).equals("full loss") || GeneCopyNumberDataSource.type(copyNumber)
                    .equals("partial loss")) {
                geneCopyNumbersLoss.add(copyNumber.gene());
            }
        }
        return String.join(", ", geneCopyNumbersLoss);
    }

    @NotNull
    private static String geneFusions(@NotNull List<ReportableGeneFusion> fusions) {
        List<String> geneFusions = Lists.newArrayList();
        for (ReportableGeneFusion fusion : GeneFusionDataSource.sort(fusions)) {
            geneFusions.add(GeneFusionDataSource.name(fusion));
        }
        return String.join(", ", geneFusions);
    }
}
