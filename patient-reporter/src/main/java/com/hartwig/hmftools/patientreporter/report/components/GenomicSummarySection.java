package com.hartwig.hmftools.patientreporter.report.components;

import static com.hartwig.hmftools.patientreporter.report.Commons.dataTableStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.fontStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.summaryHeaderStyle;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;

import java.text.DecimalFormat;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.report.data.TumorImpliedTumorCharacteristicsDataSource;
import com.hartwig.hmftools.patientreporter.report.data.TumorReportedGenomicAlterationsDataSource;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

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
                        cmp.text(TumorReportedGenomicAlterationsDataSource.somaticVariantsWithDriver(report.somaticVariants()))
                                .setStyle(dataTableStyle())),
                cmp.horizontalList(cmp.text("Number of reported variants")
                                .setStyle(summaryHeaderStyle())
                                .setFixedWidth(genomicAlterationHeaderWidth),
                        cmp.text(TumorReportedGenomicAlterationsDataSource.countSomaticVariants(report.somaticVariants()))
                                .setStyle(dataTableStyle())),
                cmp.horizontalList(cmp.text("Genes with copy-gain")
                                .setStyle(summaryHeaderStyle())
                                .setFixedWidth(genomicAlterationHeaderWidth),
                        cmp.text(TumorReportedGenomicAlterationsDataSource.amplificationGenes(report.geneCopyNumbers()))
                                .setStyle(dataTableStyle())),
                cmp.horizontalList(cmp.text("Genes with copy-loss")
                                .setStyle(summaryHeaderStyle())
                                .setFixedWidth(genomicAlterationHeaderWidth),
                        cmp.text(TumorReportedGenomicAlterationsDataSource.lossGenes(report.geneCopyNumbers())).setStyle(dataTableStyle())),
                cmp.horizontalList(cmp.text("Gene fusions").setStyle(summaryHeaderStyle()).setFixedWidth(genomicAlterationHeaderWidth),
                        cmp.text(TumorReportedGenomicAlterationsDataSource.geneFusions(report.geneFusions())).setStyle(dataTableStyle())));

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
                        cmp.text(PatientReportFormat.correctValueForFitReliability(TumorImpliedTumorCharacteristicsDataSource.interpretMutationalLoad(
                                report.tumorMutationalLoad()), report.hasReliablePurityFit())).setStyle(dataTableStyle())),
                cmp.horizontalList(cmp.text("Microsatellite (in)stability").setStyle(summaryHeaderStyle()),
                        cmp.text(PatientReportFormat.correctValueForFitReliability(TumorImpliedTumorCharacteristicsDataSource.interpretMSI(
                                report.microsatelliteIndelsPerMb()),report.hasReliablePurityFit())).setStyle(dataTableStyle())));

        return cmp.horizontalList(genomicAlterationSummary, cmp.horizontalGap(40), tumorCharacteristicsPart);
    }

    @NotNull
    private static String impliedPurityString(@NotNull AnalysedPatientReport report) {
        return report.hasReliablePurityFit() ? PatientReportFormat.formatPercent(report.impliedPurity()) : "[below detection threshold]";
    }
}
