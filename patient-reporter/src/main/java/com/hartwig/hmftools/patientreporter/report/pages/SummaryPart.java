package com.hartwig.hmftools.patientreporter.report.pages;

import static com.hartwig.hmftools.patientreporter.report.Commons.dataTableStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.fontStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.tableHeaderStyle;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;

import java.text.DecimalFormat;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.report.components.ChordSection;
import com.hartwig.hmftools.patientreporter.report.components.MicrosatelliteSection;
import com.hartwig.hmftools.patientreporter.report.components.MutationalBurdenSection;
import com.hartwig.hmftools.patientreporter.report.components.MutationalLoadSection;
import com.hartwig.hmftools.patientreporter.report.data.TumorReportedGenomicAlterationsDataSource;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;

public final class SummaryPart {

    @NotNull
    public static ComponentBuilder<?, ?> summaryData(@NotNull AnalysedPatientReport report) {
        String avgTumorPloidyString = new DecimalFormat("#.##").format(report.averageTumorPloidy());

        final ComponentBuilder<?, ?> summaryReportedGenomicAlterations =
                cmp.verticalList(cmp.horizontalList(cmp.text("Summary implied tumor characteristics")
                                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)),
                        cmp.horizontalList(cmp.text(
                                "Genomic alterations (somatic variants, copy number changes, gene disruptions and gene fusions). Details can ve found on page 2.")
                                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                                .setStyle(fontStyle().setFontSize(8))),
                        cmp.verticalGap(10),
                        cmp.horizontalList(cmp.text("Tumor Purity").setStyle(tableHeaderStyle()),
                                cmp.text(EvidencePage.impliedPurityString(report)).setStyle(dataTableStyle())),
                        cmp.horizontalList(cmp.text("Average tumor ploidy").setStyle(tableHeaderStyle()),
                                cmp.text(PatientReportFormat.correctValueForFitReliability(avgTumorPloidyString,
                                        report.hasReliablePurityFit())).setStyle(dataTableStyle())),
                        cmp.horizontalList(cmp.text("Tumor mutational load").setStyle(tableHeaderStyle()),
                                cmp.text(MutationalLoadSection.interpret(report.tumorMutationalLoad(), report.hasReliablePurityFit()))
                                        .setStyle(dataTableStyle())),
                        cmp.horizontalList(cmp.text("Tumor mutational burden").setStyle(tableHeaderStyle()),
                                cmp.text(MutationalBurdenSection.interpret(report.tumorMutationalBurden(), report.hasReliablePurityFit()))
                                        .setStyle(dataTableStyle())),
                        cmp.horizontalList(cmp.text("Microsatellite (in)stability").setStyle(tableHeaderStyle()),
                                cmp.text(MicrosatelliteSection.interpret(report.microsatelliteIndelsPerMb(), report.hasReliablePurityFit()))
                                        .setStyle(dataTableStyle())),
                        cmp.horizontalList(cmp.text("BRCAness signature").setStyle(tableHeaderStyle()),
                                cmp.text(ChordSection.interpret(report.chordAnalysis().hrdValue(), report.hasReliablePurityFit()))
                                        .setStyle(dataTableStyle())));

        final ComponentBuilder<?, ?> summaryImpliedTumorCharacteristics =
                cmp.verticalList(cmp.horizontalList(cmp.text("Summary reported genomic alterations")
                                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)),
                        cmp.horizontalList(cmp.text("Whole genome sequencing based tumor characteristics. Details can be found on page 2. ")
                                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                                .setStyle(fontStyle().setFontSize(8))),
                        cmp.verticalGap(10),
                        cmp.text("SomaticVariants").setStyle(tableHeaderStyle()),
                        cmp.horizontalList(cmp.text("Number").setStyle(tableHeaderStyle()),
                                cmp.text(TumorReportedGenomicAlterationsDataSource.countSomaticVariants(report.somaticVariants()))
                                        .setStyle(dataTableStyle())),
                        cmp.horizontalList(cmp.text("Genes").setStyle(tableHeaderStyle()),
                                cmp.text(TumorReportedGenomicAlterationsDataSource.somaticVariantsWithDriver(report.somaticVariants()))
                                        .setStyle(dataTableStyle())),
                        cmp.text("Genes with copy number amplifications").setStyle(tableHeaderStyle()),
                        cmp.text(TumorReportedGenomicAlterationsDataSource.amplificationGenes(report.geneCopyNumbers()))
                                .setStyle(dataTableStyle()),
                        cmp.text("Genes with copy number losses").setStyle(tableHeaderStyle()),
                        cmp.text(TumorReportedGenomicAlterationsDataSource.lossGenes(report.geneCopyNumbers())).setStyle(dataTableStyle()),
                        cmp.text("Genes with fusions").setStyle(tableHeaderStyle()),
                        cmp.text(TumorReportedGenomicAlterationsDataSource.geneFusions(report.geneFusions())).setStyle(dataTableStyle()),
                        cmp.text("Genes with disruptions").setStyle(tableHeaderStyle()),
                        cmp.text(TumorReportedGenomicAlterationsDataSource.geneDisruptions(report.geneDisruptions()))
                                .setStyle(dataTableStyle()));

        return cmp.horizontalList(summaryImpliedTumorCharacteristics,
                cmp.horizontalGap(40),
                cmp.verticalGap(3),
                summaryReportedGenomicAlterations);
    }
}
