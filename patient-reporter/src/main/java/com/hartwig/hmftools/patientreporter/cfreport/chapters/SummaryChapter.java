package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.ClinicalTrial;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.cfreport.MathUtil;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.InlineBarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.LineDivider;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.components.TumorLocationAndTypeTable;
import com.hartwig.hmftools.patientreporter.cfreport.data.ClinicalTrials;
import com.hartwig.hmftools.patientreporter.cfreport.data.DataUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.EvidenceItems;
import com.hartwig.hmftools.patientreporter.cfreport.data.GainsAndLosses;
import com.hartwig.hmftools.patientreporter.cfreport.data.GeneFusions;
import com.hartwig.hmftools.patientreporter.cfreport.data.GeneUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.HrDeficiency;
import com.hartwig.hmftools.patientreporter.cfreport.data.SomaticVariants;
import com.hartwig.hmftools.patientreporter.cfreport.data.TumorPurity;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.Style;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;
import com.itextpdf.layout.property.VerticalAlignment;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SummaryChapter implements ReportChapter {

    private static final float TABLE_SPACER_HEIGHT = 5;

    @NotNull
    private final AnalysedPatientReport patientReport;

    public SummaryChapter(@NotNull final AnalysedPatientReport patientReport) {
        this.patientReport = patientReport;
    }

    @NotNull
    @Override
    public String name() {
        return patientReport.isCorrectedReport() ? "DNA Analysis Report (Corrected)" : "DNA Analysis Report";
    }

    @Override
    public boolean isFullWidth() {
        return false;
    }

    @Override
    public boolean hasCompleteSidebar() {
        return true;
    }

    @Override
    public final void render(@NotNull Document reportDocument) {
        reportDocument.add(new Paragraph("Summary").addStyle(ReportResources.chapterTitleStyle()));
        reportDocument.add(TumorLocationAndTypeTable.createTumorLocationAndType(patientReport.sampleReport().primaryTumorLocationString(),
                patientReport.sampleReport().cancerSubTypeString(),
                contentWidth()));
        reportDocument.add(new Paragraph("\nThe information regarding 'primary tumor location' and 'cancer subtype'  is based on "
                + "information received \nfrom the originating hospital.").addStyle(ReportResources.subTextStyle()));

        renderSummaryText(reportDocument, patientReport.clinicalSummary());
        renderTreatmentIndications(reportDocument, patientReport.tumorSpecificEvidence(), patientReport.clinicalTrials());
        renderTumorCharacteristics(reportDocument, patientReport);
        renderGenomicAlterations(reportDocument, patientReport);
    }

    private void renderSummaryText(@NotNull final Document reportDocument, @Nullable final String text) {
        if (text == null || text.isEmpty()) {
            return;
        }

        Div div = createSectionStartDiv(contentWidth());
        div.add(new Paragraph("Conclusion").addStyle(ReportResources.sectionTitleStyle()));

        div.add(new Paragraph(text).setWidth(contentWidth()).addStyle(ReportResources.bodyTextStyle()).setFixedLeading(11));

        reportDocument.add(div);
    }

    private void renderTreatmentIndications(@NotNull Document reportDocument, @NotNull List<EvidenceItem> tumorSpecificEvidence,
            @NotNull List<ClinicalTrial> trials) {
        Div div = createSectionStartDiv(contentWidth());

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell()
                .add(new Paragraph("Treatment indications (tumor-type specific)").addStyle(ReportResources.sectionTitleStyle())));

        table.addCell(TableUtil.createLayoutCell(1, 2).setHeight(TABLE_SPACER_HEIGHT));

        int therapyEventCount = EvidenceItems.uniqueEventCount(tumorSpecificEvidence);
        int therapyCount = EvidenceItems.uniqueTherapyCount(tumorSpecificEvidence);
        table.addCell(createMiddleAlignedCell().add(new Paragraph("Number of alterations with therapy indication").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createTreatmentIndicationCell(therapyEventCount, therapyCount, "treatment(s)"));

        int trialEventCount = ClinicalTrials.uniqueEventCount(trials);
        int trialCount = ClinicalTrials.uniqueTrialCount(trials);
        table.addCell(createMiddleAlignedCell().add(new Paragraph("Number of alterations with clinical trial eligibility").addStyle(
                ReportResources.bodyTextStyle())));
        table.addCell(createTreatmentIndicationCell(trialEventCount, trialCount, "trial(s)"));

        div.add(table);

        reportDocument.add(div);
    }

    private void renderTumorCharacteristics(@NotNull Document reportDocument, @NotNull AnalysedPatientReport patientReport) {
        boolean hasReliablePurity = patientReport.hasReliablePurity();

        Div div = createSectionStartDiv(contentWidth());

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, .33f, .66f }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell()
                .add(new Paragraph("Tumor characteristics").addStyle(ReportResources.sectionTitleStyle())));
        table.addCell(TableUtil.createLayoutCell(1, 3).setHeight(TABLE_SPACER_HEIGHT));

        double impliedPurity = patientReport.impliedPurity();
        double impliedPurityPercentage = MathUtil.mapPercentage(impliedPurity, TumorPurity.RANGE_MIN, TumorPurity.RANGE_MAX);
        renderTumorPurity(hasReliablePurity,
                DataUtil.formatPercentage(impliedPurityPercentage),
                impliedPurity,
                TumorPurity.RANGE_MIN,
                TumorPurity.RANGE_MAX,
                table);

        Style dataStyle = hasReliablePurity ? ReportResources.dataHighlightStyle() : ReportResources.dataHighlightNaStyle();

        table.addCell(createMiddleAlignedCell().add(new Paragraph("Average tumor ploidy").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell(2).add(createHighlightParagraph(GeneUtil.ploidyToCopiesString(patientReport.averageTumorPloidy(),
                hasReliablePurity)).addStyle(dataStyle)));

        String mutationalLoadString =
                hasReliablePurity ? patientReport.tumorMutationalLoadStatus().display() : DataUtil.NA_STRING;
        table.addCell(createMiddleAlignedCell().add(new Paragraph("Tumor mutational load").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell(2).add(createHighlightParagraph(mutationalLoadString).addStyle(dataStyle)));

        String microSatelliteStabilityString =
                hasReliablePurity ? patientReport.microsatelliteStatus().display() : DataUtil.NA_STRING;
        table.addCell(createMiddleAlignedCell().add(new Paragraph("Microsatellite (in)stability").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell(2).add(createHighlightParagraph(microSatelliteStabilityString).addStyle(dataStyle)));
        div.add(table);

        String hrdString = hasReliablePurity ? HrDeficiency.interpretChordStatusToString(patientReport.chordAnalysis().hrdValue()) : DataUtil.NA_STRING;
        table.addCell(createMiddleAlignedCell().add(new Paragraph("HR Status").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell(2).add(createHighlightParagraph(hrdString).addStyle(dataStyle)));

        reportDocument.add(div);
    }

    private static void renderTumorPurity(boolean hasReliablePurity, @NotNull String valueLabel, double value, double min, double max,
            @NotNull Table table) {
        final String label = "Tumor purity of biopsy";
        table.addCell(createMiddleAlignedCell().add(new Paragraph(label).addStyle(ReportResources.bodyTextStyle())));

        if (hasReliablePurity) {
            table.addCell(createMiddleAlignedCell().add(createHighlightParagraph(valueLabel).addStyle(ReportResources.dataHighlightStyle())));
            table.addCell(createMiddleAlignedCell().add(createInlineBarChart(value, min, max)));
        } else {
            table.addCell(createMiddleAlignedCell(2).add(createHighlightParagraph(DataUtil.BELOW_DETECTION_STRING).addStyle(ReportResources.dataHighlightNaStyle())));
        }
    }

    private void renderGenomicAlterations(@NotNull Document report, @NotNull AnalysedPatientReport patientReport) {
        Div div = createSectionStartDiv(contentWidth());

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell().add(new Paragraph("Genomic alterations").addStyle(ReportResources.sectionTitleStyle())));
        table.addCell(TableUtil.createLayoutCell(1, 2).setHeight(TABLE_SPACER_HEIGHT));

        Set<String> driverVariantGenes = SomaticVariants.driverGenesWithVariant(patientReport.reportableVariants());

        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Driver genes with variant(s)").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(sortGenes(driverVariantGenes)));

        int reportedVariants = SomaticVariants.countReportableVariants(patientReport.reportableVariants());
        Style reportedVariantsStyle =
                (reportedVariants > 0) ? ReportResources.dataHighlightStyle() : ReportResources.dataHighlightNaStyle();
        table.addCell(createMiddleAlignedCell().add(new Paragraph("Number of reported variants").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell().add(createHighlightParagraph(String.valueOf(reportedVariants)).addStyle(
                reportedVariantsStyle)));

        Set<String> amplifiedGenes = GainsAndLosses.amplifiedGenes(patientReport.gainsAndLosses());
        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Genes with copy-gain").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(sortGenes(amplifiedGenes)));

        Set<String> copyLossGenes = GainsAndLosses.lostGenes(patientReport.gainsAndLosses());
        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Genes with copy-loss").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(sortGenes(copyLossGenes)));

        // TODO First do DEV-1013 - disable in final report.
        //        Set<String> disruptedGenes = HomozygousDisruptions.disruptedGenes(patientReport.reportableHomozygousDisruptions());
        //        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
        //                .add(new Paragraph("Disrupted genes").addStyle(ReportResources.bodyTextStyle())));
        //        table.addCell(createGeneListCell(sortGenes(disruptedGenes)));

        Set<String> fusionGenes = GeneFusions.uniqueGeneFusions(patientReport.geneFusions());
        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Gene fusions").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(sortGenes(fusionGenes)));

        div.add(table);

        report.add(div);
    }

    @NotNull
    @VisibleForTesting
    static Set<String> sortGenes(@NotNull Set<String> driverVariantGenes) {
        List<String> genesList = Lists.newArrayList(driverVariantGenes);
        Collections.sort(genesList);
        return Sets.newHashSet(genesList);
    }

    @NotNull
    private static Div createSectionStartDiv(float width) {
        return new Div().setKeepTogether(true).setWidth(width).add(LineDivider.createLineDivider(width));
    }

    @NotNull
    private static Cell createMiddleAlignedCell() {
        return createMiddleAlignedCell(1);
    }

    @NotNull
    private static Cell createMiddleAlignedCell(int colSpan) {
        return TableUtil.createLayoutCell(1, colSpan).setVerticalAlignment(VerticalAlignment.MIDDLE);
    }

    @NotNull
    private static Cell createGeneListCell(@NotNull Set<String> genes) {
        String geneString = (genes.size() > 0) ? String.join(", ", genes) : DataUtil.NONE_STRING;

        Style style = (genes.size() > 0) ? ReportResources.dataHighlightStyle() : ReportResources.dataHighlightNaStyle();

        return createMiddleAlignedCell().add(createHighlightParagraph(geneString)).addStyle(style);
    }

    @NotNull
    private static Paragraph createHighlightParagraph(@NotNull String text) {
        return new Paragraph(text).setFixedLeading(14);
    }

    @NotNull
    private static Cell createTreatmentIndicationCell(int eventCount, int treatmentCount, @NotNull String treatmentsName) {
        String treatmentText;
        Style style;
        if (eventCount > 0) {
            treatmentText = String.format("%d | %d %s", eventCount, treatmentCount, treatmentsName);
            style = ReportResources.dataHighlightStyle();
        } else {
            treatmentText = DataUtil.NONE_STRING;
            style = ReportResources.dataHighlightNaStyle();
        }

        return createMiddleAlignedCell().add(createHighlightParagraph(treatmentText)).addStyle(style);
    }

    @NotNull
    private static InlineBarChart createInlineBarChart(double value, double min, double max) {
        InlineBarChart chart = new InlineBarChart(value, min, max);
        chart.setWidth(41);
        chart.setHeight(6);
        return chart;
    }
}
