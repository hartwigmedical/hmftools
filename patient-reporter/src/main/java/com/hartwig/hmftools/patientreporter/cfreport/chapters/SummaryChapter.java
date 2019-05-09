package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import java.util.List;
import java.util.Set;

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
import com.hartwig.hmftools.patientreporter.cfreport.data.GeneCopyNumbers;
import com.hartwig.hmftools.patientreporter.cfreport.data.GeneFusions;
import com.hartwig.hmftools.patientreporter.cfreport.data.GeneUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.MicroSatelliteStatus;
import com.hartwig.hmftools.patientreporter.cfreport.data.MutationalLoad;
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
        return "Summary";
    }

    @Override
    public String pageNumberPrefix() {
        return name();
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
        reportDocument.add(TumorLocationAndTypeTable.createTumorLocationAndType(patientReport.sampleReport().primaryTumorLocationString(),
                patientReport.sampleReport().cancerSubTypeString(),
                contentWidth()));

        final String summaryContent = patientReport.summarySample();
        renderSummaryText(summaryContent, reportDocument);

        renderTreatmentIndications(patientReport.tumorSpecificEvidence(), patientReport.clinicalTrials(), reportDocument);
        renderTumorCharacteristics(patientReport, reportDocument);
        renderGenomicAlterations(patientReport, reportDocument);
    }

    private void renderSummaryText(@Nullable final String text, @NotNull final Document reportDocument) {
        if (text == null || text.isEmpty()) {
            return;
        }

        Div div = createSectionStartDiv(contentWidth());
        div.add(new Paragraph("Summary").addStyle(ReportResources.sectionTitleStyle()));

        div.add(new Paragraph(text).setWidth(contentWidth()).addStyle(ReportResources.bodyTextStyle()).setFixedLeading(11));

        reportDocument.add(div);
    }

    private void renderTreatmentIndications(@NotNull List<EvidenceItem> tumorSpecificEvidence, @NotNull List<ClinicalTrial> trials,
            @NotNull Document reportDocument) {
        Div div = createSectionStartDiv(contentWidth());

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell()
                .add(new Paragraph("Treatment indications").addStyle(ReportResources.sectionTitleStyle())));

        table.addCell(TableUtil.createLayoutCell()
                .add(new Paragraph("Summary of alterations with number of treatment indication and/or clinical trials").addStyle(
                        ReportResources.bodyTextStyle()).setFixedLeading(ReportResources.BODY_TEXT_LEADING)));
        table.addCell(TableUtil.createLayoutCell(1, 2).setHeight(TABLE_SPACER_HEIGHT));

        int therapyEventCount = EvidenceItems.uniqueEventCount(tumorSpecificEvidence);
        int therapyCount = EvidenceItems.uniqueTherapyCount(tumorSpecificEvidence);
        table.addCell(createMiddleAlignedCell().add(new Paragraph("Gene alteration(s) with therapy indication(s)").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createTreatmentIndicationCell(therapyEventCount, therapyCount, "treatments"));

        int trialEventCount = ClinicalTrials.uniqueEventsCount(trials);
        int trialCount = ClinicalTrials.uniqueTrialCount(trials);
        table.addCell(createMiddleAlignedCell().add(new Paragraph("Gene alteration(s) with clinical trial eligibility").addStyle(
                ReportResources.bodyTextStyle())));
        table.addCell(createTreatmentIndicationCell(trialEventCount, trialCount, "trials"));

        div.add(table);

        reportDocument.add(div);
    }

    private void renderTumorCharacteristics(@NotNull AnalysedPatientReport patientReport, @NotNull Document reportDocument) {
        final boolean hasReliablePurityFit = patientReport.hasReliablePurityFit();

        Div div = createSectionStartDiv(contentWidth());

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, .33f, .66f }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell()
                .add(new Paragraph("Tumor characteristics summary").addStyle(ReportResources.sectionTitleStyle())));
        table.addCell(TableUtil.createLayoutCell(1, 2)
                .add(new Paragraph("Whole genome sequencing based tumor characteristics.").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(TableUtil.createLayoutCell(1, 3).setHeight(TABLE_SPACER_HEIGHT));

        final double impliedPurity = patientReport.impliedPurity();
        final double impliedPurityPercentage = MathUtil.mapPercentage(impliedPurity, TumorPurity.RANGE_MIN, TumorPurity.RANGE_MAX);
        renderTumorPurity(hasReliablePurityFit,
                DataUtil.formatPercentage(impliedPurityPercentage),
                impliedPurity,
                TumorPurity.RANGE_MIN,
                TumorPurity.RANGE_MAX,
                table);

        Style dataStyle = hasReliablePurityFit ? ReportResources.dataHighlightStyle() : ReportResources.dataHighlightNaStyle();

        table.addCell(createMiddleAlignedCell().add(new Paragraph("Average tumor ploidy").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell(2).add(createHighlightParagraph(GeneUtil.ploidyToCopiesString(patientReport.averageTumorPloidy(),
                hasReliablePurityFit)).addStyle(dataStyle)));

        final String mutationalLoadString =
                hasReliablePurityFit ? MutationalLoad.interpretToString(patientReport.tumorMutationalLoad()) : DataUtil.NA_STRING;
        table.addCell(createMiddleAlignedCell().add(new Paragraph("Tumor mutational load").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell(2).add(createHighlightParagraph(mutationalLoadString).addStyle(dataStyle)));

        final String microSatelliteStabilityString = hasReliablePurityFit
                ? MicroSatelliteStatus.interpretToString(patientReport.microsatelliteIndelsPerMb())
                : DataUtil.NA_STRING;
        table.addCell(createMiddleAlignedCell().add(new Paragraph("Microsatellite (in)stability").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell(2).add(createHighlightParagraph(microSatelliteStabilityString).addStyle(dataStyle)));
        div.add(table);
        reportDocument.add(div);
    }

    private static void renderTumorPurity(boolean hasReliablePurityFit, @NotNull String valueLabel, double value, double min, double max,
            @NotNull Table table) {
        final String label = "Tumor purity of biopsy";
        table.addCell(createMiddleAlignedCell().add(new Paragraph(label).addStyle(ReportResources.bodyTextStyle())));

        if (hasReliablePurityFit) {
            table.addCell(createMiddleAlignedCell().add(createHighlightParagraph(valueLabel).addStyle(ReportResources.dataHighlightStyle())));
            table.addCell(createMiddleAlignedCell().add(createInlineBarChart(value, min, max)));
        } else {
            table.addCell(createMiddleAlignedCell(2).add(createHighlightParagraph(DataUtil.NA_STRING).addStyle(ReportResources.dataHighlightNaStyle())));
        }
    }

    private void renderGenomicAlterations(@NotNull AnalysedPatientReport patientReport, @NotNull Document report) {
        final Div div = createSectionStartDiv(contentWidth());

        final Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.createLayoutCell()
                .add(new Paragraph("Genomic alterations \nsummary").addStyle(ReportResources.sectionTitleStyle())));
        table.addCell(TableUtil.createLayoutCell()
                .add(new Paragraph("Summary on genomic alterations "
                        + "(somatic variants, copy number changes, gene disruptions and gene fusions).").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(TableUtil.createLayoutCell(1, 2).setHeight(TABLE_SPACER_HEIGHT));

        final Set<String> driverVariantGenes = SomaticVariants.driverGenesWithVariant(patientReport.reportableVariants());
        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Driver genes with variant").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(driverVariantGenes));

        final int reportedVariants = SomaticVariants.countReportableVariants(patientReport.reportableVariants());
        Style reportedVariantsStyle =
                (reportedVariants > 0) ? ReportResources.dataHighlightStyle() : ReportResources.dataHighlightNaStyle();
        table.addCell(createMiddleAlignedCell().add(new Paragraph("Nr. of reported variants").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell().add(createHighlightParagraph(String.valueOf(reportedVariants)).addStyle(
                reportedVariantsStyle)));

        final Set<String> amplifiedGenes = GeneCopyNumbers.amplifiedGenes(patientReport.geneCopyNumbers());
        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Genes with copy-gain").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(amplifiedGenes));

        final Set<String> copyLossGenes = GeneCopyNumbers.lossGenes(patientReport.geneCopyNumbers());
        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Genes with copy-loss").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(copyLossGenes));

        final Set<String> fusionGenes = GeneFusions.uniqueGeneFusions(patientReport.geneFusions());
        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Gene fusions").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(fusionGenes));

        div.add(table);

        report.add(div);
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
            treatmentText = String.format("%d (%d %s)", eventCount, treatmentCount, treatmentsName);
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
