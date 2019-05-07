package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.common.actionability.ClinicalTrial;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.cfreport.MathUtil;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.InlineBarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.LineDivider;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
import com.hartwig.hmftools.patientreporter.cfreport.components.TumorLocationAndTypeTable;
import com.hartwig.hmftools.patientreporter.cfreport.data.*;
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

import java.util.List;

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

        // @TODO Replace this fixed text with the patientReport.summaryText method.
        // Return value from that method can be null which is gracefully handled by renderSummaryText :

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

    private void renderTreatmentIndications(@NotNull final List<EvidenceItem> tumorSpecificEvidence,
            @NotNull final List<ClinicalTrial> trials, @NotNull Document reportDocument) {
        Div div = createSectionStartDiv(contentWidth());

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, 1 }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.getLayoutCell().add(new Paragraph("Treatment indications").addStyle(ReportResources.sectionTitleStyle())));

        table.addCell(TableUtil.getLayoutCell()
                .add(new Paragraph("Summary of number of alterations with number of treatment indication and/or clinical studies").addStyle(
                        ReportResources.bodyTextStyle()).setFixedLeading(ReportResources.BODY_TEXT_LEADING)));
        table.addCell(TableUtil.getLayoutCell(1, 2).setHeight(TABLE_SPACER_HEIGHT)); // Spacer

        int therapyGeneCount = EvidenceItems.uniqueEventCount(tumorSpecificEvidence);
        int therapyCount = EvidenceItems.uniqueTherapyCount(tumorSpecificEvidence);
        table.addCell(createMiddleAlignedCell().add(new Paragraph("Gene alteration(s) with therapy indication(s)")
                .addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createTreatmentIndicationCell(therapyGeneCount, therapyCount, "treatments"));

        int studyGeneCount = ClinicalTrials.uniqueOnLabelEventCount(trials);
        int studyCount = ClinicalTrials.uniqueOnLabelStudies(trials);
        table.addCell(createMiddleAlignedCell().add(new Paragraph("Gene alteration(s) with clinical study eligibility").addStyle(
                ReportResources.bodyTextStyle())));
        table.addCell(createTreatmentIndicationCell(studyGeneCount, studyCount, "studies"));

        div.add(table);

        reportDocument.add(div);
    }

    private void renderTumorCharacteristics(@NotNull final AnalysedPatientReport patientReport, @NotNull final Document reportDocument) {
        final boolean hasReliablePurityFit = patientReport.hasReliablePurityFit();

        Div div = createSectionStartDiv(contentWidth());

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1, .33f, .66f }));
        table.setWidth(contentWidth());
        table.addCell(TableUtil.getLayoutCell()
                .add(new Paragraph("Tumor characteristics summary").addStyle(ReportResources.sectionTitleStyle())));
        table.addCell(TableUtil.getLayoutCell(1, 2)
                .add(new Paragraph("Whole genome sequencing based tumor characteristics.").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(TableUtil.getLayoutCell(1, 3).setHeight(TABLE_SPACER_HEIGHT)); // Spacer

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

        table.addCell(createMiddleAlignedCell().add(new Paragraph("Tumor mutational load").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell(2).add(createHighlightParagraph(MutationalLoad.interpretToString(patientReport.tumorMutationalLoad(),
                hasReliablePurityFit)).addStyle(dataStyle)));

        table.addCell(createMiddleAlignedCell().add(new Paragraph("Microsatellite (in)stability").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell(2).add(createHighlightParagraph(MicroSatelliteStatus.interpretToString(patientReport.microsatelliteIndelsPerMb()))
                .addStyle(dataStyle)));
        div.add(table);
        reportDocument.add(div);
    }

    private static void renderTumorPurity(boolean hasReliablePurityFit, @NotNull final String valueLabel, double value, double min,
            double max, final Table table) {
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
        table.addCell(TableUtil.getLayoutCell()
                .add(new Paragraph("Genomic alterations \nsummary").addStyle(ReportResources.sectionTitleStyle())));
        table.addCell(TableUtil.getLayoutCell()
                .add(new Paragraph("Summary on genomic alterations "
                        + "(somatic variants, copy number changes, gene disruptions and gene fusions).").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(TableUtil.getLayoutCell(1, 2).setHeight(TABLE_SPACER_HEIGHT)); // Spacer

        final String[] driverVariantGenes = SomaticVariants.variantsWithDriver(patientReport.reportableVariants());
        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Genes with driver variant").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(driverVariantGenes));

        final int reportedVariants = SomaticVariants.countSomaticVariants(patientReport.reportableVariants());
        Style reportedVariantsStyle =
                (reportedVariants > 0) ? ReportResources.dataHighlightStyle() : ReportResources.dataHighlightNaStyle();
        table.addCell(createMiddleAlignedCell().add(new Paragraph("Nr. of reported variants").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createMiddleAlignedCell().add(createHighlightParagraph(String.valueOf(reportedVariants)).addStyle(
                reportedVariantsStyle)));

        final String[] copyGainGenes = GeneCopyNumbers.amplificationGenes(patientReport.geneCopyNumbers());
        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Genes with copy-gain").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(copyGainGenes));

        final String[] copyLossGenes = GeneCopyNumbers.lossGenes(patientReport.geneCopyNumbers());
        table.addCell(createMiddleAlignedCell().setVerticalAlignment(VerticalAlignment.TOP)
                .add(new Paragraph("Genes with copy-loss").addStyle(ReportResources.bodyTextStyle())));
        table.addCell(createGeneListCell(copyLossGenes));

        final String[] fusionGenes = GeneFusions.geneFusions(patientReport.geneFusions());
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
        return TableUtil.getLayoutCell(1, colSpan).setVerticalAlignment(VerticalAlignment.MIDDLE);
    }

    @NotNull
    private static Cell createGeneListCell(@NotNull String[] genes) {
        String geneString = (genes.length > 0) ? String.join(", ", genes) : DataUtil.NONE_STRING;

        Style style = (genes.length > 0) ? ReportResources.dataHighlightStyle() : ReportResources.dataHighlightNaStyle();

        return createMiddleAlignedCell().add(createHighlightParagraph(geneString)).addStyle(style);
    }

    @NotNull
    private static Paragraph createHighlightParagraph(String text) {
        return new Paragraph(text).setFixedLeading(14);
    }

    @NotNull
    private static Cell createTreatmentIndicationCell(int geneCount, int treatmentCount, @NotNull String treatmentsName) {
        String treatmentText;
        Style style;
        if (geneCount > 0) {
            treatmentText = String.format("%d (%d %s)", geneCount, treatmentCount, treatmentsName);
            style = ReportResources.dataHighlightStyle();
        } else {
            treatmentText = DataUtil.NONE_STRING;
            style = ReportResources.dataHighlightNaStyle();
        }

        return createMiddleAlignedCell().add(createHighlightParagraph(treatmentText)).addStyle(style);
    }

    @NotNull
    private static InlineBarChart createInlineBarChart(double v, double min, double max) {
        InlineBarChart chart = new InlineBarChart(v, min, max);
        chart.setWidth(41);
        chart.setHeight(6);
        return chart;
    }
}
