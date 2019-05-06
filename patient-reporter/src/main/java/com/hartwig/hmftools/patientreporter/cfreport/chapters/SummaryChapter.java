package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.common.actionability.ClinicalTrial;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.cfreport.MathUtil;
import com.hartwig.hmftools.patientreporter.cfreport.data.*;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.DataLabel;
import com.hartwig.hmftools.patientreporter.cfreport.components.InlineBarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.LineDivider;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableUtil;
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

    private final Style BODY_TEXT_STYLE = ReportResources.bodyTextStyle();

    private final static float TABLE_SPACER_HEIGHT = 5;

    private final AnalysedPatientReport patientReport;

    public SummaryChapter(@NotNull final AnalysedPatientReport patientReport) {
        this.patientReport = patientReport;
    }

    @Override
    public final String getName() {
        return "Summary";
    }

    @Override
    public String getPageNumberPrefix() {
        return getName();
    }

    public boolean isFullWidth() {
        return false;
    }

    public boolean hasCompleteSidebar() {
        return true;
    }

    @Override
    public final void render(@NotNull final Document reportDocument) {
        renderTumorLocationAndType(patientReport, reportDocument);

        // @TODO Replace this fixed text with the patientReport.summaryText method.
        // Return value from that method can be null which is gracefully handled by renderSummaryText :
        // final String summaryContent = patientReport.summaryText();
        final String summaryContent = "Melanoma sample with an activating BRAF mutation that is associated with " +
                "response to BRAF-inhibitors (in combination with an MEK-inhibitor). The tumor shows a complete " +
                "inactivation of CDKN2A, indicating potential benefit of CDK4/6 inhibitors (e.g. palbociclib). The " +
                "observed complete loss of PTEN likely results in an activation of the PI3K-AKT-mTOR pathway and " +
                "suggests eligibility for treatment (study) using mTOR/PI3K inhibitors. In addition, the tumor samples " +
                "shows a high mutational burden that is associated with an increased response rate to checkpoint " +
                "inhibitor immunotherapy.";
        renderSummaryText(summaryContent, reportDocument);

        renderTreatmentIndications(patientReport.tumorSpecificEvidence(), patientReport.clinicalTrials(), reportDocument);
        renderTumorCharacteristics(patientReport, reportDocument);
        renderGenomicAlterations(patientReport, reportDocument);
    }

    private void renderTumorLocationAndType(@NotNull final AnalysedPatientReport patientReport, @NotNull final Document reportDocument) {

        Div div = new Div();
        div.setKeepTogether(true);
        div.setWidth(getContentWidth());

        // Initialize table
        Table table = new Table(UnitValue.createPercentArray(new float[] {1, 1}));
        table.setWidth(getContentWidth());

        patientReport.sampleReport().patientTumorLocation();

        table.addCell(TableUtil.getLayoutCell()
                .add(new Paragraph("PRIMARY TUMOR LOCATION")
                        .addStyle(ReportResources.subTextStyle())));
        table.addCell(TableUtil.getLayoutCell()
                .add(new Paragraph("CANCER SUBTYPE")
                        .addStyle(ReportResources.subTextStyle())));
        table.addCell(TableUtil.getLayoutCell()
                .add(DataLabel.createDataLabel(patientReport.sampleReport().primaryTumorLocationString())));
        table.addCell(TableUtil.getLayoutCell()
                .add(DataLabel.createDataLabel(patientReport.sampleReport().cancerSubTypeString())));

        div.add(table);
        reportDocument.add(div);

    }

    private void renderSummaryText(@Nullable final String text, @NotNull final Document reportDocument) {

        if (text == null) {
            return;
        }

        Div div = createSectionStartDiv(getContentWidth());
        div.add(new Paragraph("Summary")
                .addStyle(ReportResources.sectionTitleStyle()));

        div.add(new Paragraph(text)
                .setWidth(getContentWidth())
                .addStyle(BODY_TEXT_STYLE).setFixedLeading(11));

        reportDocument.add(div);

    }

    private void renderTreatmentIndications(@NotNull final List<EvidenceItem> tumorSpecificEvidence,
                                            @NotNull final List<ClinicalTrial> trials, @NotNull Document reportDocument) {

        // Initialize div
        Div div = createSectionStartDiv(getContentWidth());

        // Initialize table
        Table table = new Table(UnitValue.createPercentArray(new float[] {1, 1}));
        table.setWidth(getContentWidth());
        table.addCell(TableUtil.getLayoutCell()
                .add(new Paragraph("Treatment indications")
                        .addStyle(ReportResources.sectionTitleStyle())));

        table.addCell(TableUtil.getLayoutCell()
                .add(new Paragraph("Summary of number of alterations with number of treatment indication and/or clinical studies")
                        .addStyle(BODY_TEXT_STYLE)
                        .setFixedLeading(ReportResources.BODY_TEXT_LEADING)));
        table.addCell(TableUtil.getLayoutCell(1, 2).setHeight(TABLE_SPACER_HEIGHT)); // Spacer

        // Alterations/therapies
        int therapyGeneCount = EvidenceItems.uniqueEventCount(tumorSpecificEvidence);
        int therapyCount = EvidenceItems.uniqueTherapyCount(tumorSpecificEvidence);
        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph("Gene alteration(s) with therapy indication(s)")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(createTreatmentIndicationCell(therapyGeneCount, therapyCount, "treatments"));

        // Alterations/clinical study
        int studyGeneCount = ClinicalTrials.uniqueOnLabelEventCount(trials);
        int studyCount = ClinicalTrials.uniqueOnLabelStudies(trials);
        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph("Gene alteration(s) with clinical study eligibility")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(createTreatmentIndicationCell(studyGeneCount, studyCount, "studies"));

        div.add(table);

        reportDocument.add(div);

    }

    private void renderTumorCharacteristics(@NotNull final AnalysedPatientReport patientReport, @NotNull final Document reportDocument) {

        final boolean hasReliablePurityFit = patientReport.hasReliablePurityFit();

        // Initialize div
        Div div = createSectionStartDiv(getContentWidth());

        // Initialize table
        Table table = new Table(UnitValue.createPercentArray(new float[] {1, .33f, .66f}));
        table.setWidth(getContentWidth());
        table.addCell(TableUtil.getLayoutCell()
                .add(new Paragraph("Tumor characteristics summary")
                        .addStyle(ReportResources.sectionTitleStyle())));
        table.addCell(TableUtil.getLayoutCell(1, 2).add(
                new Paragraph("Whole genome sequencing based tumor characteristics.")
                    .addStyle(BODY_TEXT_STYLE)));
        table.addCell(TableUtil.getLayoutCell(1, 3).setHeight(TABLE_SPACER_HEIGHT)); // Spacer

        // Tumor purity
        final double impliedPurity = patientReport.impliedPurity();
        final double impliedPurityPercentage = MathUtil.mapPercentage(impliedPurity, TumorPurity.RANGE_MIN, TumorPurity.RANGE_MAX);
        renderTumorCharacteristicBarChartRow(
                hasReliablePurityFit,
                "Tumor purity of biopsy",
                DataUtil.formatPercentage(impliedPurityPercentage),
                impliedPurity,
                TumorPurity.RANGE_MIN,
                TumorPurity.RANGE_MAX,
                table
        );

        Style dataStyle = hasReliablePurityFit ? ReportResources.dataHighlightStyle()
                : ReportResources.dataHighlightNaStyle();

        // Tumor ploidy
        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph("Average tumor ploidy")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(createMiddleAlignedCell(1, 2)
                .add(createHighlightParagraph(GeneUtil.getPloidyToCopiesString(patientReport.averageTumorPloidy(), hasReliablePurityFit))
                        .addStyle(dataStyle)));

        // Tumor mutational load
        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph( "Tumor mutational load")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(createMiddleAlignedCell(1, 2)
                .add(createHighlightParagraph(MutationalLoad.interpretToString(patientReport.tumorMutationalLoad(),
                        hasReliablePurityFit))
                        .addStyle(dataStyle)));

        // Microsatellite stability
         table.addCell(createMiddleAlignedCell()
                .add(new Paragraph( "Microsatellite (in)stability")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(createMiddleAlignedCell(1, 2)
                .add(createHighlightParagraph(MicroSatelliteStatus.interpretToString(patientReport.microsatelliteIndelsPerMb(),
                        hasReliablePurityFit))
                        .addStyle(dataStyle)));

        div.add(table);
        reportDocument.add(div);

    }

    private void renderTumorCharacteristicBarChartRow(boolean hasReliablePurityFit, @NotNull String label, @NotNull final String valueLabel, double value, double min, double max, final Table table) {

        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph(label)
                        .addStyle(BODY_TEXT_STYLE)));

        if (hasReliablePurityFit) {

            table.addCell(createMiddleAlignedCell()
                    .add(createHighlightParagraph(valueLabel)
                            .addStyle(ReportResources.dataHighlightStyle())));

            table.addCell(createMiddleAlignedCell()
                    .add(createInlineBarChart(value, min, max)));

        } else {

            table.addCell(createMiddleAlignedCell(1, 2)
                    .add(createHighlightParagraph(DataUtil.NAString)
                            .addStyle(ReportResources.dataHighlightNaStyle())));

        }

    }

    private void renderGenomicAlterations(@NotNull final AnalysedPatientReport patientReport, @NotNull Document report) {

        // Initialize div
        final Div div = createSectionStartDiv(getContentWidth());

        // Initialize table
        final Table table = new Table(UnitValue.createPercentArray(new float[] {1, 1}));
        table.setWidth(getContentWidth());
        table.addCell(TableUtil.getLayoutCell()
                .add(new Paragraph("Genomic alterations \nsummary")
                        .addStyle(ReportResources.sectionTitleStyle())));
        table.addCell(TableUtil.getLayoutCell()
                .add(new Paragraph("Summary on genomic alterations " +
                "(somatic variants, copy number changes, gene disruptions and gene fusions).")
                .addStyle(BODY_TEXT_STYLE)));
        table.addCell(TableUtil.getLayoutCell(1, 2).setHeight(TABLE_SPACER_HEIGHT)); // Spacer

        // Genes with driver variant
        final String[] driverVariantGenes = SomaticVariants.somaticVariantsWithDriver(patientReport.somaticVariants());
        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph("Genes with driver variant")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(createGeneListCell(driverVariantGenes));

        // Reported variants
        final int reportedVariants = SomaticVariants.countSomaticVariants(patientReport.somaticVariants());
        Style reportedVariantsStyle = (reportedVariants > 0) ? ReportResources.dataHighlightStyle() : ReportResources.dataHighlightNaStyle();
        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph("Nr. of reported variants")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(createMiddleAlignedCell()
                .add(createHighlightParagraph(String.valueOf(reportedVariants))
                .addStyle(reportedVariantsStyle)));

        // Copy gain genes
        final String[] copyGainGenes = GeneCopyNumbers.amplificationGenes(patientReport.geneCopyNumbers());
        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph("Genes with copy-gain")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(createGeneListCell(copyGainGenes));

        // Copy loss genes
        final String[] copyLossGenes = GeneCopyNumbers.lossGenes(patientReport.geneCopyNumbers());
        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph("Genes with copy-loss")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(createGeneListCell(copyLossGenes));

        // Gene fusions
        final String[] fusionGenes = GeneFusions.geneFusions(patientReport.geneFusions());
        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph("Gene fusions")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(createGeneListCell(fusionGenes));

        div.add(table);

        report.add(div);


    }

    @NotNull
    private static Div createSectionStartDiv(float width) {

        Div div = new Div();
        div.setKeepTogether(true);
        div.setWidth(width);

        // Add divider and section title
        div.add(LineDivider
                .createLineDivider(width)
                .setMarginBottom(4));

        return div;

    }

    @NotNull
    private static Cell createMiddleAlignedCell() {
        return createMiddleAlignedCell(1, 1);
    }

    @NotNull
    private static Cell createMiddleAlignedCell(int rowspan, int colspan) {
        return TableUtil.getLayoutCell(rowspan, colspan)
                .setVerticalAlignment(VerticalAlignment.MIDDLE);
    }

    @NotNull
    private static Cell createGeneListCell(@NotNull String[] genes) {

        String geneString = (genes.length > 0)
                ? String.join(", ", genes)
                : DataUtil.NoneString;

        Style style = (genes.length > 0)
                ? ReportResources.dataHighlightStyle()
                : ReportResources.dataHighlightNaStyle();

        return createMiddleAlignedCell()
                .add(createHighlightParagraph(geneString))
                .addStyle(style);

    }

    @NotNull
    private static Paragraph createHighlightParagraph(String text) {
        return new Paragraph(text)
                .setFixedLeading(14);
    }

    @NotNull
    private static Cell createTreatmentIndicationCell(int geneCount, int treatmentCount, @NotNull String treatmentsName) {

        String treatmentText;
        Style style;
        if (geneCount > 0) {
            treatmentText = String.format("%d (%d %s)", geneCount, treatmentCount, treatmentsName);
            style = ReportResources.dataHighlightStyle();
        } else {
            treatmentText = DataUtil.NoneString;
            style = ReportResources.dataHighlightNaStyle();
        }

        return createMiddleAlignedCell()
                .add(createHighlightParagraph(treatmentText))
                .addStyle(style);

    }

    @NotNull
    private static InlineBarChart createInlineBarChart(double v, double min, double max) {
        InlineBarChart chart = new InlineBarChart(v, min, max);
        chart.setWidth(41);
        chart.setHeight(6);
        return chart;
    }

}
