package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.cfreport.DataUtility;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.DataLabel;
import com.hartwig.hmftools.patientreporter.cfreport.components.InlineBarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.LineDivider;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableHelper;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.Style;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;
import com.itextpdf.layout.property.VerticalAlignment;
import org.jetbrains.annotations.NotNull;

import java.text.DecimalFormat;
import java.util.StringJoiner;

public class SummaryChapter extends ReportChapter {

    private final Style BODY_TEXT_STYLE = ReportResources.bodyTextStyle();

    private final static float TABLE_SPACER_HEIGHT = 5;

    @Override
    public final String getName() {
        return "Summary";
    }

    @Override
    public final ChapterType getChapterType() {
        return ChapterType.SummaryChapter;
    }

    @Override
    protected final void renderChapterContent(@NotNull final AnalysedPatientReport patientReport, @NotNull final Document reportDocument) {
        renderTumorLocationAndType(patientReport, reportDocument);
        renderSummaryText(patientReport, reportDocument);
        renderTreatmentIndications(patientReport, reportDocument);
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

        table.addCell(TableHelper.getLayoutCell()
                .add(new Paragraph("PRIMARY TUMOR LOCATION")
                        .addStyle(ReportResources.subTextStyle())));
        table.addCell(TableHelper.getLayoutCell()
                .add(new Paragraph("CANCER SUBTYPE")
                        .addStyle(ReportResources.subTextStyle())));
        table.addCell(TableHelper.getLayoutCell()
                .add(DataLabel.createDataLabel(patientReport.sampleReport().primaryTumorLocationString())));
        table.addCell(TableHelper.getLayoutCell()
                .add(DataLabel.createDataLabel(patientReport.sampleReport().cancerSubTypeString())));

        div.add(table);
        reportDocument.add(div);

    }

    private void renderSummaryText(@NotNull final AnalysedPatientReport patientReport, @NotNull final Document reportDocument) {

        Div div = createSectionStartDiv(getContentWidth());
        div.add(new Paragraph("Summary")
                .addStyle(ReportResources.sectionTitleStyle()));

        // @TODO Report summary not in data
        // Add content to div
        String content = "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Vivamus eget porta turpis. Lorem ipsum dolor sit amet, " +
                "consectetur adipiscing elit. Nullam interdum sodales ullamcorper. Nulla vestibulum ipsum quis ipsum congue, quis commodo velit " +
                "condimentum. Suspendisse eget nulla egestas, fermentum urna ut, bibendum ipsum. Nulla varius, dui elementum faucibus ultricies, " +
                "nisi velit dignissim arcu, nec feugiat dui magna eu felis. Maecenas at odio pharetra, sodales velit vitae, gravida mauris. Pellentesque " +
                "id ultrices diam. Integer non ex ut neque auctor pellentesque. Ut et nibh faucibus, pretium erat efficitur, vehicula lorem.";

        div.add(new Paragraph(content)
                .setWidth(getContentWidth())
                .addStyle(BODY_TEXT_STYLE).setFixedLeading(11));

        reportDocument.add(div);

    }

    private void renderTreatmentIndications(@NotNull final AnalysedPatientReport patientReport, @NotNull Document reportDocument) {

        // Initialize div
        Div div = createSectionStartDiv(getContentWidth());

        // Initialize table
        Table table = new Table(UnitValue.createPercentArray(new float[] {1, 1}));
        table.setWidth(getContentWidth());
        table.addCell(TableHelper.getLayoutCell()
                .add(new Paragraph("Treatment indications")
                        .addStyle(ReportResources.sectionTitleStyle())));

        table.addCell(TableHelper.getLayoutCell()
                .add(new Paragraph("Summary of number of alterations with number of treatment indication and/or clinical studies")
                        .addStyle(BODY_TEXT_STYLE)
                        .setFixedLeading(ReportResources.BODY_TEXT_LEADING)));
        table.addCell(TableHelper.getLayoutCell(1, 2).setHeight(TABLE_SPACER_HEIGHT)); // Spacer

        // @TODO Number of alterations/treatments not directly in the data?
        // Alterations/therapy
        int therapyGeneCount = 1;
        int therapyCount = 8;
        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph("Gene alteration(s) with therapy indication(s)")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(createTreatmentIndicationCell(therapyGeneCount, therapyCount, "treatments"));

        // @TODO Number of alterations/studies not directly in the data?
        // Alterations/clinical study
        int studyGeneCount = 2;
        int studyCount = 7;
        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph("Gene alteration(s) with clinical study eligibility")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(createTreatmentIndicationCell(studyGeneCount, studyCount, "studies"));

        div.add(table);

        reportDocument.add(div);

    }

    private void renderTumorCharacteristics(@NotNull final AnalysedPatientReport patientReport, @NotNull final Document reportDocument) {

        // Initialize div
        Div div = createSectionStartDiv(getContentWidth());

        // Initialize table
        Table table = new Table(UnitValue.createPercentArray(new float[] {1, .33f, .66f}));
        table.setWidth(getContentWidth());
        table.addCell(TableHelper.getLayoutCell()
                .add(new Paragraph("Tumor characteristics summary")
                        .addStyle(ReportResources.sectionTitleStyle())));
        table.addCell(TableHelper.getLayoutCell(1, 2).add(
                new Paragraph("Whole genome sequencing based tumor characteristics.")
                    .addStyle(BODY_TEXT_STYLE)));
        table.addCell(TableHelper.getLayoutCell(1, 3).setHeight(TABLE_SPACER_HEIGHT)); // Spacer


        // Tumor purity
        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph("Tumor purity of biopsy")
                        .addStyle(BODY_TEXT_STYLE)));
        if (patientReport.hasReliablePurityFit()) {

            final double impliedPurity = patientReport.impliedPurity();
            final double impliedPurityPercentage = DataUtility.mapPercentage(impliedPurity, DataUtility.IMPLIED_TUMOR_PURITY_MIN, DataUtility.IMPLIED_TUMOR_PURITY_MAX);

            table.addCell(createMiddleAlignedCell()
                    .add(createHighlightParagraph(DataUtility.formatPercentage(impliedPurityPercentage))
                            .addStyle(ReportResources.dataHighlightStyle())));

            table.addCell(createMiddleAlignedCell()
                    .add(createInlineBarChart((float) impliedPurity, (float) DataUtility.IMPLIED_TUMOR_PURITY_MIN, (float) DataUtility.IMPLIED_TUMOR_PURITY_MAX)));

        } else {

            // Purity below detection threshold
            table.addCell(createMiddleAlignedCell(1, 2)
                    .add(createHighlightParagraph("N/A (Below detection threshold)")
                            .addStyle(ReportResources.dataHighlightNaStyle())));

        }


        // Tumor ploidy
        String ploidyString;
        Style ploidyStyle;
        if (patientReport.hasReliablePurityFit()) {
            ploidyString = new DecimalFormat("#.#")
                    .format(patientReport.averageTumorPloidy());
            ploidyStyle = ReportResources.dataHighlightStyle();
        } else {
            ploidyString = "N/A";
            ploidyStyle = ReportResources.dataHighlightNaStyle();
        }
        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph("Average tumor ploidy")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(createMiddleAlignedCell(1, 2)
                .add(createHighlightParagraph(ploidyString)
                        .addStyle(ploidyStyle)));


        // Tumor mutational load
        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph("Tumor mutational load")
                        .addStyle(BODY_TEXT_STYLE)));
        if (patientReport.hasReliablePurityFit()) {

            final int mutationalLoad = patientReport.tumorMutationalLoad();
            final String mutationalLoadString = DataUtility.MutationalLoad.interpretToString(mutationalLoad);

            table.addCell(createMiddleAlignedCell()
                    .add(createHighlightParagraph(mutationalLoadString)
                            .addStyle(ReportResources.dataHighlightStyle())));

            // @TODO : Check chart mapping
            table.addCell(createMiddleAlignedCell()
                    .add(createInlineBarChart((float) mutationalLoad,
                            (float) DataUtility.MutationalLoad.RANGE_MIN,
                            (float) DataUtility.MutationalLoad.RANGE_MAX)));

        } else {

            table.addCell(createMiddleAlignedCell(1, 2)
                    .add(createHighlightParagraph("N/A")
                            .addStyle(ReportResources.dataHighlightNaStyle())));

        }


        // Microsatellite stability
        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph("Microsatellite (in)stability")
                        .addStyle(BODY_TEXT_STYLE)));
        if (patientReport.hasReliablePurityFit()) {

            final double microSatelliteIndels = patientReport.microsatelliteIndelsPerMb();
            final String microSatelliteStabilityString = DataUtility.MicroSatellite.interpretToString(microSatelliteIndels);

            table.addCell(createMiddleAlignedCell()
                    .add(createHighlightParagraph(microSatelliteStabilityString)
                            .addStyle(ReportResources.dataHighlightStyle())));

            // @TODO : Check chart mapping
            table.addCell(createMiddleAlignedCell()
                    .add(createInlineBarChart((float) microSatelliteIndels, 0f, 10f)));

        } else {

            table.addCell(createMiddleAlignedCell(1, 2)
                    .add(createHighlightParagraph("N/A")
                            .addStyle(ReportResources.dataHighlightNaStyle())));

        }

        div.add(table);

        reportDocument.add(div);

    }

    private void renderGenomicAlterations(@NotNull final AnalysedPatientReport patientReport, @NotNull Document report) {

        // Initialize div
        Div div = createSectionStartDiv(getContentWidth());

        // Initialize table
        Table table = new Table(UnitValue.createPercentArray(new float[] {1, 1}));
        table.setWidth(getContentWidth());
        table.addCell(TableHelper.getLayoutCell()
                .add(new Paragraph("Genomic alterations \nsummary")
                        .addStyle(ReportResources.sectionTitleStyle())));
        table.addCell(TableHelper.getLayoutCell()
                .add(new Paragraph("Summary on genomic alterations " +
                "(somatic variants, copy number changes, gene disruptions and gene fusions).")
                .addStyle(BODY_TEXT_STYLE)));
        table.addCell(TableHelper.getLayoutCell(1, 2).setHeight(TABLE_SPACER_HEIGHT)); // Spacer

        // Genes with driver variant
        String[] driverVariantGenes = {"CDKN2A", "BRAF"};
        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph("Genes with driver variant")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(createGeneListCell(driverVariantGenes));

        // Reported variants
        int reportedVariants = 4;
        Style reportedVariantsStyle = (reportedVariants > 0) ? ReportResources.dataHighlightStyle() : ReportResources.dataHighlightNaStyle();
        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph("Nr. of reported variants")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(createMiddleAlignedCell()
                .add(createHighlightParagraph(String.valueOf(reportedVariants))
                .addStyle(reportedVariantsStyle)));

        // Copy gain genes
        String[] copyGainGenes = {};
        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph("Genes with copy-gain")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(createGeneListCell(copyGainGenes));

        // Copy loss genes
        String[] copyLossGenes = {"PTEN"};
        table.addCell(createMiddleAlignedCell()
                .add(new Paragraph("Genes with copy-loss")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(createGeneListCell(copyLossGenes));

        // Gene fusions
        String[] fusionGenes = {};
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
        return TableHelper.getLayoutCell(rowspan, colspan)
                .setVerticalAlignment(VerticalAlignment.MIDDLE);
    }

    @NotNull
    private static Cell createGeneListCell(@NotNull String[] genes) {

        // Concatenate genes
        String geneString;
        if (genes.length == 0) {
            geneString = "NONE";
        } else {

            StringJoiner joiner = new StringJoiner(", ");
            for (String s: genes) {
                joiner.add(s);
            }
            geneString = joiner.toString();

        }

        // Fetch style
        Style style = genes.length > 0 ? ReportResources.dataHighlightStyle() : ReportResources.dataHighlightNaStyle();

        // Build table
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
            treatmentText = "NONE";
            style = ReportResources.dataHighlightNaStyle();
        }

        return createMiddleAlignedCell()
                .add(createHighlightParagraph(treatmentText))
                .addStyle(style);

    }

    @NotNull
    private static InlineBarChart createInlineBarChart(float v, float min, float max) {
        InlineBarChart chart = new InlineBarChart(v, min, max);
        chart.setWidth(41);
        chart.setHeight(6);
        return chart;
    }

}
