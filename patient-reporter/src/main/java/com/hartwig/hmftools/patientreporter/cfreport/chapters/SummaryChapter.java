package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.*;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.Style;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;
import com.itextpdf.layout.property.VerticalAlignment;
import org.jetbrains.annotations.NotNull;

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
    protected final void renderChapterContent(@NotNull Document report) {
        renderTumorLocationAndType(report);
        renderSummaryText(report);
        renderTreatmentIndications(report);
        renderTumorCharacteristics(report);
        renderGenomicAlterations(report);
    }

    private final void renderTumorLocationAndType(@NotNull Document report) {

        Div div = new Div();
        div.setKeepTogether(true);
        div.setWidth(getContentWidth());

        // Initialize table
        Table table = new Table(UnitValue.createPercentArray(new float[] {1, 1}));
        table.setWidth(getContentWidth());

        table.addCell(TableHelper.getLayoutCell()
                .add(new Paragraph("PRIMARY TUMOR LOCATION")
                        .addStyle(ReportResources.subTextStyle())));
        table.addCell(TableHelper.getLayoutCell()
                .add(new Paragraph("CANCER SUBTYPE")
                        .addStyle(ReportResources.subTextStyle())));
        table.addCell(TableHelper.getLayoutCell()
                .add(DataLabel.createDataLabel("Skin")));
        table.addCell(TableHelper.getLayoutCell()
                .add(DataLabel.createDataLabel("Melanoma")));

        div.add(table);
        report.add(div);

    }

    private final void renderSummaryText(@NotNull Document report) {

        Div div = initializeSummarySectionDiv(getContentWidth());
        div.add(new Paragraph("Summary")
                .addStyle(ReportResources.sectionTitleStyle()));

        // Add content to div
        String content = "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Vivamus eget porta turpis. Lorem ipsum dolor sit amet, " +
                "consectetur adipiscing elit. Nullam interdum sodales ullamcorper. Nulla vestibulum ipsum quis ipsum congue, quis commodo velit " +
                "condimentum. Suspendisse eget nulla egestas, fermentum urna ut, bibendum ipsum. Nulla varius, dui elementum faucibus ultricies, " +
                "nisi velit dignissim arcu, nec feugiat dui magna eu felis. Maecenas at odio pharetra, sodales velit vitae, gravida mauris. Pellentesque " +
                "id ultrices diam. Integer non ex ut neque auctor pellentesque. Ut et nibh faucibus, pretium erat efficitur, vehicula lorem.";

        div.add(new Paragraph(content)
                .setWidth(getContentWidth())
                .addStyle(BODY_TEXT_STYLE));

        report.add(div);

    }

    private final void renderTreatmentIndications(@NotNull Document report) {

        // Initialize div
        Div div = initializeSummarySectionDiv(getContentWidth());

        // Initialize table
        Table table = new Table(UnitValue.createPercentArray(new float[] {1, 1}));
        table.setWidth(getContentWidth());
        table.addCell(TableHelper.getLayoutCell()
                .add(new Paragraph("Treatment indications")
                        .addStyle(ReportResources.sectionTitleStyle())));

        table.addCell(TableHelper.getLayoutCell()
                .add(new Paragraph("Summary of number of alterations with number of treatment indication and/or clinical studies")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(TableHelper.getLayoutCell(1, 2).setHeight(TABLE_SPACER_HEIGHT)); // Spacer

        // Alterations/therapy
        int therapyGeneCount = 1;
        int therapyCount = 8;
        table.addCell(getMiddleAlignedCell()
                .add(new Paragraph("Gene alteration(s) with therapy indication(s)").addStyle(BODY_TEXT_STYLE)));
        table.addCell(getTreatmentIndicationCell(therapyGeneCount, therapyCount, "treatments"));

        // Alterations/clinical study
        int studyGeneCount = 2;
        int studyCount = 7;
        table.addCell(getMiddleAlignedCell()
                .add(new Paragraph("Gene alteration(s) with clinical study eligibility").addStyle(BODY_TEXT_STYLE)));
        table.addCell(getTreatmentIndicationCell(studyGeneCount, studyCount, "studies"));

        div.add(table);

        report.add(div);

    }

    private final void renderTumorCharacteristics(@NotNull Document report) {

        // Initialize div
        Div div = initializeSummarySectionDiv(getContentWidth());

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
        float tumorPurity = 100; //74.4f;
        table.addCell(getMiddleAlignedCell()
                .add(new Paragraph("Tumor purity of biopsy")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(getMiddleAlignedCell()
                .add(new Paragraph(String.format(java.util.Locale.US,"%.0f%%", tumorPurity))
                        .addStyle(ReportResources.dataHighlightStyle())));
        table.addCell(getMiddleAlignedCell()
                .add(getInlineBarChart(tumorPurity, 0f, 100f)));

        // Tumor characteristics
        float ploidy = 3.1f;
        table.addCell(getMiddleAlignedCell()
                .add(new Paragraph("Average tumor ploidy")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(getMiddleAlignedCell()
                .add(new Paragraph(String.format(java.util.Locale.US, "%.1f", ploidy))
                        .addStyle(ReportResources.dataHighlightStyle())));
        table.addCell(getMiddleAlignedCell());

        // Tumor mutational load
        String mutationalLoad = "High";
        table.addCell(getMiddleAlignedCell()
                .add(new Paragraph("Tumor mutational load")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(getMiddleAlignedCell()
                .add(new Paragraph(mutationalLoad))
                .addStyle(ReportResources.dataHighlightStyle()));
        table.addCell(getMiddleAlignedCell()
                .add(getInlineBarChart(90, 0f, 100f)));

        // Microsatellite stability
        String microsatelliteStability = "Stable";
        table.addCell(getMiddleAlignedCell()
                .add(new Paragraph("Microsatellite (in)stability")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(getMiddleAlignedCell()
                .add(new Paragraph(microsatelliteStability)
                        .addStyle(ReportResources.dataHighlightStyle())));
        table.addCell(getMiddleAlignedCell()
                .add(getInlineBarChart(.6f, 0f, 10f)));

        div.add(table);

        report.add(div);

    }

    private final void renderGenomicAlterations(@NotNull Document report) {

        // Initialize div
        Div div = initializeSummarySectionDiv(getContentWidth());

        // Initialize table
        Table table = new Table(UnitValue.createPercentArray(new float[] {1, 1}));
        table.setWidth(getContentWidth());
        table.addCell(TableHelper.getLayoutCell()
                .add(new Paragraph("Genomic alterations r\nsummary")
                        .addStyle(ReportResources.sectionTitleStyle())));
        table.addCell(TableHelper.getLayoutCell()
                .add(new Paragraph("Summary on genomic alterations " +
                "(somatic variants, copy number changes, gene disruptions and gene fusions).")
                .addStyle(BODY_TEXT_STYLE)));
        table.addCell(TableHelper.getLayoutCell(1, 2).setHeight(TABLE_SPACER_HEIGHT)); // Spacer

        // Genes with driver variant
        String[] driverVariantGenes = {"CDKN2A", "BRAF"};
        table.addCell(getMiddleAlignedCell()
                .add(new Paragraph("Genes with driver variant")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(getGeneListCell(driverVariantGenes));

        // Reported variants
        int reportedVariants = 4;
        Style reportedVariantsStyle = (reportedVariants > 0) ? ReportResources.dataHighlightStyle() : ReportResources.dataHighlightNaStyle();
        table.addCell(getMiddleAlignedCell()
                .add(new Paragraph("Nr. of reported variants")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(getMiddleAlignedCell()
                .add(new Paragraph(String.valueOf(reportedVariants))
                .addStyle(reportedVariantsStyle)));

        // Copy gain genes
        String[] copyGainGenes = {};
        table.addCell(getMiddleAlignedCell()
                .add(new Paragraph("Genes with copy-gain")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(getGeneListCell(copyGainGenes));

        // Copy loss genes
        String[] copyLossGenes = {"PTEN"};
        table.addCell(getMiddleAlignedCell()
                .add(new Paragraph("Genes with copy-loss")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(getGeneListCell(copyLossGenes));

        // Gene fusions
        String[] fusionGenes = {};
        table.addCell(getMiddleAlignedCell()
                .add(new Paragraph("Gene fusions")
                        .addStyle(BODY_TEXT_STYLE)));
        table.addCell(getGeneListCell(fusionGenes));

        div.add(table);

        report.add(div);


    }

    @NotNull
    private static final Div initializeSummarySectionDiv(float width) {

        Div div = new Div();
        div.setKeepTogether(true);
        div.setWidth(width);

        // Add divider and section title
        div.add(LineDivider.createLineDivider(width));

        return div;

    }

    @NotNull
    private static final Cell getMiddleAlignedCell() {
        Cell c = TableHelper.getLayoutCell()
                .setVerticalAlignment(VerticalAlignment.MIDDLE);
        return c;
    }

    @NotNull
    private static final Cell getGeneListCell(@NotNull String[] genes) {

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
        Cell c = getMiddleAlignedCell()
                .add(new Paragraph(geneString))
                .addStyle(style);
                return c;

    }

    @NotNull
    private static final Cell getTreatmentIndicationCell(int geneCount, int treatmentCount, @NotNull String treatmentsName) {

        String treatmentText;
        Style style;
        if (geneCount > 0) {
            treatmentText = String.format("%d (%d %s)", geneCount, treatmentCount, treatmentsName);
            style = ReportResources.dataHighlightStyle();
        } else {
            treatmentText = "NONE";
            style = ReportResources.dataHighlightNaStyle();
        }

        return getMiddleAlignedCell()
                .add(new Paragraph(treatmentText))
                .addStyle(style);

    }

    @NotNull
    private static final InlineBarChart getInlineBarChart(float v, float min, float max) {
        InlineBarChart tumorPurityChart = new InlineBarChart(v, min, max);
        tumorPurityChart.setWidth(41);
        tumorPurityChart.setHeight(6);
        return tumorPurityChart;
    }

}
