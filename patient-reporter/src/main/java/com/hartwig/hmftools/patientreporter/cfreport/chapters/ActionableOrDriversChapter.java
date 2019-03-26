package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.cfreport.components.InlineBarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableHelper;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.Style;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.TextAlignment;
import org.jetbrains.annotations.NotNull;

public class ActionableOrDriversChapter extends ReportChapter {

    private final static float TABLE_BOTTOM_SPACER_HEIGHT = 30;

    @Override
    public String getName() {
        return "Actionable or drivers";
    }

    @Override
    public ChapterType getChapterType() {
        return ChapterType.ContentChapter;
    }

    @Override
    protected void renderChapterContent(Document report) {
        renderTumorVariants(report);
        renderGainsAndLosses(report);
        renderSomaticFusions(report);
        renderDisruptions(report);
    }

    private final void renderTumorVariants(@NotNull Document report) {

        final boolean isAvailable = true;

        // Section title
        report.add(getSectionTitle("Tumor type specific variants", isAvailable));

        // Section content
        if (isAvailable) {

            Table table = TableHelper.createReportContentTable(new float[]{1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
            table.addHeaderCell(TableHelper.getHeaderCell("Gene"));
            table.addHeaderCell(TableHelper.getHeaderCell("Variant"));
            table.addHeaderCell(TableHelper.getHeaderCell("Impact"));
            table.addHeaderCell(TableHelper.getHeaderCell("Read depth").setTextAlignment(TextAlignment.RIGHT));
            table.addHeaderCell(TableHelper.getHeaderCell("Hotspot"));
            table.addHeaderCell(TableHelper.getHeaderCell("Ploidy (VAF)"));
            table.addHeaderCell(TableHelper.getHeaderCell()); // Spacer for graph
            table.addHeaderCell(TableHelper.getHeaderCell("Clonality"));
            table.addHeaderCell(TableHelper.getHeaderCell("Biallelic"));
            table.addHeaderCell(TableHelper.getHeaderCell("Driver"));

            for (int i = 0; i < 4; i++) {
                table.addCell(TableHelper.getContentCell("BRAF*"));
                table.addCell(TableHelper.getContentCell("c.1799T>A"));
                table.addCell(TableHelper.getContentCell("p.Val6000Glu"));
                table.addCell(TableHelper.getContentCell("107 / 161").setTextAlignment(TextAlignment.RIGHT));
                table.addCell(TableHelper.getContentCell("Yes"));
                table.addCell(TableHelper.getContentCell("AAAABB (65%)"));
                table.addCell(TableHelper.getContentCell(new InlineBarChart(.6f, 0.0f, 1.0f)));
                table.addCell(TableHelper.getContentCell("Clonal"));
                table.addCell(TableHelper.getContentCell("-"));
                table.addCell(TableHelper.getContentCell("High"));
            }

            table.addCell(TableHelper.getLayoutCell(10, 1).setHeight(TABLE_BOTTOM_SPACER_HEIGHT));
            report.add(table);

        } else {
            report.add(getNoneTable());
        }
    }

    private final void renderGainsAndLosses(@NotNull Document report) {

        final boolean isAvailable = true;

        // Section title
        report.add(getSectionTitle("Tumor type specific gains and losses", isAvailable));

        // Section content
        if (isAvailable) {

            Table table = TableHelper.createReportContentTable(new float[]{4, 4, 4, 2, 1});
            table.addHeaderCell(TableHelper.getHeaderCell("Chromosome"));
            table.addHeaderCell(TableHelper.getHeaderCell("Chromosome band"));
            table.addHeaderCell(TableHelper.getHeaderCell("Gene"));
            table.addHeaderCell(TableHelper.getHeaderCell("Type"));
            table.addHeaderCell(TableHelper.getHeaderCell("Copies").setTextAlignment(TextAlignment.RIGHT));

            for (int i = 0; i < 4; i++) {
                table.addCell(TableHelper.getContentCell("10"));
                table.addCell(TableHelper.getContentCell("q23.31"));
                table.addCell(TableHelper.getContentCell("PTEN"));
                table.addCell(TableHelper.getContentCell("partial loss"));
                table.addCell(TableHelper.getContentCell("0").setTextAlignment(TextAlignment.RIGHT));
            }

            table.addCell(TableHelper.getLayoutCell(10, 1).setHeight(TABLE_BOTTOM_SPACER_HEIGHT));
            report.add(table);

        } else {
            report.add(getNoneTable());
        }

    }

    private final void renderSomaticFusions(@NotNull Document report) {
        final boolean isAvailable = true;

        // Section title
        report.add(getSectionTitle("Somatic gene fusions", isAvailable));

        // Section content
        if (isAvailable) {

            Table table = TableHelper.createReportContentTable(new float[]{2, .8f, 1.5f, 1.5f, 3f});
            table.addHeaderCell(TableHelper.getHeaderCell("Fusion"));
            table.addHeaderCell(TableHelper.getHeaderCell("Transcript"));
            table.addHeaderCell(TableHelper.getHeaderCell("Context"));
            table.addHeaderCell(TableHelper.getHeaderCell("Copies").setTextAlignment(TextAlignment.RIGHT));
            table.addHeaderCell(TableHelper.getHeaderCell("Source"));

            for (int i = 0; i < 4; i++) {
                table.addCell(TableHelper.getContentCell("---"));
                table.addCell(TableHelper.getContentCell("## - ##"));
                table.addCell(TableHelper.getContentCell("## - ##"));
                     table.addCell(TableHelper.getContentCell("0").setTextAlignment(TextAlignment.RIGHT));
                table.addCell(TableHelper.getContentCell("http://oncokb.org/"));
            }

            table.addCell(TableHelper.getLayoutCell(10, 1).setHeight(TABLE_BOTTOM_SPACER_HEIGHT));
            report.add(table);


        } else {
            report.add(getNoneTable());
        }

    }

    private final void renderDisruptions(@NotNull Document report) {

        final boolean isAvailable = true;

        // Section title
        report.add(getSectionTitle("Somatic gene disruptions", isAvailable));

        // Section content
        if (isAvailable) {

            Table table = TableHelper.createReportContentTable(new float[]{1, 1, 2, .8f, 1.5f, 1.5f, 1.5f});
            table.addHeaderCell(TableHelper.getHeaderCell("Location"));
            table.addHeaderCell(TableHelper.getHeaderCell("Gene"));
            table.addHeaderCell(TableHelper.getHeaderCell("Disrupted range"));
            table.addHeaderCell(TableHelper.getHeaderCell("Type"));
            table.addHeaderCell(TableHelper.getHeaderCell("Copies").setTextAlignment(TextAlignment.RIGHT));
            table.addHeaderCell(TableHelper.getHeaderCell("Gene min copies")).setTextAlignment(TextAlignment.RIGHT);
            table.addHeaderCell(TableHelper.getHeaderCell("Gene max copies").setTextAlignment(TextAlignment.RIGHT));

            for (int i = 0; i < 4; i++) {
                table.addCell(TableHelper.getContentCell("q23.31"));
                table.addCell(TableHelper.getContentCell("PTEN"));
                table.addCell(TableHelper.getContentCell("Intron 5 -> Intron 6"));
                table.addCell(TableHelper.getContentCell("DEL"));
                table.addCell(TableHelper.getContentCell("1.8").setTextAlignment(TextAlignment.RIGHT));
                table.addCell(TableHelper.getContentCell("0").setTextAlignment(TextAlignment.RIGHT));
                table.addCell(TableHelper.getContentCell("2").setTextAlignment(TextAlignment.RIGHT));
            }

            table.addCell(TableHelper.getLayoutCell(10, 1).setHeight(TABLE_BOTTOM_SPACER_HEIGHT));
            report.add(table);

        } else {
            report.add(getNoneTable());
        }

    }

    @NotNull
    private static final Paragraph getSectionTitle(@NotNull String title, boolean isAvailable) {

        Style style = ReportResources.sectionTitleStyle();
        if (!isAvailable) {
            style.setFontColor(ReportResources.PALETTE_LIGHT_GREY);
        }

        return new Paragraph(title)
                .addStyle(style);

    }

    @NotNull
    private static final Table getNoneTable() {

        Table table = TableHelper.createReportContentTable(new float[] {1});
        table.addCell(TableHelper.getDisabledContentCell(new Paragraph("NONE")));
        table.addCell(TableHelper.getLayoutCell(10, 1).setHeight(TABLE_BOTTOM_SPACER_HEIGHT));

        return table;

    }

}
