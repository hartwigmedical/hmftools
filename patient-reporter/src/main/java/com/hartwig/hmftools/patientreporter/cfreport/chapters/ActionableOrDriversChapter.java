package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.cfreport.components.InlineBarChart;
import com.hartwig.hmftools.patientreporter.cfreport.components.TableHelper;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Cell;
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
        report.add(createTumorVariantsTable());
        report.add(createGainsAndLossesTable());
        report.add(createSomaticFusionsTable());
        report.add(createDisruptionsTable());
    }

    @NotNull
    private final Table createTumorVariantsTable() {

        final String chapterTitle = "Tumor specific variants";
        final boolean isAvailable = true;

        if (!isAvailable) {
            return createNoneTable(chapterTitle);
        }

        // Create content table
        Table contentTable = TableHelper.createReportContentTable(new float[] {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, new Cell[]  {
                TableHelper.getHeaderCell("Gene"),
                TableHelper.getHeaderCell("Variant"),
                TableHelper.getHeaderCell("Impact"),
                TableHelper.getHeaderCell("Read depth").setTextAlignment(TextAlignment.RIGHT),
                TableHelper.getHeaderCell("Hotspot"),
                TableHelper.getHeaderCell("Ploidy (VAF)"),
                TableHelper.getHeaderCell(), // Spacer for graph
                TableHelper.getHeaderCell("Clonality"),
                TableHelper.getHeaderCell("Biallelic"),
                TableHelper.getHeaderCell("Driver")
        });

        for (int i = 0; i < 4; i++) {
            contentTable.addCell(TableHelper.getContentCell("BRAF*"));
            contentTable.addCell(TableHelper.getContentCell("c.1799T>A*"));
            contentTable.addCell(TableHelper.getContentCell("p.Val6000Glu"));
            contentTable.addCell(TableHelper.getContentCell("107 / 161").setTextAlignment(TextAlignment.RIGHT));
            contentTable.addCell(TableHelper.getContentCell("Yes"));
            contentTable.addCell(TableHelper.getContentCell("AAAABB (65%)"));
            contentTable.addCell(TableHelper.getContentCell(new InlineBarChart(.6f, 0.0f, 1.0f)));
            contentTable.addCell(TableHelper.getContentCell("Clonal"));
            contentTable.addCell(TableHelper.getContentCell("-"));
            contentTable.addCell(TableHelper.getContentCell("High"));
        }

        // Create report table that handles page breaks
        return TableHelper.createWrappingReportTable(chapterTitle, contentTable);

    }

    @NotNull
    private final Table createGainsAndLossesTable() {

        final String chapterTitle = "Tumor specific gains & losses";
        final boolean isAvailable = false;

        if (!isAvailable) {
            return createNoneTable(chapterTitle);
        }

        // Create content table
        Table contentTable = TableHelper.createReportContentTable(new float[] {4, 4, 4, 2, 1}, new Cell[]  {
                TableHelper.getHeaderCell("Chromosome"),
                TableHelper.getHeaderCell("Chromosome band"),
                TableHelper.getHeaderCell("Gene"),
                TableHelper.getHeaderCell("Type"),
                TableHelper.getHeaderCell("Copies").setTextAlignment(TextAlignment.RIGHT)
        });

        for (int i = 0; i < 4; i++) {
            contentTable.addCell(TableHelper.getContentCell("10"));
            contentTable.addCell(TableHelper.getContentCell("q23.31"));
            contentTable.addCell(TableHelper.getContentCell("PTEN"));
            contentTable.addCell(TableHelper.getContentCell("partial loss"));
            contentTable.addCell(TableHelper.getContentCell("0").setTextAlignment(TextAlignment.RIGHT));
        }

        // Create report table that handles page breaks
        return TableHelper.createWrappingReportTable(chapterTitle, contentTable);

    }

    @NotNull
    private final Table createSomaticFusionsTable() {

        final String chapterTitle = "Somatic gene fusions";
        final boolean isAvailable = true;

        if (!isAvailable) {
            return createNoneTable(chapterTitle);
        }

        // Create content table
        Table contentTable = TableHelper.createReportContentTable(new float[] {2, .8f, 1.5f, 1.5f, 3f}, new Cell[]  {
                TableHelper.getHeaderCell("Fusion"),
                TableHelper.getHeaderCell("Transcript"),
                TableHelper.getHeaderCell("Context"),
                TableHelper.getHeaderCell("Copies").setTextAlignment(TextAlignment.RIGHT),
                TableHelper.getHeaderCell("Source")
        });

        for (int i = 0; i < 4; i++) {
            contentTable.addCell(TableHelper.getContentCell("---"));
            contentTable.addCell(TableHelper.getContentCell("## - ##"));
            contentTable.addCell(TableHelper.getContentCell("## - ##"));
            contentTable.addCell(TableHelper.getContentCell("0").setTextAlignment(TextAlignment.RIGHT));
            contentTable.addCell(TableHelper.getContentCell("http://oncokb.org/"));
        }

        // Create report table that handles page breaks
        return TableHelper.createWrappingReportTable(chapterTitle, contentTable);

    }

    @NotNull
    private final Table createDisruptionsTable() {

        final String chapterTitle = "Tumor specific gene disruptions";
        final boolean isAvailable = true;

        if (!isAvailable) {
            return createNoneTable(chapterTitle);
        }

        // Create content table
        Table contentTable = TableHelper.createReportContentTable(new float[] {1, 1, 2, .8f, 1.5f, 1.5f, 1.5f}, new Cell[]  {
                TableHelper.getHeaderCell("Location"),
                TableHelper.getHeaderCell("Gene"),
                TableHelper.getHeaderCell("Disrupted range"),
                TableHelper.getHeaderCell("Type"),
                TableHelper.getHeaderCell("Copies").setTextAlignment(TextAlignment.RIGHT),
                TableHelper.getHeaderCell("Gene min copies").setTextAlignment(TextAlignment.RIGHT),
                TableHelper.getHeaderCell("Gene max copies").setTextAlignment(TextAlignment.RIGHT)
        });

        for (int i = 0; i < 4; i++) {
            contentTable.addCell(TableHelper.getContentCell("q23.31"));
            contentTable.addCell(TableHelper.getContentCell("PTEN"));
            contentTable.addCell(TableHelper.getContentCell("Intron 5 -> Intron 6"));
            contentTable.addCell(TableHelper.getContentCell("DEL"));
            contentTable.addCell(TableHelper.getContentCell("1.8").setTextAlignment(TextAlignment.RIGHT));
            contentTable.addCell(TableHelper.getContentCell("0").setTextAlignment(TextAlignment.RIGHT));
            contentTable.addCell(TableHelper.getContentCell("2").setTextAlignment(TextAlignment.RIGHT));
        }

        // Create report table that handles page breaks
        return TableHelper.createWrappingReportTable(chapterTitle, contentTable);

    }

    @NotNull
    private static Table createNoneTable(String sectionTitle) {
//
//        Table table = TableHelper.createReportContentTable(new float[] {1}, new Cell[] {});
//        table.addCell(TableHelper.getDisabledContentCell(new Paragraph("NONE")));
//        table.addCell(TableHelper.getLayoutCell(10, 1).setHeight(TABLE_BOTTOM_SPACER_HEIGHT));
//

        return TableHelper.createNoneReportTable(sectionTitle);


    }

}
