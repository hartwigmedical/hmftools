package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.data.Util;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.layout.borders.Border;
import com.itextpdf.layout.borders.SolidBorder;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.IBlockElement;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.VerticalAlignment;
import org.jetbrains.annotations.NotNull;

public final class TableHelper {

    private final static float TABLE_BOTTOM_MARGIN = 20;

    /**
     * Get a table that implements the visual style for the main report tables
     */
    @NotNull
    public static Table createReportContentTable(float[] columnPercentageWidths, Cell[] headerCells) {

        Table table = new Table(columnPercentageWidths)
                .setWidth(ReportResources.CONTENT_WIDTH_WIDE);

        for (Cell headerCell: headerCells) {
            table.addHeaderCell(headerCell);
        }

        return table;

    }

    @NotNull
    public static Table createNoneReportTable(@NotNull String tableTitle) {

        Cell headerCell = new Cell()
                .setBorder(Border.NO_BORDER)
                .add(new Paragraph(tableTitle)
                        .addStyle(ReportResources.sectionTitleStyle()
                                .setFontColor(ReportResources.PALETTE_LIGHT_GREY)));

        Table table = TableHelper.createReportContentTable(new float[] {1}, new Cell[] { headerCell });
        table.setKeepTogether(true);
        table.setMarginBottom(TABLE_BOTTOM_MARGIN);
        table.addCell(TableHelper.getDisabledContentCell(new Paragraph(Util.NoneString)));

        return table;

    }

    public static Table createWrappingReportTable(String tableTitle, Table contentTable)  {

        // Add "continues on next page" footer to the content table
        contentTable.addFooterCell(new Cell(1, contentTable.getNumberOfColumns())
                .setBorder(Border.NO_BORDER)
                .setPaddingTop(5)
                .setPaddingBottom(5)
                .add(new Paragraph("The table continues on the next page".toUpperCase())
                        .addStyle(ReportResources.subTextStyle())))
                .setSkipLastFooter(true);

        // Wrap content table with a table that shows "continued from the previous table" after page break
        Table continuedWrapTable = new Table(1)
                .setWidth(contentTable.getWidth())
                .addHeaderCell(new Cell()
                        .setBorder(Border.NO_BORDER)
                        .add(new Paragraph("Continued from the previous page".toUpperCase())
                                .addStyle(ReportResources.subTextStyle())))
                .setSkipFirstHeader(true)
                .addCell(new Cell()
                        .add(contentTable)
                        .setPadding(0)
                        .setBorder(Border.NO_BORDER));

        // Wrap continuedWrapTable with table that shows the table title
        return new Table(1)
                .setWidth(contentTable.getWidth())
                .setMarginBottom(TABLE_BOTTOM_MARGIN)
                .addHeaderCell(new Cell()
                        .setBorder(Border.NO_BORDER)
                        .setPadding(0)
                        .add(new Paragraph(tableTitle)
                                .addStyle(ReportResources.sectionTitleStyle())))
                .addCell(new Cell()
                        .add(continuedWrapTable)
                        .setPadding(0)
                        .setBorder(Border.NO_BORDER));

    }


    /**
     * Get a cell that implements the visual header style for the main report tables
     */
    @NotNull
    public static Cell getHeaderCell(@NotNull String text) {
        return getHeaderCell(text, 1);
    }

    /**
     * Get a cell that implements the visual header style for the main report tables
     */
    @NotNull
    public static Cell getHeaderCell(@NotNull String text, int colspan) {
        return getHeaderCell(colspan)
                .add(new Paragraph(text.toUpperCase()));
    }

    /**
     * Get a cell that implements the visual header style for the main report tables
     */
    @NotNull
    public static Cell getHeaderCell() {
        return getHeaderCell(1);
    }

    /**
     * Get a cell that implements the visual header style for the main report tables
     */
    @NotNull
    public static Cell getHeaderCell(int colspan) {
        Cell c = new Cell(1, colspan);
        c.setHeight(23); // Set fixed height to create consistent spacing between table title and header
        c.setBorder(Border.NO_BORDER);
        c.setVerticalAlignment(VerticalAlignment.BOTTOM);
        c.addStyle(ReportResources.tableHeaderStyle());
        return c;
    }

    /**
     * Get a cell that implements the visual content style for the main report tables
     */
    @NotNull
    public static Cell getContentCell(@NotNull String text) {
        return getContentCell(new Paragraph(text));
    }

    /**
     * Get a cell that implements the visual content style for the main report tables
     */
    @NotNull
    public static Cell getContentCell(@NotNull IBlockElement element) {
        Cell c = new Cell();
        c.setBorder(Border.NO_BORDER);
        c.setBorderBottom(new SolidBorder(ReportResources.PALETTE_MID_GREY, 0.25f));
        c.addStyle(ReportResources.tableContentStyle());
        c.setKeepTogether(true);
        c.add(element);
        return c;
    }

    /**
     * Get a cell that implements the visual content disabled style for the main report tables
     */
    @NotNull
    private static Cell getDisabledContentCell(@NotNull IBlockElement element) {
        Cell c = new Cell();
        c.setBorder(Border.NO_BORDER);
        c.setBorderBottom(new SolidBorder(ReportResources.PALETTE_LIGHT_GREY, 0.25f));
        c.addStyle(ReportResources.tableContentStyle().setFontColor(ReportResources.PALETTE_LIGHT_GREY));
        c.setKeepTogether(true);
        c.add(element);
        return c;
    }

    /**
     * Get a cell that has no borders (just for positioning elements)
     */
    @NotNull
    public static Cell getLayoutCell() {
        return getLayoutCell(1, 1);
    }

    /**
     * Get a cell that has no borders (just for positioning elements)
     */
    @NotNull
    public static Cell getLayoutCell(int rowspan, int colspan) {

        Cell c = new Cell(rowspan, colspan);
        c.setBorder(Border.NO_BORDER);
        c.setKeepTogether(true);
        c.setPadding(0);
        c.setMargin(0);
        return c;
    }

}
