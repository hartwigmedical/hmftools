package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.layout.borders.Border;
import com.itextpdf.layout.borders.SolidBorder;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.IBlockElement;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.HorizontalAlignment;
import com.itextpdf.layout.property.TextAlignment;
import com.itextpdf.layout.property.UnitValue;
import com.itextpdf.layout.property.VerticalAlignment;
import org.jetbrains.annotations.NotNull;

public final class TableHelper {

    /**
     * Get a table that implements the visual style for the main report tables
     *
     * @param columnPercentageWidths
     * @return
     */
    @NotNull
    public static Table getReportTable(float[] columnPercentageWidths) {
        return new Table(UnitValue.createPercentArray(columnPercentageWidths))
                .setWidth(ReportResources.CONTENT_WIDTH_WIDE);
    }

    /**
     * Get a table that implements the visual style for the main report tables
     *
     * @param columnPercentageWidths
     * @param headerNames
     * @return
     */
    @NotNull
    public static Table getReportTable(float[] columnPercentageWidths, @NotNull String[] headerNames) {

        // Initialize table
        Table table = getReportTable(columnPercentageWidths);

        // Add header columns
        for (String name: headerNames) {
            table.addHeaderCell(getHeaderCell(name));
        }

        return table;

    }

    /**
     * Get a cell that implements the visual header style for the main report tables
     *
     * @return
     */
    @NotNull
    public static Cell getHeaderCell(@NotNull String text) {
        Cell c = getHeaderCell();
        c.add(new Paragraph(text.toUpperCase()));
        return c;
    }

    /**
     * Get a cell that implements the visual header style for the main report tables
     *
     * @return
     */
    @NotNull
    public static Cell getHeaderCell() {
        Cell c = new Cell();
        c.setBorder(Border.NO_BORDER);
        c.setVerticalAlignment(VerticalAlignment.BOTTOM);
        c.addStyle(ReportResources.tableHeaderStyle());
        return c;
    }

    /**
     * Get a cell that implements the visual content style for the main report tables
     *
     * @param text
     * @return
     */
    @NotNull
    public static Cell getContentCell(@NotNull String text) {
        return getContentCell(new Paragraph(text));
    }

    /**
     * Get a cell that implements the visual content style for the main report tables
     *
     * @param element
     * @return
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
     *
     * @param element
     * @return
     */
    @NotNull
    public static Cell getDisabledContentCell(@NotNull IBlockElement element) {
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
     *
     * @return
     */
    @NotNull
    public static Cell getLayoutCell() {
        return getLayoutCell(1, 1);
    }

    /**
     * Get a cell that has no borders (just for positioning elements)
     *
     * @return
     */
    @NotNull
    public static Cell getLayoutCell(int rowspan, int colspan) {

        Cell c = new Cell(rowspan, colspan);
        c.setBorder(Border.NO_BORDER);
        c.setKeepTogether(true);
        return c;
    }

}
