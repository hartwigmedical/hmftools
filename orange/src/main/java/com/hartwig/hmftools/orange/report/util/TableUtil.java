package com.hartwig.hmftools.orange.report.util;

import com.hartwig.hmftools.common.utils.DataUtil;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.layout.borders.Border;
import com.itextpdf.layout.borders.SolidBorder;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.IBlockElement;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;
import com.itextpdf.layout.property.VerticalAlignment;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class TableUtil {

    private static final float TABLE_BOTTOM_MARGIN = 20;

    private TableUtil() {
    }

    @NotNull
    public static Table createEmptyTable(@NotNull String title, float width) {
        Cell headerCell = new Cell().setBorder(Border.NO_BORDER).add(new Paragraph(title).addStyle(ReportResources.tableTitleStyle()));

        Table table = TableUtil.createReportContentTable(width, new float[] { 1 }, new Cell[] { headerCell });
        table.setKeepTogether(true);
        table.setMarginBottom(TABLE_BOTTOM_MARGIN);
        table.addCell(TableUtil.createContentCell(new Paragraph(DataUtil.NONE_STRING)));

        return table;
    }

    @NotNull
    public static Table createReportContentTable(float width, @NotNull float[] columnPercentageWidths, @NotNull Cell[] headerCells) {
        Table table = new Table(UnitValue.createPercentArray(columnPercentageWidths)).setWidth(width);
        table.setFixedLayout();

        for (Cell headerCell : headerCells) {
            table.addHeaderCell(headerCell);
        }

        return table;
    }

    @NotNull
    public static Table createWrappingReportTable(@NotNull Table contentTable) {
        return createWrappingReportTable(contentTable, null);
    }

    @NotNull
    public static Table createWrappingReportTable(@NotNull Table contentTable, @Nullable String title) {
        contentTable.addFooterCell(new Cell(1, contentTable.getNumberOfColumns()).setBorder(Border.NO_BORDER)
                .setPaddingTop(5)
                .setPaddingBottom(5)
                .add(new Paragraph("The table continues on the next page".toUpperCase()).addStyle(ReportResources.subTextStyle())))
                .setSkipLastFooter(true);

        Table continuedWrapTable = new Table(1).setMinWidth(contentTable.getWidth())
                .addHeaderCell(new Cell().setBorder(Border.NO_BORDER)
                        .add(new Paragraph("Continued from the previous page".toUpperCase()).addStyle(ReportResources.subTextStyle())))
                .setSkipFirstHeader(true)
                .addCell(new Cell().add(contentTable).setPadding(0).setBorder(Border.NO_BORDER));

        Table table = new Table(1).setMinWidth(contentTable.getWidth()).setMarginBottom(TABLE_BOTTOM_MARGIN);
        if (title != null) {
            table.addHeaderCell(new Cell().setBorder(Border.NO_BORDER)
                    .setPadding(0)
                    .add(new Paragraph(title).addStyle(ReportResources.tableTitleStyle())));
        }

        table.addCell(new Cell().add(continuedWrapTable).setPadding(0).setBorder(Border.NO_BORDER));
        return table;
    }

    @NotNull
    public static Cell createHeaderCell(@NotNull String text) {
        return createHeaderCell(text, 1);
    }

    @NotNull
    public static Cell createHeaderCell(@NotNull String text, int colSpan) {
        return createHeaderCell(colSpan).add(new Paragraph(text.toUpperCase()));
    }

    @NotNull
    private static Cell createHeaderCell(int colSpan) {
        Cell cell = new Cell(1, colSpan);
        cell.setHeight(23); // Set fixed height to create consistent spacing between table title and header
        cell.setBorder(Border.NO_BORDER);
        cell.setVerticalAlignment(VerticalAlignment.BOTTOM);
        cell.addStyle(ReportResources.tableHeaderStyle());
        return cell;
    }

    @NotNull
    public static Cell createTransparentCell(@NotNull String text) {
        return createTransparentCell(new Paragraph(text));
    }

    @NotNull
    public static Cell createTransparentCell(@NotNull IBlockElement element) {
        Cell cell = new Cell();
        cell.setBorder(Border.NO_BORDER);
        cell.setBorderBottom(Border.NO_BORDER);
        cell.addStyle(ReportResources.tableContentStyle());
        cell.setKeepTogether(true);
        cell.add(element);
        return cell;
    }

    @NotNull
    public static Cell createContentCell(@NotNull String text) {
        return createContentCell(new Paragraph(text));
    }

    @NotNull
    public static Cell createContentCell(@NotNull IBlockElement element) {
        Cell cell = new Cell();
        cell.setBorder(Border.NO_BORDER);
        cell.setBorderBottom(new SolidBorder(ReportResources.PALETTE_MID_GREY, 0.25F));
        cell.addStyle(ReportResources.tableContentStyle());
        cell.setKeepTogether(true);
        cell.add(element);
        return cell;
    }

    @NotNull
    public static Cell createImageCell(@NotNull Image image) {
        Cell cell = new Cell();
        cell.setBorder(Border.NO_BORDER);
        cell.setBorderBottom(new SolidBorder(ReportResources.PALETTE_MID_GREY, 0.25F));
        cell.addStyle(ReportResources.tableContentStyle());
        cell.setKeepTogether(true);
        cell.add(image);
        return cell;
    }

    @NotNull
    public static Cell createKeyCell(@NotNull String text) {
        return createKeyCell(new Paragraph(text));
    }

    @NotNull
    public static Cell createKeyCell(@NotNull IBlockElement element) {
        Cell cell = new Cell();
        cell.setBorder(Border.NO_BORDER);
        cell.addStyle(ReportResources.tableHeaderStyle());
        cell.setKeepTogether(true);
        cell.add(element);
        return cell;
    }

    @NotNull
    public static Cell createValueCell(@NotNull String text) {
        return createValueCell(new Paragraph(text));
    }

    @NotNull
    public static Cell createValueCell(@NotNull IBlockElement element) {
        Cell cell = new Cell();
        cell.setBorder(Border.NO_BORDER);
        cell.addStyle(ReportResources.tableContentStyle());
        cell.setKeepTogether(true);
        cell.add(element);
        return cell;
    }
}
