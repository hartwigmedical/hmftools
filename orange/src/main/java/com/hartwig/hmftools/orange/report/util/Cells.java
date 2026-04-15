package com.hartwig.hmftools.orange.report.util;

import com.hartwig.hmftools.orange.report.ReportResources;

import org.apache.pdfbox.pdmodel.PDPage;

import be.quodlibet.boxable.BaseTable;
import be.quodlibet.boxable.Cell;
import be.quodlibet.boxable.Row;
import be.quodlibet.boxable.VerticalAlignment;
import be.quodlibet.boxable.line.LineStyle;

public class Cells
{
    // Compact padding to match iText's tight row spacing
    private static final float CELL_PADDING_TOP = 1f;
    private static final float CELL_PADDING_BOTTOM = 1f;
    private static final float CELL_PADDING_LEFT = 2f;
    private static final float CELL_PADDING_RIGHT = 2f;

    public static final float COMPACT_ROW_HEIGHT = 10f;

    private final ReportResources mReportResources;

    public Cells(final ReportResources reportResources)
    {
        mReportResources = reportResources;
    }

    private static void setCompactPadding(final Cell<PDPage> cell)
    {
        cell.setTopPadding(CELL_PADDING_TOP);
        cell.setBottomPadding(CELL_PADDING_BOTTOM);
        cell.setLeftPadding(CELL_PADDING_LEFT);
        cell.setRightPadding(CELL_PADDING_RIGHT);
    }

    public void applyHeaderStyle(final Cell<PDPage> cell)
    {
        ReportResources.TextStyle style = mReportResources.tableHeaderStyle();
        cell.setFont(style.font());
        cell.setFontSize(style.fontSize());
        cell.setTextColor(style.color());
        cell.setTopBorderStyle(null);
        cell.setBottomBorderStyle(null);
        cell.setLeftBorderStyle(null);
        cell.setRightBorderStyle(null);
        cell.setValign(VerticalAlignment.BOTTOM);
        setCompactPadding(cell);
    }

    public void applyContentStyle(final Cell<PDPage> cell)
    {
        ReportResources.TextStyle style = mReportResources.tableContentStyle();
        cell.setFont(style.font());
        cell.setFontSize(style.fontSize());
        cell.setTextColor(style.color());
        cell.setTopBorderStyle(null);
        cell.setLeftBorderStyle(null);
        cell.setRightBorderStyle(null);
        cell.setBottomBorderStyle(new LineStyle(ReportResources.PALETTE_MID_GREY, 0.25f));
        setCompactPadding(cell);
    }

    public void applyKeyStyle(final Cell<PDPage> cell)
    {
        ReportResources.TextStyle style = mReportResources.keyStyle();
        cell.setFont(style.font());
        cell.setFontSize(style.fontSize());
        cell.setTextColor(style.color());
        cell.setTopBorderStyle(null);
        cell.setBottomBorderStyle(null);
        cell.setLeftBorderStyle(null);
        cell.setRightBorderStyle(null);
        setCompactPadding(cell);
    }

    public void applyValueStyle(final Cell<PDPage> cell)
    {
        ReportResources.TextStyle style = mReportResources.valueStyle();
        cell.setFont(style.font());
        cell.setFontSize(style.fontSize());
        cell.setTextColor(style.color());
        cell.setTopBorderStyle(null);
        cell.setBottomBorderStyle(null);
        cell.setLeftBorderStyle(null);
        cell.setRightBorderStyle(null);
        setCompactPadding(cell);
    }

    public void applyUrlStyle(final Cell<PDPage> cell)
    {
        ReportResources.TextStyle style = mReportResources.urlStyle();
        cell.setFont(style.font());
        cell.setFontSize(style.fontSize());
        cell.setTextColor(style.color());
        cell.setTopBorderStyle(null);
        cell.setLeftBorderStyle(null);
        cell.setRightBorderStyle(null);
        cell.setBottomBorderStyle(new LineStyle(ReportResources.PALETTE_MID_GREY, 0.25f));
        setCompactPadding(cell);
    }

    public void applyWarningStyle(final Cell<PDPage> cell)
    {
        ReportResources.TextStyle style = mReportResources.qcWarningStyle();
        cell.setFont(style.font());
        cell.setFontSize(style.fontSize());
        cell.setTextColor(style.color());
        cell.setTopBorderStyle(null);
        cell.setBottomBorderStyle(null);
        cell.setLeftBorderStyle(null);
        cell.setRightBorderStyle(null);
        setCompactPadding(cell);
    }

    public void applyTitleStyle(final Cell<PDPage> cell)
    {
        ReportResources.TextStyle style = mReportResources.tableTitleStyle();
        cell.setFont(style.font());
        cell.setFontSize(style.fontSize());
        cell.setTextColor(style.color());
        cell.setTopBorderStyle(null);
        cell.setBottomBorderStyle(null);
        cell.setLeftBorderStyle(null);
        cell.setRightBorderStyle(null);
        setCompactPadding(cell);
    }

    // ----- Convenience methods for adding styled cells to rows -----
    public Cell<PDPage> addHeaderCell(final Row<PDPage> row, float widthPct, final String text)
    {
        Cell<PDPage> cell = row.createCell(widthPct, text.toUpperCase());
        applyHeaderStyle(cell);
        return cell;
    }

    public Cell<PDPage> addContentCell(final Row<PDPage> row, float widthPct, final String text)
    {
        Cell<PDPage> cell = row.createCell(widthPct, text);
        applyContentStyle(cell);
        return cell;
    }

    public Cell<PDPage> addKeyCell(final Row<PDPage> row, float widthPct, final String text)
    {
        Cell<PDPage> cell = row.createCell(widthPct, text);
        applyKeyStyle(cell);
        return cell;
    }

    public Cell<PDPage> addValueCell(final Row<PDPage> row, float widthPct, final String text)
    {
        Cell<PDPage> cell = row.createCell(widthPct, text);
        applyValueStyle(cell);
        return cell;
    }

    public Cell<PDPage> addUrlCell(final Row<PDPage> row, float widthPct, final String text, final String url)
    {
        Cell<PDPage> cell = row.createCell(widthPct, text);
        applyUrlStyle(cell);
        // Note: Boxable doesn't natively support hyperlinks on cells.
        // The URL text is displayed; actual link annotation would need custom handling.
        return cell;
    }

    public Cell<PDPage> addWarningCell(final Row<PDPage> row, final String text)
    {
        Cell<PDPage> cell = row.createCell(100, text);
        applyWarningStyle(cell);
        return cell;
    }

    public java.util.List<Cell<PDPage>> addRow(
            final BaseTable table, final float[] pcts, final java.util.List<String> values)
    {
        Row<PDPage> row = table.createRow(COMPACT_ROW_HEIGHT);
        java.util.List<Cell<PDPage>> createdCells = new java.util.ArrayList<>();
        for(int i = 0; i < values.size(); i++)
        {
            float p = i < pcts.length ? pcts[i] : pcts[pcts.length - 1];
            createdCells.add(addContentCell(row, p, values.get(i)));
        }
        return createdCells;
    }

    public String createContent(final String text)
    {
        return text;
    }
}
