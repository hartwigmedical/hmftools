package com.hartwig.hmftools.orange.report.util;

import java.io.IOException;

import com.hartwig.hmftools.orange.report.DocumentContext;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.apache.pdfbox.pdmodel.PDPage;

import be.quodlibet.boxable.BaseTable;
import be.quodlibet.boxable.Cell;
import be.quodlibet.boxable.Row;

public class Tables
{
    private static final float ROW_HEIGHT = Cells.COMPACT_ROW_HEIGHT;
    private static final float HEADER_ROW_HEIGHT = 16f;

    private final ReportResources mReportResources;

    public Tables(final ReportResources reportResources)
    {
        mReportResources = reportResources;
    }

    public BaseTable createEmpty(final DocumentContext docCtx, final String title, float width) throws IOException
    {
        return createNonContent(docCtx, title, width, ReportResources.NONE);
    }

    public BaseTable createNonContent(final DocumentContext docCtx, final String title, float width, final String value) throws IOException
    {
        Cells cells = new Cells(mReportResources);

        BaseTable table = docCtx.createTable(width, null);

        // Title row
        Row<PDPage> titleRow = table.createRow(HEADER_ROW_HEIGHT);
        Cell<PDPage> titleCell = titleRow.createCell(100, title);
        cells.applyTitleStyle(titleCell);

        // Content row
        Row<PDPage> row = table.createRow(ROW_HEIGHT);
        cells.addContentCell(row, 100, value);

        return table;
    }

    public BaseTable createContent(
            final DocumentContext docCtx, float width, final float[] columnPercentageWidths,
            final String[] headerTexts) throws IOException
    {
        Cells cells = new Cells(mReportResources);

        BaseTable table = docCtx.createTable(width, columnPercentageWidths);

        // Normalize widths to percentages
        float totalWidth = 0;
        for(float w : columnPercentageWidths)
        {
            totalWidth += w;
        }

        Row<PDPage> headerRow = table.createRow(HEADER_ROW_HEIGHT);
        for(int i = 0; i < headerTexts.length; i++)
        {
            float pct = (columnPercentageWidths[i] / totalWidth) * 100f;
            cells.addHeaderCell(headerRow, pct, headerTexts[i]);
        }
        table.addHeaderRow(headerRow);

        return table;
    }

    public BaseTable createWithTitle(
            final DocumentContext docCtx, final String title, float width, final float[] columnPercentageWidths,
            final String[] headerTexts) throws IOException
    {
        Cells cells = new Cells(mReportResources);

        BaseTable table = docCtx.createTable(width, columnPercentageWidths);

        // Title row spanning full width
        Row<PDPage> titleRow = table.createRow(HEADER_ROW_HEIGHT);
        Cell<PDPage> titleCell = titleRow.createCell(100, title);
        cells.applyTitleStyle(titleCell);

        // Header row
        float totalWidth = 0;
        for(float w : columnPercentageWidths)
        {
            totalWidth += w;
        }

        Row<PDPage> headerRow = table.createRow(HEADER_ROW_HEIGHT);
        for(int i = 0; i < headerTexts.length; i++)
        {
            float pct = (columnPercentageWidths[i] / totalWidth) * 100f;
            cells.addHeaderCell(headerRow, pct, headerTexts[i]);
        }
        table.addHeaderRow(headerRow);

        return table;
    }

    public static float rowHeight()
    {
        return ROW_HEIGHT;
    }
}
