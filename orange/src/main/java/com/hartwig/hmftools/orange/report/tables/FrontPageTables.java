package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;

import org.jetbrains.annotations.Nullable;

public class FrontPageTables
{
    public static Table buildSampleSummary(
            final Map<String, String> summary, @Nullable final String qcWarning, float width, final ReportResources reportResources)
    {
        Cells cells = new Cells(reportResources);
        Table table = buildSummaryTable(summary, width, cells);

        if(qcWarning != null)
        {
            table.addCell(cells.createSpanningWarning(table, qcWarning));
        }

        return table;
    }

    public static Table buildTechnicalSummary(
            final Map<String, String> summary, float width, final ReportResources reportResources)
    {
        return buildSummaryTable(summary, width, new Cells(reportResources));
    }

    public static Table buildDriverSummary(
            final Map<String, String> driverSummary, float width, final ReportResources reportResources)
    {
        return buildKeyValueTable(driverSummary, new float[] { 1, 2 }, "Driver Summary", reportResources);
    }

    public static Table buildGenomeWideFeatures(
            final Map<String, String> features, float width, final ReportResources reportResources)
    {
        return buildKeyValueTable(features, new float[] { 1.9f, 1.0f }, "Genome Wide Biomarkers", reportResources);
    }

    private static Table buildSummaryTable(final Map<String, String> summary, float width, Cells cells)
    {
        List<Integer> widths = new ArrayList<>();
        List<Cell> cellEntries = new ArrayList<>();

        for(String key : summary.keySet())
        {
            addEntry(cells, widths, cellEntries, 2, key);
        }

        Table table = Tables.createContent(width, intToFloatArray(widths), cellArray(cellEntries));

        for(String value : summary.values())
        {
            table.addCell(cells.createContent(value));
        }

        return table;
    }

    private static Table buildKeyValueTable(
            final Map<String, String> entries, float[] columnWidths, String title, final ReportResources reportResources)
    {
        Cells cells = new Cells(reportResources);
        Table table = new Table(UnitValue.createPercentArray(columnWidths));

        for(Map.Entry<String, String> entry : entries.entrySet())
        {
            table.addCell(cells.createKey(entry.getKey()));
            table.addCell(cells.createValue(entry.getValue()));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }
}
