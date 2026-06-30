package com.hartwig.hmftools.orange.report.tables;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.stl;
import static net.sf.dynamicreports.report.builder.DynamicReports.type;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.orange.report.OrangeColors;
import com.hartwig.hmftools.orange.report.OrangeFonts;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.column.TextColumnBuilder;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;
import net.sf.jasperreports.engine.data.JRMapCollectionDataSource;

import org.jetbrains.annotations.Nullable;

public final class FrontPageTables
{
    public static JasperReportBuilder buildSampleSummary(
            final Map<String, String> summary, @Nullable final String qcWarning)
    {
        JasperReportBuilder table = buildHorizontalHeaderTable(summary);

        if(qcWarning != null)
        {
            return report()
                    .summary(
                            cmp.subreport(table),
                            cmp.text(qcWarning).setStyle(OrangeFonts.QC_WARNING_STYLE)
                    )
                    .setDataSource(new JRMapCollectionDataSource(List.of()));
        }

        return table;
    }

    public static JasperReportBuilder buildTechnicalSummary(final Map<String, String> summary)
    {
        return buildHorizontalHeaderTable(summary);
    }

    public static VerticalListBuilder buildDriverSummaryBox(final Map<String, String> driverSummary, int rows)
    {
        return buildKeyValueBox("Driver Summary", driverSummary, rows);
    }

    public static VerticalListBuilder buildGenomeWideFeaturesBox(final Map<String, String> features, int rows)
    {
        return buildKeyValueBox("Genome Wide Biomarkers", features, rows);
    }

    private static JasperReportBuilder buildHorizontalHeaderTable(final Map<String, String> summary)
    {
        List<TextColumnBuilder<String>> columns = new ArrayList<>();
        Map<String, Object> rowData = new LinkedHashMap<>();

        for(Map.Entry<String, String> entry : summary.entrySet())
        {
            String key = entry.getKey().toUpperCase();
            String fieldKey = key.replace(" ", "_").toLowerCase();
            columns.add(col.column(key, fieldKey, type.stringType()));
            rowData.put(fieldKey, entry.getValue());
        }

        JasperReportBuilder table = report()
                .setColumnTitleStyle(OrangeFonts.TABLE_HEADER_STYLE_UNPADDED)
                .setColumnStyle(OrangeFonts.TABLE_CONTENT_STYLE_UNPADDED);

        for(TextColumnBuilder<String> column : columns)
        {
            table.addColumn(column);
        }

        return table.setDataSource(new JRMapCollectionDataSource(List.of(rowData)));
    }

    // Matches Kotlin's createTwoColumnTablePaddedWithEmptyRows in Helpers.kt.
    // Both driver summary and genome wide features boxes are padded to the same row count
    // so that their heights match and borders look symmetrical.
    private static VerticalListBuilder buildKeyValueBox(final String title, final Map<String, String> data, int rows)
    {
        VerticalListBuilder container = cmp.verticalList();
        container.add(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE));
        container.add(cmp.verticalGap(5));

        List<Map.Entry<String, String>> entries = new ArrayList<>(data.entrySet());
        for(int i = 0; i < rows; i++)
        {
            String key = i < entries.size() ? entries.get(i).getKey() : "";
            String value = i < entries.size() ? entries.get(i).getValue() : "";
            container.add(cmp.horizontalList(
                    cmp.text(key).setStyle(OrangeFonts.TABLE_CONTENT_STYLE),
                    cmp.text(value).setStyle(OrangeFonts.TABLE_CONTENT_STYLE)
            ));
        }

        container.setStyle(stl.style()
                .setBorder(stl.border(stl.penThin().setLineColor(OrangeColors.PALETTE_LIGHT_GREY)))
                .setPadding(stl.padding(10)));

        return container;
    }
}
