package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatPercentage;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.type;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.hartwig.hmftools.datamodel.isofox.RnaQCStatus;
import com.hartwig.hmftools.datamodel.isofox.RnaStatistics;
import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.jasperreports.engine.data.JRMapCollectionDataSource;

public final class RnaStatisticsTable
{
    public static JasperReportBuilder build(final String title, final RnaStatistics rnaStatistics,
            final ReportResources reportResources)
    {
        StringJoiner qcSj = new StringJoiner(", ");
        for(RnaQCStatus status : rnaStatistics.qcStatus())
        {
            qcSj.add(status.name());
        }

        double duplicateRate = rnaStatistics.duplicateFragments() / (double) rnaStatistics.totalFragments();

        Map<String, Object> row = new HashMap<>();
        row.put("qc", qcSj.toString());
        row.put("totalfrags", String.valueOf(rnaStatistics.totalFragments()));
        row.put("duprate", formatPercentage(duplicateRate));
        row.put("splicedrate", formatPercentage(rnaStatistics.splicedFragmentPerc()));
        row.put("unsplicedrate", formatPercentage(rnaStatistics.unsplicedFragmentPerc()));
        row.put("altrate", formatPercentage(rnaStatistics.altFragmentPerc()));
        row.put("chimericrate", formatPercentage(rnaStatistics.chimericFragmentPerc()));

        List<Map<String, ?>> rows = new ArrayList<>();
        rows.add(row);

        return report()
                .setColumnTitleStyle(OrangeFonts.TABLE_HEADER_STYLE_UNPADDED)
                .setColumnStyle(OrangeFonts.TABLE_CONTENT_STYLE_UNPADDED)
                .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE_WITH_GAP))
                .columns(
                        col.column("QC", "qc", type.stringType()).setWidth(10),
                        col.column("TOTAL FRAGMENTS", "totalfrags", type.stringType()).setWidth(10),
                        col.column("DUPLICATE RATE", "duprate", type.stringType()).setWidth(10),
                        col.column("SPLICED RATE", "splicedrate", type.stringType()).setWidth(10),
                        col.column("UNSPLICED RATE", "unsplicedrate", type.stringType()).setWidth(10),
                        col.column("ALT-SLICED RATE", "altrate", type.stringType()).setWidth(10),
                        col.column("CHIMERIC RATE", "chimericrate", type.stringType()).setWidth(10)
                )
                .setDataSource(new JRMapCollectionDataSource(rows));
    }
}
