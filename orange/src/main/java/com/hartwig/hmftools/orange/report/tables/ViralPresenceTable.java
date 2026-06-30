package com.hartwig.hmftools.orange.report.tables;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.type;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;
import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.jasperreports.engine.JREmptyDataSource;
import net.sf.jasperreports.engine.data.JRMapCollectionDataSource;

public final class ViralPresenceTable
{
    public static JasperReportBuilder build(
            final String title, final List<VirusInterpreterEntry> viruses, final ReportResources reportResources)
    {
        if(viruses.isEmpty())
        {
            return report()
                    .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE),
                            cmp.text(ReportResources.NONE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE))
                    .setDataSource(new JREmptyDataSource());
        }

        List<Map<String, ?>> rows = new ArrayList<>();
        for(VirusInterpreterEntry virus : viruses)
        {
            Map<String, Object> row = new LinkedHashMap<>();
            row.put("virus", virus.name());
            row.put("qcstatus", virus.qcStatus().toString());
            row.put("type", virus.interpretation() != null ? virus.interpretation().name() : "");
            row.put("integrations", String.valueOf(virus.integrations()));
            row.put("perccov", TableCommon.formatPercentage(virus.percentageCovered(), false));
            row.put("meancov", TableCommon.formatSingleDigitDecimal(virus.meanCoverage()));
            row.put("clonalcov", virus.expectedClonalCoverage() != null
                    ? TableCommon.formatSingleDigitDecimal(virus.expectedClonalCoverage()) : "");
            row.put("driver", virus.driverInterpretation().toString());
            rows.add(row);
        }

        return report()
                .setColumnTitleStyle(OrangeFonts.TABLE_HEADER_STYLE_UNPADDED)
                .setColumnStyle(OrangeFonts.TABLE_CONTENT_STYLE_UNPADDED)
                .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE_WITH_GAP))
                .columns(
                        col.column("VIRUS", "virus", type.stringType()).setWidth(25),
                        col.column("QC STATUS", "qcstatus", type.stringType()).setWidth(25),
                        col.column("TYPE", "type", type.stringType()).setWidth(10),
                        col.column("INTEGRATIONS", "integrations", type.stringType()).setWidth(20),
                        col.column("% COVERED", "perccov", type.stringType()).setWidth(15),
                        col.column("MEAN COV", "meancov", type.stringType()).setWidth(15),
                        col.column("EXP CLONAL COV", "clonalcov", type.stringType()).setWidth(20),
                        col.column("DRIVER", "driver", type.stringType()).setWidth(10)
                )
                .setDataSource(new JRMapCollectionDataSource(rows));
    }
}
