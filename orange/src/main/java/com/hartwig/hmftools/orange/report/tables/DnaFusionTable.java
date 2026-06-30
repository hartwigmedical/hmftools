package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.type;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.linx.FusionPhasedType;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.jasperreports.engine.JREmptyDataSource;
import net.sf.jasperreports.engine.data.JRMapCollectionDataSource;

public final class DnaFusionTable
{
    public static JasperReportBuilder build(
            final String title, final List<LinxFusion> fusions, final ReportResources reportResources)
    {
        if(fusions.isEmpty())
        {
            return report()
                    .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE),
                            cmp.text(ReportResources.NONE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE))
                    .setDataSource(new JREmptyDataSource());
        }

        boolean hasRna = fusions.stream().anyMatch(f -> f.rnaSupport() != null);

        List<Map<String, ?>> rows = new ArrayList<>();
        for(LinxFusion fusion : sort(fusions))
        {
            Map<String, Object> row = new LinkedHashMap<>();
            row.put("fusion", fusionDisplay(fusion));
            row.put("junctions", transcriptJunctions(fusion));
            row.put("jcn", TableCommon.formatSingleDigitDecimal(fusion.junctionCopyNumber()));
            row.put("phasing", display(fusion.phased()));
            row.put("type", fusion.reportedType().toString());
            if(hasRna)
            {
                row.put("rnafrags", fusion.rnaSupport() != null
                        ? String.valueOf(fusion.rnaSupport().alleleReadCount())
                        : ReportResources.NOT_AVAILABLE);
            }
            row.put("driver", fusion.driverInterpretation().toString());
            rows.add(row);
        }

        JasperReportBuilder table = report()
                .setColumnTitleStyle(OrangeFonts.TABLE_HEADER_STYLE_UNPADDED)
                .setColumnStyle(OrangeFonts.TABLE_CONTENT_STYLE_UNPADDED)
                .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE_WITH_GAP))
                .columns(
                        col.column("FUSION", "fusion", type.stringType()).setWidth(20),
                        col.column("JUNCTIONS", "junctions", type.stringType()).setWidth(50),
                        col.column("JCN", "jcn", type.stringType()).setWidth(10),
                        col.column("PHASING", "phasing", type.stringType()).setWidth(10),
                        col.column("TYPE", "type", type.stringType()).setWidth(30)
                );

        if(hasRna)
        {
            table.addColumn(col.column("RNA FRAGS", "rnafrags", type.stringType()).setWidth(10));
        }

        table.addColumn(col.column("DRIVER", "driver", type.stringType()).setWidth(10));
        table.setDataSource(new JRMapCollectionDataSource(rows));
        return table;
    }

    private static String fusionDisplay(final LinxFusion fusion)
    {
        return format("%s::%s", fusion.geneUp(), fusion.geneDown());
    }

    private static String transcriptJunctions(final LinxFusion fusion)
    {
        return format("%s (%s) - %s (%s)",
                fusion.contextUp(), fusion.transcriptUp(), fusion.contextDown(), fusion.transcriptDown());
    }

    private static String display(final FusionPhasedType type)
    {
        switch(type)
        {
            case INFRAME:
                return "Inframe";
            case SKIPPED_EXONS:
                return "Skipped exons";
            case OUT_OF_FRAME:
                return "Out of frame";
        }
        throw new IllegalStateException("Unknown FusionPhasedType: " + type);
    }

    private static List<LinxFusion> sort(final List<LinxFusion> fusions)
    {
        return fusions.stream().sorted((f1, f2) ->
        {
            if(f1.driverInterpretation() == f2.driverInterpretation())
            {
                if(f1.geneUp().equals(f2.geneUp()))
                {
                    return f1.geneDown().compareTo(f2.geneDown());
                }
                return f1.geneUp().compareTo(f2.geneUp());
            }
            return f1.driverInterpretation() == DriverInterpretation.HIGH ? -1 : 1;
        }).collect(Collectors.toList());
    }
}
