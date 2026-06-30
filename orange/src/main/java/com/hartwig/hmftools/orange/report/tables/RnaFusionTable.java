package com.hartwig.hmftools.orange.report.tables;

import static java.lang.Math.round;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSupportField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.zeroPrefixed;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.type;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.isofox.RnaFusion;
import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.jasperreports.engine.JREmptyDataSource;
import net.sf.jasperreports.engine.data.JRMapCollectionDataSource;

public final class RnaFusionTable
{
    public static JasperReportBuilder build(final String title, final List<RnaFusion> fusions,
            final ReportResources reportResources)
    {
        if(fusions.isEmpty())
        {
            return report()
                    .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE),
                            cmp.text(ReportResources.NONE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE))
                    .setDataSource(new JREmptyDataSource());
        }

        List<Map<String, ?>> rows = new ArrayList<>();
        for(RnaFusion fusion : sort(fusions))
        {
            int averageDepth = (int) round((fusion.depthStart() + fusion.depthEnd()) * 0.5);

            Map<String, Object> row = new HashMap<>();
            row.put("fusion", fusion.display());
            row.put("type", fusion.knownType().toString());
            row.put("svtype", fusion.svType().toString());
            row.put("juncstart", fusion.junctionTypeStart());
            row.put("juncend", fusion.junctionTypeEnd());
            row.put("support", formatSupportField(fusion.splitFragments(), averageDepth));
            row.put("cohortfreq", String.valueOf(fusion.cohortFrequency()));
            rows.add(row);
        }

        return report()
                .setColumnTitleStyle(OrangeFonts.TABLE_HEADER_STYLE_UNPADDED)
                .setColumnStyle(OrangeFonts.TABLE_CONTENT_STYLE_UNPADDED)
                .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE_WITH_GAP))
                .columns(
                        col.column("FUSION", "fusion", type.stringType()).setWidth(20),
                        col.column("TYPE", "type", type.stringType()).setWidth(20),
                        col.column("SV TYPE", "svtype", type.stringType()).setWidth(10),
                        col.column("JUNC START", "juncstart", type.stringType()).setWidth(10),
                        col.column("JUNC END", "juncend", type.stringType()).setWidth(10),
                        col.column("SUPPORT", "support", type.stringType()).setWidth(10),
                        col.column("COHORT FREQ", "cohortfreq", type.stringType()).setWidth(10)
                )
                .setDataSource(new JRMapCollectionDataSource(rows));
    }

    private static List<RnaFusion> sort(final List<RnaFusion> fusions)
    {
        return fusions.stream().sorted((f1, f2) ->
        {
            String loc1 = zeroPrefixed(f1.chromosomeStart());
            String loc2 = zeroPrefixed(f2.chromosomeStart());
            if(loc1.equals(loc2))
            {
                return Integer.compare(f1.positionStart(), f2.positionStart());
            }
            return loc1.compareTo(loc2);
        }).collect(Collectors.toList());
    }
}
