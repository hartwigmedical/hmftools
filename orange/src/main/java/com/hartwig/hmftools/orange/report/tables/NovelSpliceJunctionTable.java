package com.hartwig.hmftools.orange.report.tables;

import static java.lang.Math.round;
import static java.lang.String.format;

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

import com.hartwig.hmftools.datamodel.isofox.AltSpliceJunctionType;
import com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction;
import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.jasperreports.engine.JREmptyDataSource;
import net.sf.jasperreports.engine.data.JRMapCollectionDataSource;

public final class NovelSpliceJunctionTable
{
    public static JasperReportBuilder build(final String title, final List<NovelSpliceJunction> junctions,
            final ReportResources reportResources)
    {
        if(junctions.isEmpty())
        {
            return report()
                    .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE),
                            cmp.text(ReportResources.NONE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE))
                    .setDataSource(new JREmptyDataSource());
        }

        List<Map<String, ?>> rows = new ArrayList<>();
        for(NovelSpliceJunction junction : sort(junctions))
        {
            int fragments = junction.fragmentCount();
            int averageDepth = (int) round((junction.depthStart() + junction.depthEnd()) * 0.5);

            Map<String, Object> row = new HashMap<>();
            row.put("gene", junction.gene());
            row.put("junctions", junctionsDisplay(junction));
            row.put("type", junction.type().toString());
            row.put("juncstart", String.valueOf(junction.regionStart()));
            row.put("juncend", String.valueOf(junction.regionEnd()));
            row.put("support", formatSupportField(fragments, averageDepth));
            row.put("cohortfreq", String.valueOf(junction.cohortFrequency()));
            rows.add(row);
        }

        return report()
                .setColumnTitleStyle(OrangeFonts.TABLE_HEADER_STYLE_UNPADDED)
                .setColumnStyle(OrangeFonts.TABLE_CONTENT_STYLE_UNPADDED)
                .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE_WITH_GAP))
                .columns(
                        col.column("GENE", "gene", type.stringType()).setWidth(10),
                        col.column("JUNCTIONS", "junctions", type.stringType()).setWidth(30),
                        col.column("TYPE", "type", type.stringType()).setWidth(20),
                        col.column("JUNC START", "juncstart", type.stringType()).setWidth(10),
                        col.column("JUNC END", "juncend", type.stringType()).setWidth(10),
                        col.column("SUPPORT", "support", type.stringType()).setWidth(10),
                        col.column("COHORT FREQ", "cohortfreq", type.stringType()).setWidth(10)
                )
                .setDataSource(new JRMapCollectionDataSource(rows));
    }

    private static String junctionsDisplay(final NovelSpliceJunction junction)
    {
        boolean dupType = junction.type() == AltSpliceJunctionType.CIRCULAR;
        int positionStart = dupType ? junction.junctionEnd() : junction.junctionStart();
        int positionEnd = dupType ? junction.junctionStart() : junction.junctionEnd();
        int exonStart = dupType ? junction.exonEnd() : junction.exonStart();
        int exonEnd = dupType ? junction.exonStart() : junction.exonEnd();
        return format("Exon %d (%s:%d) - Exon %d (%s:%d)",
                exonStart, junction.chromosome(), positionStart, exonEnd, junction.chromosome(), positionEnd);
    }

    private static List<NovelSpliceJunction> sort(final List<NovelSpliceJunction> junctions)
    {
        return junctions.stream().sorted((j1, j2) ->
        {
            String loc1 = zeroPrefixed(j1.chromosome());
            String loc2 = zeroPrefixed(j2.chromosome());
            if(loc1.equals(loc2))
            {
                return Integer.compare(j1.junctionStart(), j2.junctionStart());
            }
            return loc1.compareTo(loc2);
        }).collect(Collectors.toList());
    }
}
