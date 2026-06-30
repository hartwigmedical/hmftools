package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.type;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.zeroPrefixed;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxGeneOrientation;
import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.jasperreports.engine.JREmptyDataSource;
import net.sf.jasperreports.engine.data.JRMapCollectionDataSource;

import org.jetbrains.annotations.Nullable;

public final class DisruptionTable
{
    public static JasperReportBuilder build(
            final String title, final List<LinxBreakend> breakends, final ReportResources reportResources)
    {
        if(breakends.isEmpty())
        {
            return report()
                    .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE),
                            cmp.text(ReportResources.NONE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE))
                    .setDataSource(new JREmptyDataSource());
        }

        List<LinxBreakend> sorted = sort(breakends);
        List<Map<String, ?>> rows = new ArrayList<>();

        for(int i = 0; i < sorted.size(); i++)
        {
            LinxBreakend breakend = sorted.get(i);
            LinxBreakend other = null;
            if(i + 1 < sorted.size() && sorted.get(i + 1).svId() == breakend.svId())
            {
                other = sorted.get(i + 1);
                i++;
            }

            Map<String, Object> row = new LinkedHashMap<>();
            row.put("gene", breakend.gene());
            row.put("position", locationDisplay(breakend, other));
            row.put("zygosity", breakend.undisruptedCopyNumber() < 0.5 ? "HOM" : "HET");
            row.put("context", contextDisplay(breakend, other));
            row.put("type", breakend.type().toString());
            row.put("jcn", formatSingleDigitDecimal(breakend.junctionCopyNumber()));
            row.put("undisruptedcn", formatSingleDigitDecimal(breakend.undisruptedCopyNumber()));
            row.put("driver", DriverInterpretation.interpret(breakend.driverLikelihood()).toString());
            rows.add(row);
        }

        return report()
                .setColumnTitleStyle(OrangeFonts.TABLE_HEADER_STYLE_UNPADDED)
                .setColumnStyle(OrangeFonts.TABLE_CONTENT_STYLE_UNPADDED)
                .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE_WITH_GAP))
                .columns(
                        col.column("GENE", "gene", type.stringType()).setWidth(10),
                        col.column("POSITION", "position", type.stringType()).setWidth(30),
                        col.column("ZYGOSITY", "zygosity", type.stringType()).setWidth(10),
                        col.column("CONTEXT", "context", type.stringType()).setWidth(20),
                        col.column("TYPE", "type", type.stringType()).setWidth(10),
                        col.column("JCN", "jcn", type.stringType()).setWidth(10),
                        col.column("UNDISRUPTED CN", "undisruptedcn", type.stringType()).setWidth(20),
                        col.column("DRIVER", "driver", type.stringType()).setWidth(10)
                )
                .setDataSource(new JRMapCollectionDataSource(rows));
    }

    private static String locationDisplay(final LinxBreakend breakend)
    {
        return format("%s:%d", breakend.chromosome(), breakend.position());
    }

    private static String locationDisplay(final LinxBreakend lower, @Nullable final LinxBreakend upper)
    {
        if(upper == null)
        {
            return locationDisplay(lower);
        }
        return format("%s - %s", locationDisplay(lower), locationDisplay(upper));
    }

    private static String contextDisplay(final LinxBreakend lower, @Nullable final LinxBreakend upper)
    {
        if(upper == null)
        {
            return contextStr(lower, true);
        }
        return format("%s - %s", contextStr(lower, false), contextStr(upper, false));
    }

    private static String contextStr(final LinxBreakend breakend, boolean includeOrientation)
    {
        String exonRange = null;
        if(breakend.exonUp() > 0)
        {
            if(breakend.exonUp() == breakend.exonDown())
            {
                exonRange = format("Exon %d", breakend.exonUp());
            }
            else if(breakend.exonDown() - breakend.exonUp() == 1)
            {
                exonRange = format("Intron %d", breakend.exonUp());
            }
        }
        else if(breakend.exonUp() == 0 && (breakend.exonDown() == 1 || breakend.exonDown() == 2))
        {
            exonRange = "Promoter Region";
        }

        return includeOrientation ? format("%s %s", exonRange, orientationStr(breakend.geneOrientation())) : exonRange;
    }

    private static String orientationStr(final LinxGeneOrientation orientation)
    {
        switch(orientation)
        {
            case UPSTREAM:
                return "Upstream";
            case DOWNSTREAM:
                return "Downstream";
        }
        return null;
    }

    private static List<LinxBreakend> sort(final List<LinxBreakend> breakends)
    {
        return breakends.stream().sorted((b1, b2) ->
        {
            if(b1.svId() == b2.svId())
            {
                return Integer.compare(b1.id(), b2.id());
            }
            String loc1 = zeroPrefixed(b1.chromosome()) + b1.chromosomeBand();
            String loc2 = zeroPrefixed(b2.chromosome()) + b2.chromosomeBand();
            int lc = loc1.compareTo(loc2);
            if(lc != 0)
            {
                return lc;
            }
            int gc = b1.gene().compareTo(b2.gene());
            if(gc != 0)
            {
                return gc;
            }
            return b1.exonUp() - b2.exonUp();
        }).collect(Collectors.toList());
    }
}
