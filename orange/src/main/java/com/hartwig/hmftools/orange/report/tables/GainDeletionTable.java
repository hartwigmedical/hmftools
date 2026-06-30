package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.type;

import static com.hartwig.hmftools.datamodel.purple.PurpleDriverType.GERMLINE_AMP;
import static com.hartwig.hmftools.datamodel.purple.PurpleDriverType.GERMLINE_DELETION;
import static com.hartwig.hmftools.orange.report.ReportResources.NOT_AVAILABLE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatFoldChangeField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatPercentileField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatTpmField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.zeroPrefixed;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleGermlineStatus;
import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.column.TextColumnBuilder;
import net.sf.jasperreports.engine.JREmptyDataSource;
import net.sf.jasperreports.engine.data.JRMapCollectionDataSource;

public final class GainDeletionTable
{
    public static JasperReportBuilder build(
            final String title, final List<PurpleGainDeletion> gainsDels, final ReportResources reportResources, boolean hasRna)
    {
        if(gainsDels.isEmpty())
        {
            return report()
                    .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE),
                            cmp.text(ReportResources.NONE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE))
                    .setDataSource(new JREmptyDataSource());
        }

        boolean anyEventsHaveRna = hasRna && gainsDels.stream().anyMatch(x -> x.tpm() != null);

        List<Map<String, ?>> rows = new ArrayList<>();
        for(PurpleGainDeletion gainDel : sort(gainsDels))
        {
            Map<String, Object> row = new LinkedHashMap<>();
            row.put("location", gainDel.chromosome() + gainDel.chromosomeBand());
            row.put("gene", geneDisplay(gainDel));
            row.put("type", typeDisplay(gainDel));
            row.put("range", rangeDisplay(gainDel));
            row.put("mincn", formatSingleDigitDecimal(gainDel.minCopyNumber()));
            row.put("maxcn", formatSingleDigitDecimal(gainDel.maxCopyNumber()));
            row.put("relcn", formatSingleDigitDecimal(gainDel.relativeCopyNumber()));
            row.put("armcn", formatSingleDigitDecimal(gainDel.armCopyNumber()));
            if(anyEventsHaveRna)
            {
                if(gainDel.tpm() != null)
                {
                    row.put("tpm", formatTpmField(gainDel.tpm()));
                    row.put("percentile", formatPercentileField(gainDel.tpmPercentile()));
                    row.put("foldchange", formatFoldChangeField(gainDel.tpmFoldChange()));
                }
                else
                {
                    row.put("tpm", NOT_AVAILABLE);
                    row.put("percentile", NOT_AVAILABLE);
                    row.put("foldchange", NOT_AVAILABLE);
                }
            }
            row.put("driver", gainDel.driver().driverInterpretation().toString());
            rows.add(row);
        }

        List<TextColumnBuilder<String>> columns = new ArrayList<>();
        columns.add(col.column("LOCATION", "location", type.stringType()).setWidth(10));
        columns.add(col.column("GENE", "gene", type.stringType()).setWidth(10));
        columns.add(col.column("TYPE", "type", type.stringType()).setWidth(10));
        columns.add(col.column("RANGE", "range", type.stringType()).setWidth(10));
        columns.add(col.column("MIN CN", "mincn", type.stringType()).setWidth(10));
        columns.add(col.column("MAX CN", "maxcn", type.stringType()).setWidth(10));
        columns.add(col.column("REL CN", "relcn", type.stringType()).setWidth(10));
        columns.add(col.column("ARM CN", "armcn", type.stringType()).setWidth(10));
        if(anyEventsHaveRna)
        {
            columns.add(col.column("TPM", "tpm", type.stringType()).setWidth(10));
            columns.add(col.column("PERCENTILE", "percentile", type.stringType()).setWidth(10));
            columns.add(col.column("FOLD CHANGE", "foldchange", type.stringType()).setWidth(10));
        }
        columns.add(col.column("DRIVER", "driver", type.stringType()).setWidth(10));

        JasperReportBuilder table = report()
                .setColumnTitleStyle(OrangeFonts.TABLE_HEADER_STYLE_UNPADDED)
                .setColumnStyle(OrangeFonts.TABLE_CONTENT_STYLE_UNPADDED)
                .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE_WITH_GAP));

        for(TextColumnBuilder<String> c : columns)
        {
            table.addColumn(c);
        }

        return table.setDataSource(new JRMapCollectionDataSource(rows));
    }

    private static List<PurpleGainDeletion> sort(final List<PurpleGainDeletion> gainsAndDels)
    {
        return gainsAndDels.stream().sorted((gd1, gd2) ->
        {
            int lc = Double.compare(-gd1.driver().driverLikelihood(), -gd2.driver().driverLikelihood());
            if(lc != 0)
            {
                return lc;
            }
            String loc1 = zeroPrefixed(gd1.chromosome() + gd1.chromosomeBand());
            String loc2 = zeroPrefixed(gd2.chromosome() + gd2.chromosomeBand());
            if(!loc1.equals(loc2))
            {
                return loc1.compareTo(loc2);
            }
            return gd1.gene().compareTo(gd2.gene());
        }).collect(Collectors.toList());
    }

    private static String rangeDisplay(final PurpleGainDeletion gainDel)
    {
        if(gainDel.exonStart() == null && gainDel.exonEnd() == null)
        {
            return gainDel.geneRange();
        }
        return format("Exon %d - Exon %d", gainDel.exonStart(), gainDel.exonEnd());
    }

    private static String typeDisplay(final PurpleGainDeletion gainDel)
    {
        if(gainDel.driver().type() == GERMLINE_AMP)
        {
            return "AMP";
        }
        if(gainDel.driver().type() == GERMLINE_DELETION)
        {
            return gainDel.germlineAmpDelFields().germlineStatus() == PurpleGermlineStatus.HOM_DELETION ? "HOM" : "HET";
        }
        return gainDel.driver().type().toString();
    }

    private static String geneDisplay(final PurpleGainDeletion gainDel)
    {
        return gainDel.isCanonical() ? gainDel.gene() : gainDel.gene() + " (alt)";
    }
}
