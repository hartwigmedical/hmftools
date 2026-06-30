package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.type;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatPercentageField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatTwoDigitDecimal;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.common.AllelicDepth;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.column.TextColumnBuilder;
import net.sf.jasperreports.engine.JREmptyDataSource;
import net.sf.jasperreports.engine.data.JRMapCollectionDataSource;

public final class SomaticVariantTable
{
    public static JasperReportBuilder build(
            final String title, final List<PurpleVariant> variants, final ReportResources reportResources,
            boolean tumorOnly, boolean hasRna)
    {
        if(variants.isEmpty())
        {
            return report()
                    .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE),
                            cmp.text(ReportResources.NONE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE))
                    .setDataSource(new JREmptyDataSource());
        }

        List<Map<String, Object>> rows = new ArrayList<>();
        for(PurpleVariant variant : sort(variants))
        {
            List<PurpleTranscriptImpact> impacts = new ArrayList<>();
            impacts.add(variant.canonicalImpact());
            impacts.addAll(variant.otherImpacts());

            for(PurpleTranscriptImpact impact : impacts)
            {
                Map<String, Object> row = new LinkedHashMap<>();
                row.put("gene", variant.gene());
                row.put("position", locationDisplay(variant));
                row.put("hgvs", hgvsDisplay(impact));
                row.put("af", formatTwoDigitDecimal(variant.tumorDepth().alleleFrequency()));
                row.put("depth", String.valueOf(variant.tumorDepth().totalReadCount()));
                row.put("copies", copyNumberDisplay(variant));
                row.put("hotspot", hotspotDisplay(variant));
                row.put("biallelic", formatPercentageField(variant.biallelicProbability()));
                row.put("clonal", clonalLikelihood(variant));
                if(tumorOnly)
                {
                    row.put("somatic", variant.somaticLikelihood().toString());
                }
                if(hasRna)
                {
                    row.put("rna", rnaDisplay(variant));
                }
                row.put("driver", formatPercentageField(variant.driverLikelihood()));
                rows.add(row);
            }
        }

        List<TextColumnBuilder<String>> columns = new ArrayList<>();
        columns.add(col.column("GENE", "gene", type.stringType()).setWidth(15));
        columns.add(col.column("POSITION", "position", type.stringType()).setWidth(20));
        columns.add(col.column("HGVS", "hgvs", type.stringType()).setWidth(30));
        columns.add(col.column("AF", "af", type.stringType()).setWidth(10));
        columns.add(col.column("DEPTH", "depth", type.stringType()).setWidth(10));
        columns.add(col.column("COPIES", "copies", type.stringType()).setWidth(10));
        columns.add(col.column("HOTSPOT", "hotspot", type.stringType()).setWidth(12));
        columns.add(col.column("BIALL", "biallelic", type.stringType()).setWidth(10));
        columns.add(col.column("CLONAL", "clonal", type.stringType()).setWidth(12));
        if(tumorOnly)
        {
            columns.add(col.column("SOMATIC", "somatic", type.stringType()).setWidth(12));
        }
        if(hasRna)
        {
            columns.add(col.column("RNA", "rna", type.stringType()).setWidth(10));
        }
        columns.add(col.column("DRIVER", "driver", type.stringType()).setWidth(10));

        JasperReportBuilder table = report()
                .setColumnTitleStyle(OrangeFonts.TABLE_HEADER_STYLE_UNPADDED)
                .setColumnStyle(OrangeFonts.TABLE_CONTENT_STYLE_UNPADDED)
                .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE_WITH_GAP));

        for(TextColumnBuilder<String> col : columns)
        {
            table.addColumn(col);
        }

        return table.setDataSource(new JRMapCollectionDataSource((java.util.Collection<java.util.Map<String, ?>>) (java.util.Collection<?>) rows));
    }

    static List<PurpleVariant> sort(final List<PurpleVariant> variants)
    {
        return variants.stream().sorted((v1, v2) ->
        {
            int dc = Double.compare(v2.driverLikelihood(), v1.driverLikelihood());
            if(dc != 0)
            {
                return dc;
            }
            return v1.gene().compareTo(v2.gene());
        }).collect(Collectors.toList());
    }

    static String locationDisplay(final PurpleVariant variant)
    {
        boolean phased = variant.localPhaseSets() != null && !variant.localPhaseSets().isEmpty();
        String loc = format("%s:%d", variant.chromosome(), variant.position());
        return phased ? loc + " *" : loc;
    }

    static String hgvsDisplay(final PurpleTranscriptImpact impact)
    {
        if(impact.hgvsProteinImpact().isEmpty())
        {
            return impact.hgvsCodingImpact();
        }
        return format("%s [%s]", impact.hgvsCodingImpact(), impact.hgvsProteinImpact());
    }

    static String copyNumberDisplay(final PurpleVariant variant)
    {
        return format("%.1f of %.1f", variant.variantCopyNumber(), variant.adjustedCopyNumber());
    }

    static String hotspotDisplay(final PurpleVariant variant)
    {
        switch(variant.hotspot())
        {
            case HOTSPOT:
                return "Yes";
            case NEAR_HOTSPOT:
                return "Near";
            default:
                return "No";
        }
    }

    static String clonalLikelihood(final PurpleVariant variant)
    {
        return formatPercentageField(1 - variant.subclonalLikelihood());
    }

    static String rnaDisplay(final PurpleVariant variant)
    {
        AllelicDepth rnaDepth = variant.rnaDepth();
        if(rnaDepth == null)
        {
            return ReportResources.NOT_AVAILABLE;
        }
        return format("%d/%d", rnaDepth.alleleReadCount(), rnaDepth.totalReadCount());
    }
}
