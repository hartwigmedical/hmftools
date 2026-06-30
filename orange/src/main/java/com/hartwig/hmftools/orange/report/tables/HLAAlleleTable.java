package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.type;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.CommonVcfTags;
import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.column.TextColumnBuilder;
import net.sf.jasperreports.engine.JREmptyDataSource;
import net.sf.jasperreports.engine.data.JRMapCollectionDataSource;

import org.jetbrains.annotations.Nullable;

public final class HLAAlleleTable
{
    public static JasperReportBuilder build(final String title, final List<LilacAllele> alleles,
            final ReportResources reportResources, final boolean hasRna)
    {
        if(alleles.isEmpty())
        {
            return report()
                    .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE),
                            cmp.text(ReportResources.NONE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE))
                    .setDataSource(new JREmptyDataSource());
        }

        boolean hasWarnings = alleles.stream().anyMatch(x -> !x.qcStatus().equals(CommonVcfTags.PASS_FILTER));

        List<Map<String, ?>> rows = new ArrayList<>();
        for(LilacAllele allele : sort(alleles))
        {
            Map<String, Object> row = new HashMap<>();
            row.put("allele", allele.allele());
            row.put("qcstatus", allele.qcStatus());
            row.put("reffrags", fragmentString(allele.refFragments()));
            row.put("tumorfrags", fragmentString(allele.tumorFragments()));
            row.put("rnafrags", fragmentString(allele.rnaFragments()));
            row.put("tumorcn", formatSingleDigitDecimal(allele.tumorCopyNumber()));
            row.put("somaticmutations", mutationString(allele));
            rows.add(row);
        }

        List<TextColumnBuilder<String>> columns = new ArrayList<>();
        columns.add(col.column("ALLELE", "allele", type.stringType()).setWidth(10));
        columns.add(col.column("QC STATUS", "qcstatus", type.stringType()).setWidth(hasWarnings ? 30 : 10));
        columns.add(col.column("REF FRAGS", "reffrags", type.stringType()).setWidth(10));
        columns.add(col.column("TUMOR FRAGS", "tumorfrags", type.stringType()).setWidth(10));
        if(hasRna)
        {
            columns.add(col.column("RNA FRAGS", "rnafrags", type.stringType()).setWidth(10));
        }
        columns.add(col.column("TUMOR CN", "tumorcn", type.stringType()).setWidth(10));
        columns.add(col.column("SOMATIC MUTATIONS", "somaticmutations", type.stringType()).setWidth(30));

        JasperReportBuilder table = report()
                .setColumnTitleStyle(OrangeFonts.TABLE_HEADER_STYLE_UNPADDED)
                .setColumnStyle(OrangeFonts.TABLE_CONTENT_STYLE_UNPADDED)
                .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE_WITH_GAP));

        for(TextColumnBuilder<String> column : columns)
        {
            table.addColumn(column);
        }

        return table.setDataSource(new JRMapCollectionDataSource(rows));
    }

    private static List<LilacAllele> sort(final List<LilacAllele> alleles)
    {
        return alleles.stream()
                .sorted(Comparator.comparing(LilacAllele::allele)
                        .thenComparingInt(a -> a.refFragments() != null ? a.refFragments() : 0))
                .collect(Collectors.toList());
    }

    private static String fragmentString(@Nullable final Integer fragments)
    {
        return fragments == null ? ReportResources.NOT_AVAILABLE : String.valueOf(fragments);
    }

    private static String mutationString(final LilacAllele allele)
    {
        StringJoiner joiner = new StringJoiner(", ");
        if(Doubles.positive(allele.somaticMissense()))
        {
            joiner.add(formatSingleDigitDecimal(allele.somaticMissense()) + " missense");
        }
        if(Doubles.positive(allele.somaticNonsenseOrFrameshift()))
        {
            joiner.add(formatSingleDigitDecimal(allele.somaticNonsenseOrFrameshift()) + " nonsense or frameshift");
        }
        if(Doubles.positive(allele.somaticSplice()))
        {
            joiner.add(formatSingleDigitDecimal(allele.somaticSplice()) + " splice");
        }
        if(Doubles.positive(allele.somaticSynonymous()))
        {
            joiner.add(formatSingleDigitDecimal(allele.somaticSynonymous()) + " synonymous");
        }
        if(Doubles.positive(allele.somaticInframeIndel()))
        {
            joiner.add(formatSingleDigitDecimal(allele.somaticInframeIndel()) + " inframe indel");
        }
        String result = joiner.toString();
        return !result.isEmpty() ? result : ReportResources.NONE;
    }
}
