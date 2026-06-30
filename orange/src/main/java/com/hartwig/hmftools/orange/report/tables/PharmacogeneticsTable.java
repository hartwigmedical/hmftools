package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.common.peach.PeachUtil.UNKNOWN_ALLELE_STRING;
import static com.hartwig.hmftools.common.peach.PeachUtil.convertToZygosityString;

import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.type;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.jasperreports.engine.JREmptyDataSource;
import net.sf.jasperreports.engine.data.JRMapCollectionDataSource;

public final class PharmacogeneticsTable
{
    private static final String PHARMGKB_URL = "https://www.pharmgkb.org";
    private static final String PHARMGKB_NAME = "PHARMGKB";

    public static JasperReportBuilder build(final String title, final Set<PeachGenotype> genotypes,
            final ReportResources reportResources)
    {
        if(genotypes.isEmpty())
        {
            return report()
                    .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE),
                            cmp.text(ReportResources.NONE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE))
                    .setDataSource(new JREmptyDataSource());
        }

        List<Map<String, ?>> rows = new ArrayList<>();
        for(PeachGenotype genotype : sort(genotypes))
        {
            String genotypeString = genotype.allele().equals(UNKNOWN_ALLELE_STRING)
                    ? ""
                    : convertToZygosityString(genotype.alleleCount());
            Map<String, Object> row = new HashMap<>();
            row.put("gene", genotype.gene());
            row.put("haplotype", genotype.allele());
            row.put("genotype", genotypeString);
            row.put("function", genotype.function());
            row.put("linkeddrugs", genotype.linkedDrugs());
            row.put("source", sourceName(genotype.urlPrescriptionInfo()));
            rows.add(row);
        }

        return report()
                .setColumnTitleStyle(OrangeFonts.TABLE_HEADER_STYLE_UNPADDED)
                .setColumnStyle(OrangeFonts.TABLE_CONTENT_STYLE_UNPADDED)
                .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE_WITH_GAP))
                .columns(
                        col.column("GENE", "gene", type.stringType()).setWidth(10),
                        col.column("HAPLOTYPE", "haplotype", type.stringType()).setWidth(10),
                        col.column("GENOTYPE", "genotype", type.stringType()).setWidth(10),
                        col.column("FUNCTION", "function", type.stringType()).setWidth(10),
                        col.column("LINKED DRUGS", "linkeddrugs", type.stringType()).setWidth(20),
                        col.column("SOURCE", "source", type.stringType()).setWidth(10)
                )
                .setDataSource(new JRMapCollectionDataSource(rows));
    }

    private static List<PeachGenotype> sort(final Set<PeachGenotype> genotypes)
    {
        return genotypes.stream().sorted((g1, g2) ->
        {
            if(!g1.gene().equals(g2.gene()))
            {
                return g1.gene().compareTo(g2.gene());
            }
            if(!g1.allele().equals(g2.allele()))
            {
                return g1.allele().compareTo(g2.allele());
            }
            return Integer.compare(g1.alleleCount(), g2.alleleCount());
        }).collect(Collectors.toList());
    }

    private static String sourceName(final String urlPrescriptionInfo)
    {
        String url = extractUrl(urlPrescriptionInfo);
        return url.startsWith(PHARMGKB_URL) ? PHARMGKB_NAME : "";
    }

    private static String extractUrl(final String urlPrescriptionInfo)
    {
        return urlPrescriptionInfo.split(";")[0];
    }
}
