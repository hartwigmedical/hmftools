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

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.jasperreports.engine.JREmptyDataSource;
import net.sf.jasperreports.engine.data.JRMapCollectionDataSource;

public final class GermlineVariantTable
{
    public static JasperReportBuilder build(
            final String title, final List<PurpleVariant> variants, final ReportResources reportResources)
    {
        if(variants.isEmpty())
        {
            return report()
                    .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE),
                            cmp.text(ReportResources.NONE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE))
                    .setDataSource(new JREmptyDataSource());
        }

        List<Map<String, ?>> rows = new ArrayList<>();
        for(PurpleVariant variant : SomaticVariantTable.sort(variants))
        {
            List<PurpleTranscriptImpact> impacts = new ArrayList<>();
            impacts.add(variant.canonicalImpact());
            impacts.addAll(variant.otherImpacts());

            for(PurpleTranscriptImpact impact : impacts)
            {
                Map<String, Object> row = new LinkedHashMap<>();
                row.put("gene", variant.gene());
                row.put("variant", SomaticVariantTable.locationDisplay(variant));
                row.put("hgvs", SomaticVariantTable.hgvsDisplay(impact));
                row.put("zygosity", simplifiedDisplay(variant.genotypeStatus()));
                row.put("af", TableCommon.formatTwoDigitDecimal(variant.tumorDepth().alleleFrequency()));
                row.put("depth", String.valueOf(variant.tumorDepth().totalReadCount()));
                row.put("copies", SomaticVariantTable.copyNumberDisplay(variant));
                row.put("hotspot", SomaticVariantTable.hotspotDisplay(variant));
                row.put("biallelic", variant.biallelic() ? "Yes" : "No");
                row.put("gnomad", format("%4.2e", variant.gnomadFrequency()));
                row.put("clinvar", variant.clinvarPathogenicity());
                row.put("driver", DriverInterpretation.interpret(variant.driverLikelihood()).toString());
                rows.add(row);
            }
        }

        return report()
                .setColumnTitleStyle(OrangeFonts.TABLE_HEADER_STYLE_UNPADDED)
                .setColumnStyle(OrangeFonts.TABLE_CONTENT_STYLE_UNPADDED)
                .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE_WITH_GAP))
                .columns(
                        col.column("GENE", "gene", type.stringType()).setWidth(15),
                        col.column("VARIANT", "variant", type.stringType()).setWidth(20),
                        col.column("HGVS", "hgvs", type.stringType()).setWidth(20),
                        col.column("ZYGOSITY", "zygosity", type.stringType()).setWidth(15),
                        col.column("AF", "af", type.stringType()).setWidth(8),
                        col.column("DEPTH", "depth", type.stringType()).setWidth(10),
                        col.column("COPIES", "copies", type.stringType()).setWidth(10),
                        col.column("HOTSPOT", "hotspot", type.stringType()).setWidth(15),
                        col.column("BIALLELIC", "biallelic", type.stringType()).setWidth(10),
                        col.column("GNOMAD", "gnomad", type.stringType()).setWidth(15),
                        col.column("CLINVAR", "clinvar", type.stringType()).setWidth(20),
                        col.column("DRIVER", "driver", type.stringType()).setWidth(12)
                )
                .setDataSource(new JRMapCollectionDataSource(rows));
    }

    private static String simplifiedDisplay(final PurpleGenotypeStatus genotypeStatus)
    {
        switch(genotypeStatus)
        {
            case HOM_REF:
            case HOM_ALT:
                return TableCommon.VALUE_HOM;
            case HET:
                return TableCommon.VALUE_HET;
            case UNKNOWN:
                return "UNKNOWN";
        }
        throw new IllegalStateException("Unknown PurpleGenotypeStatus: " + genotypeStatus);
    }
}
