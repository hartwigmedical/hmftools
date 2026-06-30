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

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.datamodel.sigs.SignatureAllocation;
import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.jasperreports.engine.JREmptyDataSource;
import net.sf.jasperreports.engine.data.JRMapCollectionDataSource;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

public final class SignatureAllocationTable
{
    static final String MISALLOC_SIGNATURE = "MISALLOC";

    public static JasperReportBuilder build(
            final String title, final List<SignatureAllocation> signatureAllocations, final ReportResources reportResources)
    {
        if(signatureAllocations.isEmpty())
        {
            return report()
                    .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE),
                            cmp.text(ReportResources.NONE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE))
                    .setDataSource(new JREmptyDataSource());
        }

        List<Map<String, Object>> rows = new ArrayList<>();
        for(SignatureAllocation sig : sort(signatureAllocations))
        {
            if(sig.percent() < 0.01)
            {
                continue;
            }
            Map<String, Object> row = new LinkedHashMap<>();
            row.put("signature", sig.signature());
            row.put("etiology", sig.etiology());
            row.put("allocation", format("%.0f", sig.allocation()));
            row.put("percent", TableCommon.formatPercentage(sig.percent()));
            row.put("empty", "");
            rows.add(row);
        }

        if(rows.isEmpty())
        {
            return report()
                    .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE),
                            cmp.text(ReportResources.NONE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE))
                    .setDataSource(new JREmptyDataSource());
        }

        return report()
                .setColumnTitleStyle(OrangeFonts.TABLE_HEADER_STYLE_UNPADDED)
                .setColumnStyle(OrangeFonts.TABLE_CONTENT_STYLE_UNPADDED)
                .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE_WITH_GAP))
                .columns(
                        col.column("SIGNATURE", "signature", type.stringType()).setWidth(10),
                        col.column("ETIOLOGY", "etiology", type.stringType()).setWidth(20),
                        col.column("ALLOCATION", "allocation", type.stringType()).setWidth(10),
                        col.column("PERCENT", "percent", type.stringType()).setWidth(10),
                        col.column("", "empty", type.stringType()).setWidth(40)
                )
                .setDataSource(new JRMapCollectionDataSource((java.util.Collection<java.util.Map<String, ?>>) (java.util.Collection<?>) rows));
    }

    @VisibleForTesting
    static List<SignatureAllocation> sort(final List<SignatureAllocation> signatureAllocations)
    {
        return signatureAllocations.stream().sorted((a1, a2) ->
        {
            String sig1 = a1.signature().toLowerCase();
            String sig2 = a2.signature().toLowerCase();

            if(sig1.equalsIgnoreCase(MISALLOC_SIGNATURE))
            {
                return 1;
            }
            if(sig2.equalsIgnoreCase(MISALLOC_SIGNATURE))
            {
                return -1;
            }

            if(!sig1.startsWith("sig") || !sig2.startsWith("sig"))
            {
                LOGGER.warn("Signatures do not start with 'sig': {} & {}", sig1, sig2);
                return sig1.compareTo(sig2);
            }

            return Double.compare(a2.allocation(), a1.allocation());
        }).collect(Collectors.toList());
    }
}
