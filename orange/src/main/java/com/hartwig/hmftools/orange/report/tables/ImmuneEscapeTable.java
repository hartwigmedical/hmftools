package com.hartwig.hmftools.orange.report.tables;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.type;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.datamodel.immuno.ImmuneEscapeRecord;
import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.jasperreports.engine.data.JRMapCollectionDataSource;

public final class ImmuneEscapeTable
{
    public static JasperReportBuilder build(final String title, final ImmuneEscapeRecord immuneEscape,
            final ReportResources reportResources, final boolean isTumorFail)
    {
        List<Map<String, Object>> rows = new ArrayList<>();
        addRow(rows, "HLA-1 loss-of-function", toYesNoUnavailable(immuneEscape.hasHlaEscape(), isTumorFail));
        addRow(rows, "Antigen presentation pathway inactivation",
                toYesNoUnavailable(immuneEscape.hasAntigenPresentationPathwayEscape(), isTumorFail));
        addRow(rows, "IFN gamma pathway inactivation",
                toYesNoUnavailable(immuneEscape.hasIFNGammaPathwayEscape(), isTumorFail));
        addRow(rows, "(Potential) PD-L1 overexpression",
                toYesNoUnavailable(immuneEscape.hasPDL1OverexpressionEscape(), isTumorFail));
        addRow(rows, "CD58 inactivation", toYesNoUnavailable(immuneEscape.hasCD58InactivationEscape(), isTumorFail));
        addRow(rows, "Epigenetics driven immune escape via SETDB1",
                toYesNoUnavailable(immuneEscape.hasEpigeneticSETDB1Escape(), isTumorFail));

        return report()
                .setColumnTitleStyle(OrangeFonts.TABLE_HEADER_STYLE_UNPADDED)
                .setColumnStyle(OrangeFonts.TABLE_CONTENT_STYLE_UNPADDED)
                .title(cmp.text(title).setStyle(OrangeFonts.TABLE_TITLE_STYLE_WITH_GAP))
                .columns(
                        col.column("ESCAPE MECHANISM", "mechanism", type.stringType()).setWidth(20),
                        col.column("DETECTED?", "detected", type.stringType()).setWidth(10),
                        col.column("", "empty", type.stringType()).setWidth(30)
                )
                .setDataSource(new JRMapCollectionDataSource((java.util.Collection<java.util.Map<String, ?>>) (java.util.Collection<?>) rows));
    }

    private static void addRow(final List<Map<String, Object>> rows, final String mechanism, final String detected)
    {
        Map<String, Object> row = new HashMap<>();
        row.put("mechanism", mechanism);
        row.put("detected", detected);
        row.put("empty", "");
        rows.add(row);
    }

    private static String toYesNoUnavailable(final boolean value, final boolean isTumorFail)
    {
        if(isTumorFail)
        {
            return ReportResources.NOT_AVAILABLE;
        }
        return value ? "Yes" : "No";
    }
}
