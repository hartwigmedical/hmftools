package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.createStandardTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.toPercentages;

import java.io.IOException;

import com.hartwig.hmftools.datamodel.immuno.ImmuneEscapeRecord;
import com.hartwig.hmftools.orange.report.DocumentContext;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;

import be.quodlibet.boxable.BaseTable;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ImmuneEscapeTable
{
    public static BaseTable build(final DocumentContext docCtx,
            final String title, float width, final ImmuneEscapeRecord immuneEscape,
            final ReportResources reportResources, boolean isTumorFail) throws IOException
    {
        float[] colWidths = { 2, 1, 3 };
        String[] headerTexts = { "Escape Mechanism", "Detected?", Strings.EMPTY };

        BaseTable table = createStandardTable(docCtx, title, width, colWidths, headerTexts, reportResources);
        float[] pcts = toPercentages(colWidths);
        Cells cells = new Cells(reportResources);

        addEscapeRow(table, cells, pcts, "HLA-1 loss-of-function", immuneEscape.hasHlaEscape(), isTumorFail);
        addEscapeRow(table, cells, pcts, "Antigen presentation pathway inactivation", immuneEscape.hasAntigenPresentationPathwayEscape(), isTumorFail);
        addEscapeRow(table, cells, pcts, "IFN gamma pathway inactivation", immuneEscape.hasIFNGammaPathwayEscape(), isTumorFail);
        addEscapeRow(table, cells, pcts, "(Potential) PD-L1 overexpression", immuneEscape.hasPDL1OverexpressionEscape(), isTumorFail);
        addEscapeRow(table, cells, pcts, "CD58 inactivation", immuneEscape.hasCD58InactivationEscape(), isTumorFail);
        addEscapeRow(table, cells, pcts, "Epigenetics driven immune escape via SETDB1", immuneEscape.hasEpigeneticSETDB1Escape(), isTumorFail);

        return table;
    }

    private static void addEscapeRow(final BaseTable table, final Cells cells, final float[] pcts,
            final String mechanism, boolean detected, boolean isTumorFail)
    {
        java.util.List<String> rowValues = new java.util.ArrayList<>();
        rowValues.add(mechanism);
        rowValues.add(toYesNoUnavailable(detected, isTumorFail));
        rowValues.add(Strings.EMPTY);
        cells.addRow(table, pcts, rowValues);
    }

    @NotNull
    private static String toYesNoUnavailable(boolean value, boolean isTumorFail)
    {
        if(isTumorFail)
        {
            return ReportResources.NOT_AVAILABLE;
        }

        return value ? "Yes" : "No";
    }
}
