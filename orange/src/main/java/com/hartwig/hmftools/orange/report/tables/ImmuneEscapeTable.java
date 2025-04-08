package com.hartwig.hmftools.orange.report.tables;

import com.hartwig.hmftools.datamodel.immuno.ImmuneEscapeRecord;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ImmuneEscapeTable
{
    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull ImmuneEscapeRecord immuneEscape,
            @NotNull ReportResources reportResources, boolean isTumorFail)
    {
        Cells cells = new Cells(reportResources);
        Table table = Tables.createContent(width,
                new float[] { 2, 1, 3 },
                new Cell[] { cells.createHeader("Escape Mechanism"), cells.createHeader("Detected?"), cells.createHeader(Strings.EMPTY) });

        table.addCell(cells.createContent("HLA-1 loss-of-function"));
        table.addCell(cells.createContent(toYesNoUnavailable(immuneEscape.hasHlaEscape(), isTumorFail)));
        table.addCell(cells.createContent(Strings.EMPTY));

        table.addCell(cells.createContent("Antigen presentation pathway inactivation"));
        table.addCell(cells.createContent(toYesNoUnavailable(immuneEscape.hasAntigenPresentationPathwayEscape(), isTumorFail)));
        table.addCell(cells.createContent(Strings.EMPTY));

        table.addCell(cells.createContent("IFN gamma pathway inactivation"));
        table.addCell(cells.createContent(toYesNoUnavailable(immuneEscape.hasIFNGammaPathwayEscape(), isTumorFail)));
        table.addCell(cells.createContent(Strings.EMPTY));

        table.addCell(cells.createContent("(Potential) PD-L1 overexpression"));
        table.addCell(cells.createContent(toYesNoUnavailable(immuneEscape.hasPDL1OverexpressionEscape(), isTumorFail)));
        table.addCell(cells.createContent(Strings.EMPTY));

        table.addCell(cells.createContent("CD58 inactivation"));
        table.addCell(cells.createContent(toYesNoUnavailable(immuneEscape.hasCD58InactivationEscape(), isTumorFail)));
        table.addCell(cells.createContent(Strings.EMPTY));

        table.addCell(cells.createContent("Epigenetics driven immune escape via SETDB1"));
        table.addCell(cells.createContent(toYesNoUnavailable(immuneEscape.hasEpigeneticSETDB1Escape(), isTumorFail)));
        table.addCell(cells.createContent(Strings.EMPTY));

        return new Tables(reportResources).createWrapping(table, title);
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
