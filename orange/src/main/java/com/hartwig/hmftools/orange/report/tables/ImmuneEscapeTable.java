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
                new float[] { 1, 1, 5 },
                new Cell[] { cells.createHeader("Escape Mechanism"), cells.createHeader("Present?"), cells.createHeader(Strings.EMPTY) });

        table.addCell(cells.createKey("HLA-1 loss-of-function"));
        table.addCell(cells.createValue(toYesNoUnavailable(immuneEscape.hasHlaEscapePresent(), isTumorFail)));
        table.addCell(cells.createContent(Strings.EMPTY));

        table.addCell(cells.createKey("Antigen presentation pathway inactivation"));
        table.addCell(cells.createValue(toYesNoUnavailable(immuneEscape.hasAntigenPresentationPathwayEscape(), isTumorFail)));
        table.addCell(cells.createContent(Strings.EMPTY));

        table.addCell(cells.createKey("IFN gamma pathway inactivation"));
        table.addCell(cells.createValue(toYesNoUnavailable(immuneEscape.hasIFNGammaPathwayEscape(), isTumorFail)));
        table.addCell(cells.createContent(Strings.EMPTY));

        table.addCell(cells.createKey("(Potential) PD-L1 overexpression"));
        table.addCell(cells.createValue(toYesNoUnavailable(immuneEscape.hasPDL1OverexpressionEscape(), isTumorFail)));
        table.addCell(cells.createContent(Strings.EMPTY));

        table.addCell(cells.createKey("CD58 inactivation"));
        table.addCell(cells.createValue(toYesNoUnavailable(immuneEscape.hasCD58InactivationEscape(), isTumorFail)));
        table.addCell(cells.createContent(Strings.EMPTY));

        table.addCell(cells.createKey("Epigenetic driven immune escape via SETDB1"));
        table.addCell(cells.createValue(toYesNoUnavailable(immuneEscape.hasEpigeneticSETDB1Escape(), isTumorFail)));
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
