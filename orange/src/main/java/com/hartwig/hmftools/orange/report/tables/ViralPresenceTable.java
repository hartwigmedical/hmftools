package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatPercentage;
import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;

import java.util.List;
import java.util.Objects;

import com.hartwig.hmftools.datamodel.virus.AnnotatedVirus;
import com.hartwig.hmftools.datamodel.virus.VirusLikelihoodType;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ViralPresenceTable
{
    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<AnnotatedVirus> viruses,
            @NotNull ReportResources reportResources)
    {
        if(viruses.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);
        Table table = Tables.createContent(width,
                new float[] { 4, 3, 1, 1, 2, 2, 2, 2 },
                new Cell[] { cells.createHeader("Virus"), cells.createHeader("QC Status"), cells.createHeader("Type"),
                        cells.createHeader("Int"), cells.createHeader("% Covered"), cells.createHeader("Mean Cov"),
                        cells.createHeader("Exp Clon Cov"), cells.createHeader("Driver") });

        for(AnnotatedVirus virus : viruses)
        {
            table.addCell(cells.createContent(virus.name()));
            table.addCell(cells.createContent(virus.qcStatus().toString()));
            table.addCell(cells.createContent(Objects.requireNonNullElse(virus.interpretation(), Strings.EMPTY)));
            table.addCell(cells.createContent(String.valueOf(virus.integrations())));
            table.addCell(cells.createContent(formatPercentage(virus.percentageCovered(), false)));
            table.addCell(cells.createContent(formatSingleDigitDecimal(virus.meanCoverage())));
            table.addCell(cells.createContent(expectedClonalCoverageField(virus)));
            table.addCell(cells.createContent(display(virus.virusDriverLikelihoodType())));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    private static String display(VirusLikelihoodType virusLikelihoodType)
    {
        switch(virusLikelihoodType)
        {
            case HIGH:
                return "High";
            case LOW:
                return "Low";
            case UNKNOWN:
                return "Unknown";
        }
        throw new IllegalStateException();
    }

    @NotNull
    private static String expectedClonalCoverageField(@NotNull AnnotatedVirus virus)
    {
        Double expectedClonalCoverage = virus.expectedClonalCoverage();
        return expectedClonalCoverage != null ? formatSingleDigitDecimal(expectedClonalCoverage) : Strings.EMPTY;
    }
}