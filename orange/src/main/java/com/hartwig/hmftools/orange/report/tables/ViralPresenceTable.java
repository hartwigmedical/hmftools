package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatPercentage;
import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;

import java.util.List;

import com.hartwig.hmftools.datamodel.virus.AnnotatedVirus;
import com.hartwig.hmftools.datamodel.virus.VirusLikelihoodType;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ViralPresenceTable {

    private ViralPresenceTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<AnnotatedVirus> viruses) {
        if (viruses.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 4, 3, 1, 1, 2, 2, 2, 2 },
                new Cell[] { Cells.createHeader("Virus"), Cells.createHeader("QC Status"), Cells.createHeader("Type"),
                        Cells.createHeader("Int"), Cells.createHeader("% Covered"), Cells.createHeader("Mean Cov"),
                        Cells.createHeader("Exp Clon Cov"), Cells.createHeader("Driver") });

        for (AnnotatedVirus virus : viruses) {
            table.addCell(Cells.createContent(virus.name()));
            table.addCell(Cells.createContent(virus.qcStatus().toString()));
            table.addCell(Cells.createContent(virus.interpretation() != null ? virus.interpretation().name() : Strings.EMPTY));
            table.addCell(Cells.createContent(String.valueOf(virus.integrations())));
            table.addCell(Cells.createContent(formatPercentage(virus.percentageCovered(), false)));
            table.addCell(Cells.createContent(formatSingleDigitDecimal(virus.meanCoverage())));
            table.addCell(Cells.createContent(expectedClonalCoverageField(virus)));
            table.addCell(Cells.createContent(display(virus.virusDriverLikelihoodType())));
        }

        return Tables.createWrapping(table, title);
    }

    private static String display(VirusLikelihoodType virusLikelihoodType) {
        switch (virusLikelihoodType) {
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
    private static String expectedClonalCoverageField(@NotNull AnnotatedVirus virus) {
        Double expectedClonalCoverage = virus.expectedClonalCoverage();
        return expectedClonalCoverage != null ? formatSingleDigitDecimal(expectedClonalCoverage) : Strings.EMPTY;
    }
}