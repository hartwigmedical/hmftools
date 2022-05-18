package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;

import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ViralPresenceTable {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#.#");
    private static final DecimalFormat PERCENTAGE = ReportResources.decimalFormat("#'%'");

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
            table.addCell(Cells.createContent(virus.interpretation() != null ? virus.interpretation() : Strings.EMPTY));
            table.addCell(Cells.createContent(String.valueOf(virus.integrations())));
            table.addCell(Cells.createContent(PERCENTAGE.format(virus.percentageCovered())));
            table.addCell(Cells.createContent(SINGLE_DIGIT.format(virus.meanCoverage())));
            table.addCell(Cells.createContent(expectedClonalCoverageField(virus)));
            table.addCell(Cells.createContent(virus.virusDriverLikelihoodType().display()));
        }

        return Tables.createWrapping(table, title);
    }

    @NotNull
    private static String expectedClonalCoverageField(@NotNull AnnotatedVirus virus) {
        Double expectedClonalCoverage = virus.expectedClonalCoverage();
        return expectedClonalCoverage != null ? SINGLE_DIGIT.format(expectedClonalCoverage) : Strings.EMPTY;
    }
}