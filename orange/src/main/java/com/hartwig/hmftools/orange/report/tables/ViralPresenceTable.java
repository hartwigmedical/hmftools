package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;

import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.CellUtil;
import com.hartwig.hmftools.orange.report.util.TableUtil;
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
            return TableUtil.createEmpty(title, width);
        }

        Table table = TableUtil.createContent(width,
                new float[] { 4, 3, 1, 1, 2, 2, 2, 2 },
                new Cell[] { CellUtil.createHeader("Virus"), CellUtil.createHeader("QC Status"), CellUtil.createHeader("Type"),
                        CellUtil.createHeader("Int"), CellUtil.createHeader("% Covered"), CellUtil.createHeader("Mean Cov"),
                        CellUtil.createHeader("Exp Clon Cov"), CellUtil.createHeader("Driver") });

        for (AnnotatedVirus virus : viruses) {
            table.addCell(CellUtil.createContent(virus.name()));
            table.addCell(CellUtil.createContent(virus.qcStatus().toString()));
            table.addCell(CellUtil.createContent(virus.interpretation() != null ? virus.interpretation() : Strings.EMPTY));
            table.addCell(CellUtil.createContent(String.valueOf(virus.integrations())));
            table.addCell(CellUtil.createContent(PERCENTAGE.format(virus.percentageCovered())));
            table.addCell(CellUtil.createContent(SINGLE_DIGIT.format(virus.meanCoverage())));
            table.addCell(CellUtil.createContent(expectedClonalCoverageField(virus)));
            table.addCell(CellUtil.createContent(virus.virusDriverLikelihoodType().display()));
        }

        return TableUtil.createWrapping(table, title);
    }

    @NotNull
    private static String expectedClonalCoverageField(@NotNull AnnotatedVirus virus) {
        Double expectedClonalCoverage = virus.expectedClonalCoverage();
        return expectedClonalCoverage != null ? SINGLE_DIGIT.format(expectedClonalCoverage) : Strings.EMPTY;
    }
}