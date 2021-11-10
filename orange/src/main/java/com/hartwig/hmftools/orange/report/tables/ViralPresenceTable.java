package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;

import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.orange.report.ReportResources;
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
            return TableUtil.createEmptyTable(title, width);
        }

        Table table = TableUtil.createReportContentTable(width,
                new float[] { 4, 3, 1, 1, 2, 2, 2},
                new Cell[] { TableUtil.createHeaderCell("Virus"), TableUtil.createHeaderCell("QC Status"),
                        TableUtil.createHeaderCell("Type"), TableUtil.createHeaderCell("Int"),
                        TableUtil.createHeaderCell("% Covered"), TableUtil.createHeaderCell("Mean Cov"),
                        TableUtil.createHeaderCell("Clonal Cov") });

        for (AnnotatedVirus virus : viruses) {
            table.addCell(TableUtil.createContentCell(virus.name()));
            table.addCell(TableUtil.createContentCell(virus.qcStatus().toString()));
            table.addCell(TableUtil.createContentCell(virus.interpretation() != null ? virus.interpretation() : Strings.EMPTY));
            table.addCell(TableUtil.createContentCell(String.valueOf(virus.integrations())));
            table.addCell(TableUtil.createContentCell(PERCENTAGE.format(virus.percentageCovered())));
            table.addCell(TableUtil.createContentCell(SINGLE_DIGIT.format(virus.meanCoverage())));

            Double expectedClonalCoverage = virus.expectedClonalCoverage();
            String clonalCoverageString = expectedClonalCoverage != null ? SINGLE_DIGIT.format(expectedClonalCoverage) : Strings.EMPTY;
            table.addCell(TableUtil.createContentCell(clonalCoverageString));
        }

        return TableUtil.createWrappingReportTable(table, title);
    }
}