package com.hartwig.hmftools.orange.report.tables;

import java.util.List;

import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.AnnotatedVirusV1;
import com.hartwig.hmftools.orange.report.util.TableUtil;
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
            return TableUtil.createEmptyTable(title, width);
        }

        Table table = TableUtil.createReportContentTable(width,
                new float[] { 1, 1, 1, 1 },
                new Cell[] { TableUtil.createHeaderCell("Virus"), TableUtil.createHeaderCell("QC Status"),
                        TableUtil.createHeaderCell("Interpretation"), TableUtil.createHeaderCell("Integrations") });

        for (AnnotatedVirus virus : viruses) {
            table.addCell(TableUtil.createContentCell(virus.name()));
            table.addCell(TableUtil.createContentCell(virus.qcStatus().toString()));
            table.addCell(TableUtil.createContentCell(virus.interpretation() != null ? virus.interpretation().toString() : Strings.EMPTY));
            table.addCell(TableUtil.createContentCell(Integer.toString(virus.integrations())));
        }

        return TableUtil.createWrappingReportTable(table, title);
    }
}