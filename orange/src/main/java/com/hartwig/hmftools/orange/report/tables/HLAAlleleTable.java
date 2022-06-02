package com.hartwig.hmftools.orange.report.tables;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.lilac.LilacRecord;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class HLAAlleleTable {

    private HLAAlleleTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<LilacRecord> records) {
        if (records.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 1, 3 },
                new Cell[] { Cells.createHeader("Allele"), Cells.createHeader("Ref Frags"), Cells.createHeader("Tumor Frags"),
                        Cells.createHeader("RNA Frags"), Cells.createHeader("Tumor CN"), Cells.createHeader("Somatic #mutations") });

        for (LilacRecord record : sort(records)) {
            table.addCell(Cells.createContent(Strings.EMPTY));
            table.addCell(Cells.createContent(Strings.EMPTY));
            table.addCell(Cells.createContent(Strings.EMPTY));
            table.addCell(Cells.createContent(Strings.EMPTY));
            table.addCell(Cells.createContent(Strings.EMPTY));
            table.addCell(Cells.createContent(Strings.EMPTY));
        }

        return Tables.createWrapping(table, title);
    }

    @NotNull
    private static List<LilacRecord> sort(@NotNull List<LilacRecord> records) {
        return records.stream()
                .sorted(Comparator.comparing(LilacRecord::allele).thenComparingInt(LilacRecord::refFragments))
                .collect(Collectors.toList());
    }
}
