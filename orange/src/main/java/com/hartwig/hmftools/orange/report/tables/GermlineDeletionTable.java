package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public final class GermlineDeletionTable {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#.#");

    private GermlineDeletionTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<GermlineDeletion> deletions) {
        if (deletions.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader("Location"), Cells.createHeader("Gene"), Cells.createHeader("Germline status"),
                        Cells.createHeader("Tumor status"), Cells.createHeader("Germline CN"), Cells.createHeader("Tumor CN") });

        for (GermlineDeletion deletion : sort(deletions)) {
            table.addCell(Cells.createContent(deletion.Chromosome + deletion.ChromosomeBand));
            table.addCell(Cells.createContent(deletion.GeneName));
            table.addCell(Cells.createContent(deletion.NormalStatus.toString()));
            table.addCell(Cells.createContent(deletion.TumorStatus.toString()));
            table.addCell(Cells.createContent(SINGLE_DIGIT.format(deletion.GermlineCopyNumber)));
            table.addCell(Cells.createContent(SINGLE_DIGIT.format(deletion.TumorCopyNumber)));
        }

        return Tables.createWrapping(table, title);
    }

    @NotNull
    private static List<GermlineDeletion> sort(@NotNull List<GermlineDeletion> deletions) {
        return deletions.stream().sorted((deletion1, deletion2) -> {
            String location1 = Chromosomes.zeroPrefixed(deletion1.Chromosome + deletion1.ChromosomeBand);
            String location2 = Chromosomes.zeroPrefixed(deletion2.Chromosome + deletion2.ChromosomeBand);

            if (location1.equals(location2)) {
                return deletion1.GeneName.compareTo(deletion2.GeneName);
            } else {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }
}
