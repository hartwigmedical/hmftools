package com.hartwig.hmftools.orange.report.tables;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.TranscriptRegion;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class LossOfHeterozygosityTable {

    private LossOfHeterozygosityTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<GeneCopyNumber> lohGenes) {
        if (lohGenes.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 3 },
                new Cell[] { Cells.createHeader("Location"), Cells.createHeader("Gene"), Cells.createHeader("Tumor CN"),
                        Cells.createHeader(Strings.EMPTY) });

        for (GeneCopyNumber lohGene : sort(lohGenes)) {
            table.addCell(Cells.createContent(lohGene.chromosome() + lohGene.chromosomeBand()));
            table.addCell(Cells.createContent(lohGene.geneName()));
            table.addCell(Cells.createContent(String.valueOf(Math.round(Math.max(0, lohGene.minCopyNumber())))));
            table.addCell(Cells.createContent(Strings.EMPTY));
        }

        return Tables.createWrapping(table, title);
    }

    @NotNull
    private static List<GeneCopyNumber> sort(@NotNull List<GeneCopyNumber> lohGenes) {
        return lohGenes.stream().sorted(Comparator.comparing(TranscriptRegion::geneName)).collect(Collectors.toList());
    }
}
