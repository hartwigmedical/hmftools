package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
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
    public static Table build(@NotNull String title, float width, @NotNull List<PurpleGeneCopyNumber> lohGenes) {
        if (lohGenes.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 3 },
                new Cell[] { Cells.createHeader("Location"), Cells.createHeader("Gene"), Cells.createHeader("Tumor MACN"),
                        Cells.createHeader("Tumor CN"), Cells.createHeader(Strings.EMPTY) });

        for (PurpleGeneCopyNumber lohGene : sort(lohGenes)) {
            table.addCell(Cells.createContent(lohGene.chromosome() + lohGene.chromosomeBand()));
            table.addCell(Cells.createContent(lohGene.geneName()));
            table.addCell(Cells.createContent(formatSingleDigitDecimal(Math.max(0, lohGene.minMinorAlleleCopyNumber()))));
            table.addCell(Cells.createContent(formatSingleDigitDecimal(Math.max(0, lohGene.minCopyNumber()))));
            table.addCell(Cells.createContent(Strings.EMPTY));
        }

        return Tables.createWrapping(table, title);
    }

    @NotNull
    private static List<PurpleGeneCopyNumber> sort(@NotNull List<PurpleGeneCopyNumber> lohGenes) {
        return lohGenes.stream().sorted((geneCopyNumber1, geneCopyNumber2) -> {
            String location1 = Chromosomes.zeroPrefixed(geneCopyNumber1.chromosome() + geneCopyNumber1.chromosomeBand());
            String location2 = Chromosomes.zeroPrefixed(geneCopyNumber2.chromosome() + geneCopyNumber2.chromosomeBand());

            if (location1.equals(location2)) {
                return geneCopyNumber1.geneName().compareTo(geneCopyNumber2.geneName());
            } else {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }
}
