package com.hartwig.hmftools.orange.report.tables;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.orange.report.ReportResources;
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
    public static Table build(@NotNull String title, float width, @NotNull List<PurpleGeneCopyNumber> lohGenes,
            @NotNull ReportResources reportResources) {
        if (lohGenes.isEmpty()) {
            return reportResources.tables().createEmpty(title, width);
        }

        Cells cells = reportResources.cells();
        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 3 },
                new Cell[] { cells.createHeader("Location"), cells.createHeader("Gene"), cells.createHeader("Tumor MACN"),
                        cells.createHeader("Tumor CN"), cells.createHeader(Strings.EMPTY) });

        for (PurpleGeneCopyNumber lohGene : sort(lohGenes)) {
            table.addCell(cells.createContent(lohGene.chromosome() + lohGene.chromosomeBand()));
            table.addCell(cells.createContent(lohGene.geneName()));
            table.addCell(cells.createContent(String.valueOf(round(lohGene.minMinorAlleleCopyNumber()))));
            table.addCell(cells.createContent(String.valueOf(round(lohGene.minCopyNumber()))));
            table.addCell(cells.createContent(Strings.EMPTY));
        }

        return reportResources.tables().createWrapping(table, title);
    }

    private static long round(double copyNumber) {
        return Math.round(Math.max(0, copyNumber));
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
