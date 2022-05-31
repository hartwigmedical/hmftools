package com.hartwig.hmftools.orange.report.tables;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxInterpretedData;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.interpretation.Expressions;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class GeneCopyNumberTable {

    private static final Logger LOGGER = LogManager.getLogger(GeneCopyNumberTable.class);

    private GeneCopyNumberTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<ReportableGainLoss> gainLosses,
            @Nullable IsofoxInterpretedData isofox) {
        if (gainLosses.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader("Chr"), Cells.createHeader("Region"), Cells.createHeader("Gene"),
                        Cells.createHeader("Type"), Cells.createHeader("CN"), Cells.createHeader("TPM"), Cells.createHeader("Perc (Type))"),
                        Cells.createHeader("FC (Type)"), Cells.createHeader("Perc (DB)"), Cells.createHeader("FC (DB)") });

        for (ReportableGainLoss gainLoss : sort(gainLosses)) {
            table.addCell(Cells.createContent(gainLoss.chromosome()));
            table.addCell(Cells.createContent(gainLoss.chromosomeBand()));
            table.addCell(Cells.createContent(displayGene(gainLoss)));
            table.addCell(Cells.createContent(gainLoss.interpretation().display()));
            table.addCell(Cells.createContent(String.valueOf(gainLoss.minCopies())));

            GeneExpression expression = findExpressionForGene(isofox, gainLoss.gene());
            if (expression != null) {
                table.addCell(Cells.createContent(Expressions.tpm(expression)));
                table.addCell(Cells.createContent(Expressions.percentileType(expression)));
                table.addCell(Cells.createContent(Expressions.foldChangeType(expression)));
                table.addCell(Cells.createContent(Expressions.percentileDatabase(expression)));
                table.addCell(Cells.createContent(Expressions.foldChangeDatabase(expression)));
            } else {
                table.addCell(Cells.createContent(ReportResources.NOT_AVAILABLE));
                table.addCell(Cells.createContent(ReportResources.NOT_AVAILABLE));
                table.addCell(Cells.createContent(ReportResources.NOT_AVAILABLE));
                table.addCell(Cells.createContent(ReportResources.NOT_AVAILABLE));
                table.addCell(Cells.createContent(ReportResources.NOT_AVAILABLE));
            }
        }

        return Tables.createWrapping(table, title);
    }

    @Nullable
    private static GeneExpression findExpressionForGene(@Nullable IsofoxInterpretedData isofox, @NotNull String geneToFind) {
        if (isofox == null) {
            return null;
        }

        return Expressions.findByGene(isofox.allGeneExpressions(), geneToFind);
    }

    @NotNull
    private static List<ReportableGainLoss> sort(@NotNull List<ReportableGainLoss> reportableGainsAndLosses) {
        return reportableGainsAndLosses.stream().sorted((gainLoss1, gainLoss2) -> {
            String location1 = Chromosomes.zeroPrefixed(gainLoss1.chromosome() + gainLoss1.chromosomeBand());
            String location2 = Chromosomes.zeroPrefixed(gainLoss2.chromosome() + gainLoss2.chromosomeBand());

            if (location1.equals(location2)) {
                return gainLoss1.gene().compareTo(gainLoss2.gene());
            } else {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    private static String displayGene(@NotNull ReportableGainLoss gainLoss) {
        String addon = Strings.EMPTY;
        if (!gainLoss.isCanonical()) {
            addon = " (alt)";
        }
        return gainLoss.gene() + addon;
    }
}
