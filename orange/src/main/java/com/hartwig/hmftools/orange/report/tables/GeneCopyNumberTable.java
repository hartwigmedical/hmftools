package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.isofox.IsofoxData;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.orange.report.ReportResources;
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

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#.#");
    private static final DecimalFormat PERCENTAGE = ReportResources.decimalFormat("#'%'");

    private static final Logger LOGGER = LogManager.getLogger(GeneCopyNumberTable.class);

    private GeneCopyNumberTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<ReportableGainLoss> gainLosses,
            @Nullable IsofoxData isofox) {
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
                table.addCell(Cells.createContent(SINGLE_DIGIT.format(expression.tpm())));
                table.addCell(Cells.createContent(PERCENTAGE.format(expression.percentileCancer() * 100)));
                table.addCell(Cells.createContent(formatFoldChange(expression.tpm() / expression.medianTpmCancer())));
                table.addCell(Cells.createContent(PERCENTAGE.format(expression.percentileCohort() * 100)));
                table.addCell(Cells.createContent(formatFoldChange(expression.tpm() / expression.medianTpmCohort())));
            } else {
                table.addCell(Cells.createContent(Strings.EMPTY));
                table.addCell(Cells.createContent(Strings.EMPTY));
                table.addCell(Cells.createContent(Strings.EMPTY));
                table.addCell(Cells.createContent(Strings.EMPTY));
                table.addCell(Cells.createContent(Strings.EMPTY));
            }
        }

        return Tables.createWrapping(table, title);
    }

    @NotNull
    private static String formatFoldChange(double foldChange) {
        return foldChange > 1000 ? ">1000" : SINGLE_DIGIT.format(foldChange);
    }

    @Nullable
    private static GeneExpression findExpressionForGene(@Nullable IsofoxData isofox, @NotNull String gene) {
        if (isofox == null) {
            return null;
        }

        for (GeneExpression expression : isofox.geneExpressions()) {
            if (expression.geneName().equals(gene)) {
                return expression;
            }
        }

        LOGGER.warn("Could not find expression data for gene '{}'", gene);
        return null;
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
