package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class BreakendTable {

    private static final Logger LOGGER = LogManager.getLogger(BreakendTable.class);

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#0.0");

    private BreakendTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<LinxBreakend> breakends) {
        if (breakends.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader("Location"), Cells.createHeader("Gene"), Cells.createHeader("Range"),
                        Cells.createHeader("Type"), Cells.createHeader("Cluster ID"), Cells.createHeader("Junction CN"),
                        Cells.createHeader("Undisrupted CN") });

        for (LinxBreakend breakend : sort(breakends)) {
            table.addCell(Cells.createContent(location(breakend)));
            table.addCell(Cells.createContent(displayGene(breakend)));
            table.addCell(Cells.createContent(rangeField(breakend)));
            table.addCell(Cells.createContent(breakend.type().toString()));
            // TODO Populate cluster ID
            table.addCell(Cells.createContent("cl"));
            table.addCell(Cells.createContent(SINGLE_DIGIT.format(breakend.junctionCopyNumber())));
            table.addCell(Cells.createContent(SINGLE_DIGIT.format(breakend.undisruptedCopyNumber())));
        }

        return Tables.createWrapping(table, title);
    }

    @NotNull
    public static List<LinxBreakend> sort(@NotNull List<LinxBreakend> breakends) {
        return breakends.stream().sorted((breakend1, breakend2) -> {
            String locationAndGene1 = Chromosomes.zeroPrefixed(location(breakend1)) + breakend1.gene();
            String locationAndGene2 = Chromosomes.zeroPrefixed(location(breakend2)) + breakend2.gene();

            if (locationAndGene1.equals(locationAndGene2)) {
                return breakend1.exonUp() - breakend2.exonUp();
            } else {
                return locationAndGene1.compareTo(locationAndGene2);
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    private static String displayGene(@NotNull LinxBreakend breakend) {
        String addon = Strings.EMPTY;
        if (!breakend.canonical()) {
            addon = " (alt)";
        }
        return breakend.gene() + addon;
    }

    @NotNull
    private static String rangeField(@NotNull LinxBreakend breakend) {
        String exonRange = null;
        if (breakend.exonUp() > 0) {
            if (breakend.exonUp() == breakend.exonDown()) {
                exonRange = String.format("Exon %d", breakend.exonUp());
            } else if (breakend.exonDown() - breakend.exonUp() == 1) {
                exonRange = String.format("Intron %d", breakend.exonUp());
            }
        } else if (breakend.exonUp() == 0 && (breakend.exonDown() == 1 || breakend.exonDown() == 2)) {
            exonRange = "Promoter Region";
        }

        if (exonRange == null) {
            LOGGER.warn("Could not format range for breakend: {}", breakend);
            return Strings.EMPTY;
        }

        return exonRange + " " + breakend.geneOrientation();
    }

    @NotNull
    private static String location(@NotNull LinxBreakend breakend) {
        return breakend.chromosome() + breakend.chrBand();
    }
}
