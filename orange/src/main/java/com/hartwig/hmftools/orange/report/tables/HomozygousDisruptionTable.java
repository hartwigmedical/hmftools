package com.hartwig.hmftools.orange.report.tables;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.linx.HomozygousDisruption;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class HomozygousDisruptionTable {

    private HomozygousDisruptionTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<HomozygousDisruption> homozygousDisruptions) {
        if (homozygousDisruptions.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 1, 1, 4 },
                new Cell[] { Cells.createHeader("Location"), Cells.createHeader("Gene"), Cells.createHeader(Strings.EMPTY) });

        for (HomozygousDisruption homozygousDisruption : sort(homozygousDisruptions)) {
            table.addCell(Cells.createContent(homozygousDisruption.chromosome() + homozygousDisruption.chromosomeBand()));
            table.addCell(Cells.createContent(gene(homozygousDisruption)));
            table.addCell(Cells.createContent(Strings.EMPTY));
        }

        return Tables.createWrapping(table, title);
    }

    @NotNull
    private static List<HomozygousDisruption> sort(@NotNull List<HomozygousDisruption> homozygousDisruptions) {
        return homozygousDisruptions.stream().sorted((disruption1, disruption2) -> {
            String location1 = Chromosomes.zeroPrefixed(disruption1.chromosome() + disruption1.chromosomeBand());
            String location2 = Chromosomes.zeroPrefixed(disruption2.chromosome() + disruption2.chromosomeBand());

            if (location1.equals(location2)) {
                return disruption1.gene().compareTo(disruption2.gene());
            } else {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    private static String gene(@NotNull HomozygousDisruption homozygousDisruption) {
        String addon = Strings.EMPTY;
        if (!homozygousDisruption.isCanonical()) {
            addon = " (alt)";
        }
        return homozygousDisruption.gene() + addon;
    }
}
