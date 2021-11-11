package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.linx.ReportableGeneDisruption;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.CellUtil;
import com.hartwig.hmftools.orange.report.util.ChromosomeUtil;
import com.hartwig.hmftools.orange.report.util.TableUtil;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public final class GeneDisruptionTable {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#0.0");

    private GeneDisruptionTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<ReportableGeneDisruption> disruptions) {
        if (disruptions.isEmpty()) {
            return TableUtil.createEmpty(title, width);
        }

        Table table = TableUtil.createContent(width,
                new float[] { 1, 1, 1, 1, 1, 1 },
                new Cell[] { CellUtil.createHeader("Location"), CellUtil.createHeader("Gene"), CellUtil.createHeader("Range"),
                        CellUtil.createHeader("Type"), CellUtil.createHeader("Junction CN"), CellUtil.createHeader("Undisrupted CN") });

        for (ReportableGeneDisruption disruption : sort(disruptions)) {
            table.addCell(CellUtil.createContent(disruption.location()));
            table.addCell(CellUtil.createContent(disruption.gene()));
            table.addCell(CellUtil.createContent(disruption.range()));
            table.addCell(CellUtil.createContent(disruption.type()));
            table.addCell(CellUtil.createContent(SINGLE_DIGIT.format(disruption.junctionCopyNumber())));
            table.addCell(CellUtil.createContent(SINGLE_DIGIT.format(disruption.undisruptedCopyNumber())));
        }

        return TableUtil.createWrapping(table, title);
    }

    @NotNull
    public static List<ReportableGeneDisruption> sort(@NotNull List<ReportableGeneDisruption> disruptions) {
        return disruptions.stream().sorted((disruption1, disruption2) -> {
            String locationAndGene1 = ChromosomeUtil.zeroPrefixed(disruption1.location()) + disruption1.gene();
            String locationAndGene2 = ChromosomeUtil.zeroPrefixed(disruption2.location()) + disruption2.gene();

            if (locationAndGene1.equals(locationAndGene2)) {
                return disruption1.firstAffectedExon() - disruption2.firstAffectedExon();
            } else {
                return locationAndGene1.compareTo(locationAndGene2);
            }
        }).collect(Collectors.toList());
    }
}
