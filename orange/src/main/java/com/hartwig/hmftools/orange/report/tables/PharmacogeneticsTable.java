package com.hartwig.hmftools.orange.report.tables;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class PharmacogeneticsTable {

    private PharmacogeneticsTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<PeachGenotype> genotypes) {
        if (genotypes.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table contentTable = Tables.createContent(width,
                new float[] { 1, 1, 1, 2, 1 },
                new Cell[] { Cells.createHeader("Gene"), Cells.createHeader("Genotype"), Cells.createHeader("Function"),
                        Cells.createHeader("Linked drugs"), Cells.createHeader("Source") });

        for (PeachGenotype genotype : sort(genotypes)) {
            contentTable.addCell(Cells.createContent(genotype.gene()));
            contentTable.addCell(Cells.createContent(genotype.haplotype()));
            contentTable.addCell(Cells.createContent(genotype.function()));
            contentTable.addCell(Cells.createContent(genotype.linkedDrugs()));
            contentTable.addCell(Cells.createUrl(sourceName(genotype.urlPrescriptionInfo()), url(genotype.urlPrescriptionInfo())));
        }

        return Tables.createWrapping(contentTable, title);
    }

    @NotNull
    private static List<PeachGenotype> sort(@NotNull List<PeachGenotype> genotypes) {
        return genotypes.stream().sorted((genotype1, genotype2) -> {
            if (genotype1.gene().equals(genotype2.gene())) {
                return genotype1.haplotype().compareTo(genotype2.haplotype());
            } else {
                return genotype1.gene().compareTo(genotype2.gene());
            }

        }).collect(Collectors.toList());
    }

    @NotNull
    private static String sourceName(@NotNull String urlPrescriptionInfo) {
        String url = extractUrl(urlPrescriptionInfo);
        if (url.startsWith("https://www.pharmgkb.org")) {
            return "PHARMGKB";
        } else {
            return Strings.EMPTY;
        }
    }

    @NotNull
    private static String url(@NotNull String urlPrescriptionInfo) {
        String url = extractUrl(urlPrescriptionInfo);
        if (url.startsWith("https://www.pharmgkb.org")) {
            return url;
        } else {
            return Strings.EMPTY;
        }
    }

    @NotNull
    private static String extractUrl(@NotNull String urlPrescriptionInfo) {
        return urlPrescriptionInfo.split(";")[0];
    }
}
