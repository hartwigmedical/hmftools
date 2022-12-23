package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.Comparator;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public final class HLAAlleleTable {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#.#");

    private HLAAlleleTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<LilacAllele> alleles) {
        if (alleles.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 1, 3 },
                new Cell[] { Cells.createHeader("Allele"), Cells.createHeader("Ref Frags"), Cells.createHeader("Tumor Frags"),
                        Cells.createHeader("RNA Frags"), Cells.createHeader("Tumor CN"), Cells.createHeader("Somatic #mutations") });

        for (LilacAllele allele : sort(alleles)) {
            table.addCell(Cells.createContent(allele.allele()));
            table.addCell(Cells.createContent(String.valueOf(allele.refFragments())));
            table.addCell(Cells.createContent(String.valueOf(allele.tumorFragments())));
            table.addCell(Cells.createContent(String.valueOf(allele.rnaFragments())));
            table.addCell(Cells.createContent(SINGLE_DIGIT.format(allele.tumorCopyNumber())));
            table.addCell(Cells.createContent(mutationString(allele)));
        }

        return Tables.createWrapping(table, title);
    }

    @NotNull
    private static List<LilacAllele> sort(@NotNull List<LilacAllele> alleles) {
        return alleles.stream()
                .sorted(Comparator.comparing(LilacAllele::allele).thenComparingInt(LilacAllele::refFragments))
                .collect(Collectors.toList());
    }

    @NotNull
    private static String mutationString(@NotNull LilacAllele allele) {
        StringJoiner joiner = new StringJoiner(", ");
        if (Doubles.positive(allele.somaticMissense())) {
            joiner.add(SINGLE_DIGIT.format(allele.somaticMissense()) + " missense");
        }

        if (Doubles.positive(allele.somaticNonsenseOrFrameshift())) {
            joiner.add(SINGLE_DIGIT.format(allele.somaticNonsenseOrFrameshift()) + " nonsense or frameshift");
        }

        if (Doubles.positive(allele.somaticSplice())) {
            joiner.add(SINGLE_DIGIT.format(allele.somaticSplice()) + " splice");
        }

        if (Doubles.positive(allele.somaticSynonymous())) {
            joiner.add(SINGLE_DIGIT.format(allele.somaticSynonymous()) + " synonymous");
        }

        if (Doubles.positive(allele.somaticInframeIndel())) {
            joiner.add(SINGLE_DIGIT.format(allele.somaticInframeIndel()) + " inframe indel");
        }

        String result = joiner.toString();
        return !result.isEmpty() ? result : "None";
    }
}
