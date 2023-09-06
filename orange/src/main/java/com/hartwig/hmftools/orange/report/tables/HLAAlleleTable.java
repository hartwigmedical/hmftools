package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;

import java.util.Comparator;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public final class HLAAlleleTable
{
    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<LilacAllele> alleles,
            @NotNull ReportResources reportResources)
    {
        if(alleles.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);
        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 1, 3 },
                new Cell[] { cells.createHeader("Allele"), cells.createHeader("Ref Frags"), cells.createHeader("Tumor Frags"),
                        cells.createHeader("RNA Frags"), cells.createHeader("Tumor CN"), cells.createHeader("Somatic #mutations") });

        for(LilacAllele allele : sort(alleles))
        {
            table.addCell(cells.createContent(allele.allele()));
            table.addCell(cells.createContent(String.valueOf(allele.refFragments())));
            table.addCell(cells.createContent(String.valueOf(allele.tumorFragments())));
            table.addCell(cells.createContent(String.valueOf(allele.rnaFragments())));
            table.addCell(cells.createContent(formatSingleDigitDecimal(allele.tumorCopyNumber())));
            table.addCell(cells.createContent(mutationString(allele)));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    @NotNull
    private static List<LilacAllele> sort(@NotNull List<LilacAllele> alleles)
    {
        return alleles.stream()
                .sorted(Comparator.comparing(LilacAllele::allele).thenComparingInt(LilacAllele::refFragments))
                .collect(Collectors.toList());
    }

    @NotNull
    private static String mutationString(@NotNull LilacAllele allele)
    {
        StringJoiner joiner = new StringJoiner(", ");
        if(Doubles.positive(allele.somaticMissense()))
        {
            joiner.add(formatSingleDigitDecimal(allele.somaticMissense()) + " missense");
        }

        if(Doubles.positive(allele.somaticNonsenseOrFrameshift()))
        {
            joiner.add(formatSingleDigitDecimal(allele.somaticNonsenseOrFrameshift()) + " nonsense or frameshift");
        }

        if(Doubles.positive(allele.somaticSplice()))
        {
            joiner.add(formatSingleDigitDecimal(allele.somaticSplice()) + " splice");
        }

        if(Doubles.positive(allele.somaticSynonymous()))
        {
            joiner.add(formatSingleDigitDecimal(allele.somaticSynonymous()) + " synonymous");
        }

        if(Doubles.positive(allele.somaticInframeIndel()))
        {
            joiner.add(formatSingleDigitDecimal(allele.somaticInframeIndel()) + " inframe indel");
        }

        String result = joiner.toString();
        return !result.isEmpty() ? result : "None";
    }
}
