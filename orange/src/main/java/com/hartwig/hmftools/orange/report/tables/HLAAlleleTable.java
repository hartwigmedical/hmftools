package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_VARIANT;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_RNA;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.floatArray;

import java.util.Comparator;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class HLAAlleleTable
{
    public static Table build(
            final String title, float width, final List<LilacAllele> alleles, final ReportResources reportResources, boolean hasRna)
    {
        if(alleles.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        addEntry(cells, widths, cellEntries, 1, "Allele");
        addEntry(cells, widths, cellEntries, 1, "QC Status");
        addEntry(cells, widths, cellEntries, 1, "Ref Frags");
        addEntry(cells, widths, cellEntries, 1, "Tumor Frags");

        if(hasRna)
        {
            addEntry(cells, widths, cellEntries, 1, "RNA Frags");
        }

        addEntry(cells, widths, cellEntries, 1, "Tumor CN");
        addEntry(cells, widths, cellEntries, 3, "Somatic Mutations");

        Table table = Tables.createContent(width, floatArray(widths), cellArray(cellEntries));

        for(LilacAllele allele : sort(alleles))
        {
            table.addCell(cells.createContent(allele.allele()));
            table.addCell(cells.createContent(allele.qcStatus()));
            table.addCell(cells.createContent(fragmentString(allele.refFragments())));
            table.addCell(cells.createContent(fragmentString(allele.tumorFragments())));

            if(hasRna)
                table.addCell(cells.createContent(fragmentString(allele.rnaFragments())));

            table.addCell(cells.createContent(formatSingleDigitDecimal(allele.tumorCopyNumber())));
            table.addCell(cells.createContent(mutationString(allele)));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    private static List<LilacAllele> sort(final List<LilacAllele> alleles)
    {
        return alleles.stream()
                .sorted(Comparator.comparing(LilacAllele::allele)
                        .thenComparingInt(allele -> allele.refFragments() != null ? allele.refFragments() : 0))
                .collect(Collectors.toList());
    }

    private static String fragmentString(@Nullable final Integer fragments)
    {
        if(fragments == null)
        {
            return ReportResources.NOT_AVAILABLE;
        }
        else
        {
            return String.valueOf(fragments);
        }
    }

    private static String mutationString(final LilacAllele allele)
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
        return !result.isEmpty() ? result : ReportResources.NONE;
    }
}
