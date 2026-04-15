package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.stringArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createStandardTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createEmptyTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.toPercentages;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import java.util.Comparator;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.CommonVcfTags;
import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;

import be.quodlibet.boxable.BaseTable;

import com.hartwig.hmftools.orange.report.DocumentContext;

import java.io.IOException;

import org.jetbrains.annotations.Nullable;

public final class HLAAlleleTable
{
    public static BaseTable build(final DocumentContext docCtx,
            final String title, float width, final List<LilacAllele> alleles, final ReportResources reportResources, boolean hasRna)
            throws IOException
    {
        if(alleles.isEmpty())
        {
            return createEmptyTable(docCtx, title, width, reportResources);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<String> headers = Lists.newArrayList();

        boolean hasWarnings = alleles.stream().anyMatch(x -> !x.qcStatus().equals(CommonVcfTags.PASS_FILTER));

        addEntry(widths, headers, 1, "Allele");
        addEntry(widths, headers, hasWarnings ? 3 : 1, "QC Status");
        addEntry(widths, headers, 1, "Ref Frags");
        addEntry(widths, headers, 1, "Tumor Frags");

        if(hasRna)
        {
            addEntry(widths, headers, 1, "RNA Frags");
        }

        addEntry(widths, headers, 1, "Tumor CN");
        addEntry(widths, headers, 3, "Somatic Mutations");

        BaseTable table = createStandardTable(docCtx, title, width, intToFloatArray(widths), stringArray(headers), reportResources);
        float[] pcts = toPercentages(intToFloatArray(widths));

        for(LilacAllele allele : sort(alleles))
        {
            List<String> rowValues = Lists.newArrayList();
            rowValues.add(allele.allele());
            rowValues.add(allele.qcStatus());
            rowValues.add(fragmentString(allele.refFragments()));
            rowValues.add(fragmentString(allele.tumorFragments()));

            if(hasRna)
            {
                rowValues.add(fragmentString(allele.rnaFragments()));
            }

            rowValues.add(formatSingleDigitDecimal(allele.tumorCopyNumber()));
            rowValues.add(mutationString(allele));
            cells.addRow(table, pcts, rowValues);
        }

        return table;
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
