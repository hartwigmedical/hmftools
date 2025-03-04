package com.hartwig.hmftools.peach.output;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.peach.PeachUtils.GERMLINE_TOTAL_COPY_NUMBER;

import com.hartwig.hmftools.peach.HaplotypeAnalysis;
import com.hartwig.hmftools.peach.haplotype.HaplotypeCombination;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class AllHaplotypeCombinationsFile
{
    public final static String COMBINATION_SEPARATOR = ";";

    public static void write(final String filePath, final Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis) throws IOException
    {
        Files.write(new File(filePath).toPath(), toLines(geneToHaplotypeAnalysis));
    }

    public static List<String> toLines(final Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis)
    {
        List<String> lines = new ArrayList<>();
        lines.add(header());
        geneToHaplotypeAnalysis.entrySet()
                .stream()
                .sorted(Map.Entry.comparingByKey())
                .map(e -> toLines(e.getKey(), e.getValue()))
                .flatMap(Collection::stream)
                .forEach(lines::add);
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIM).add("gene").add("combination").add("count").add("nonWildTypeCount").toString();
    }

    private static List<String> toLines(final String gene, final HaplotypeAnalysis analysis)
    {
        String wildTypeName = analysis.getWildTypeHaplotypeName();
        Comparator<HaplotypeCombination> combinationComparator =
                Comparator.comparingInt((HaplotypeCombination c) -> Math.abs(c.getHaplotypeCount() - GERMLINE_TOTAL_COPY_NUMBER))
                        .thenComparing((HaplotypeCombination c) -> c.getHaplotypeCountWithout(wildTypeName))
                        .thenComparing(AllHaplotypeCombinationsFile::getHaplotypeCombinationString);
        return analysis.getHaplotypeCombinations()
                .stream()
                .sorted(combinationComparator)
                .map(c -> toLine(gene, c, analysis.getWildTypeHaplotypeName()))
                .collect(Collectors.toList());
    }

    private static String toLine(final String gene, final HaplotypeCombination combination, final String wildTypeHaplotypeName)
    {
        return new StringJoiner(TSV_DELIM).add(gene)
                .add(getHaplotypeCombinationString(combination))
                .add(Integer.toString(combination.getHaplotypeCount()))
                .add(Integer.toString(combination.getHaplotypeCountWithout(wildTypeHaplotypeName)))
                .toString();
    }

    private static String getHaplotypeCombinationString(final HaplotypeCombination combination)
    {
        return combination.getHaplotypeNameToCount()
                .entrySet()
                .stream()
                .map(e -> String.format("(%s,%s)", e.getKey(), e.getValue()))
                .sorted()
                .collect(Collectors.joining(COMBINATION_SEPARATOR));
    }
}
