package com.hartwig.hmftools.peach.output;

import com.hartwig.hmftools.peach.HaplotypeAnalysis;
import com.hartwig.hmftools.peach.HaplotypeCombination;
import org.jetbrains.annotations.NotNull;

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

import static com.hartwig.hmftools.peach.PeachUtils.GERMLINE_TOTAL_COPY_NUMBER;
import static com.hartwig.hmftools.peach.PeachUtils.TSV_DELIMITER;

public class BestHaplotypeCombinationsFile
{
    private static final String UNKNOWN_ALLELE_STRING = "UNKNOWN";

    public static void write(@NotNull String filePath, @NotNull Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis) throws IOException
    {
        Files.write(new File(filePath).toPath(), toLines(geneToHaplotypeAnalysis));
    }

    @NotNull
    private static List<String> toLines(Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis)
    {
        List<String> lines = new ArrayList<>();
        lines.add(header());
        geneToHaplotypeAnalysis.entrySet().stream()
                .sorted(Map.Entry.comparingByKey())
                .map(e -> toLines(e.getKey(), e.getValue()))
                .flatMap(Collection::stream)
                .forEach(lines::add);
        return lines;
    }

    private static String header() {
        return new StringJoiner(TSV_DELIMITER)
                .add("gene")
                .add("allele")
                .add("count")
                .toString();
    }

    private static List<String> toLines(String gene, HaplotypeAnalysis analysis)
    {
        StringJoiner joiner = new StringJoiner(TSV_DELIMITER);
        if (analysis.hasBestHaplotypeCombination())
        {
            return analysis.getBestHaplotypeCombination().getHaplotypeNameToCount().entrySet().stream()
                    .sorted(Map.Entry.comparingByKey())
                    .map(e -> toLine(gene, e.getKey(), e.getValue()))
                    .collect(Collectors.toList());
        }
        else
        {
            String unknownAlleleString = joiner
                    .add(gene)
                    .add(UNKNOWN_ALLELE_STRING)
                    .add(Integer.toString(GERMLINE_TOTAL_COPY_NUMBER))
                    .toString();
            return List.of(unknownAlleleString);
        }
    }

    private static String toLine(String gene, String haplotypeName, int count)
    {
        return new StringJoiner(TSV_DELIMITER)
                .add(gene)
                .add(haplotypeName)
                .add(Integer.toString(count))
                .toString();
    }
}
