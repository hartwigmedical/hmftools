package com.hartwig.hmftools.peach.output;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import com.hartwig.hmftools.peach.HaplotypeAnalysis;

import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

public class QcStatusFile
{
    public static void write(@NotNull String filePath, @NotNull Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis) throws IOException
    {
        Files.write(new File(filePath).toPath(), toLines(geneToHaplotypeAnalysis));
    }

    @NotNull
    public static List<String> toLines(@NotNull Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis)
    {
        List<String> lines = new ArrayList<>();
        lines.add(header());
        geneToHaplotypeAnalysis.entrySet()
                .stream()
                .sorted(Map.Entry.comparingByKey())
                .map(e -> toLine(e.getKey(), e.getValue()))
                .forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(TSV_DELIM).add("gene").add("status").toString();
    }

    @NotNull
    private static String toLine(@NotNull String gene, @NotNull HaplotypeAnalysis analysis)
    {
        return new StringJoiner(TSV_DELIM).add(gene).add(analysis.getQcStatus().toString()).toString();
    }
}
