package com.hartwig.hmftools.peach.output;

import com.hartwig.hmftools.peach.HaplotypeAnalysis;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import static com.hartwig.hmftools.peach.PeachUtils.TSV_DELIMITER;

public class QcStatusFile
{
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
                .map(e -> toLine(e.getKey(), e.getValue()))
                .forEach(lines::add);
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIMITER)
                .add("gene")
                .add("status")
                .toString();
    }

    private static String toLine(String gene, HaplotypeAnalysis analysis)
    {
        return new StringJoiner(TSV_DELIMITER)
                .add(gene)
                .add(analysis.getAnalysisStatus().toString())
                .toString();
    }
}
