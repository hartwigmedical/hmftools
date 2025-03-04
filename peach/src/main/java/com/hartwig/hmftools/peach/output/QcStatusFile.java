package com.hartwig.hmftools.peach.output;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import com.hartwig.hmftools.peach.HaplotypeAnalysis;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

public class QcStatusFile
{
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
                .map(e -> toLine(e.getKey(), e.getValue()))
                .forEach(lines::add);
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIM).add("gene").add("status").toString();
    }

    private static String toLine(final String gene, final HaplotypeAnalysis analysis)
    {
        return new StringJoiner(TSV_DELIM).add(gene).add(analysis.getQcStatus().toString()).toString();
    }
}
