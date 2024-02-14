package com.hartwig.hmftools.peach.output;

import com.hartwig.hmftools.peach.HaplotypeAnalysis;

import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import static com.hartwig.hmftools.peach.PeachUtils.TSV_DELIMITER;

public class EventsPerGeneFile
{
    public static void write(@NotNull String filePath, @NotNull Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis) throws IOException
    {
        Files.write(new File(filePath).toPath(), toLines(geneToHaplotypeAnalysis));
    }

    @NotNull
    public static List<String> toLines(Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis)
    {
        List<String> lines = new ArrayList<>();
        lines.add(header());
        geneToHaplotypeAnalysis.entrySet()
                .stream()
                .map(e -> toLines(e.getKey(), e.getValue()))
                .flatMap(Collection::stream)
                .sorted()
                .forEach(lines::add);
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIMITER).add("gene").add("event").add("count").toString();
    }

    private static List<String> toLines(String gene, HaplotypeAnalysis analysis)
    {
        return analysis.getEventIds().stream().map(e -> toLine(gene, e, analysis.getEventCount(e))).collect(Collectors.toList());
    }

    private static String toLine(String gene, String eventId, int count)
    {
        return new StringJoiner(TSV_DELIMITER).add(gene).add(eventId).add(Integer.toString(count)).toString();
    }
}
