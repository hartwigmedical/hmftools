package com.hartwig.hmftools.peach.output;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.peach.PeachUtils.convertCountToString;

import com.hartwig.hmftools.peach.HaplotypeAnalysis;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class EventsPerGeneFile
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
                .map(e -> toLines(e.getKey(), e.getValue()))
                .flatMap(Collection::stream)
                .sorted()
                .forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(TSV_DELIM).add("gene").add("event").add("count").toString();
    }

    @NotNull
    private static List<String> toLines(@NotNull String gene, @NotNull HaplotypeAnalysis analysis)
    {
        return analysis.getEventIds().stream().map(e -> toLine(gene, e, analysis.getEventCount(e))).collect(Collectors.toList());
    }

    @NotNull
    private static String toLine(@NotNull String gene, @NotNull String eventId, @Nullable Integer count)
    {
        return new StringJoiner(TSV_DELIM).add(gene).add(eventId).add(convertCountToString(count)).toString();
    }
}
