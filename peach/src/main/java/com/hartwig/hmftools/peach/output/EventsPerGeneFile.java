package com.hartwig.hmftools.peach.output;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.peach.PeachUtils.convertCountToString;

import com.hartwig.hmftools.peach.HaplotypeAnalysis;

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
                .map(e -> toLines(e.getKey(), e.getValue()))
                .flatMap(Collection::stream)
                .sorted()
                .forEach(lines::add);
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIM).add("gene").add("event").add("count").toString();
    }

    private static List<String> toLines(final String gene, final HaplotypeAnalysis analysis)
    {
        return analysis.getEventIds().stream().map(e -> toLine(gene, e, analysis.getEventCount(e))).collect(Collectors.toList());
    }

    private static String toLine(final String gene, final String eventId, @Nullable Integer count)
    {
        return new StringJoiner(TSV_DELIM).add(gene).add(eventId).add(convertCountToString(count)).toString();
    }
}
