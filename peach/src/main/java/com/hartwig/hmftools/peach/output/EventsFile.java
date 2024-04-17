package com.hartwig.hmftools.peach.output;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.peach.PeachUtils.convertCountToString;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

public class EventsFile
{
    public static void write(@NotNull String filePath, @NotNull Map<String, Integer> eventIdToCount) throws IOException
    {
        Files.write(new File(filePath).toPath(), toLines(eventIdToCount));
    }

    @NotNull
    public static List<String> toLines(@NotNull Map<String, Integer> eventIdToCount)
    {
        List<String> lines = new ArrayList<>();
        lines.add(header());
        eventIdToCount.entrySet().stream().map(e -> toLine(e.getKey(), e.getValue())).sorted().forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(TSV_DELIM).add("event").add("count").toString();
    }

    @NotNull
    private static String toLine(@NotNull String eventId, @Nullable Integer count)
    {
        return new StringJoiner(TSV_DELIM).add(eventId).add(convertCountToString(count)).toString();
    }
}
