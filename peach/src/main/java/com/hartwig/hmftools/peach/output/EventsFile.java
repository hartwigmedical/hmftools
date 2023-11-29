package com.hartwig.hmftools.peach.output;

import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import static com.hartwig.hmftools.peach.PeachUtils.TSV_DELIMITER;

public class EventsFile
{
    public static void write(@NotNull String filePath, @NotNull Map<String, Integer> eventIdToCount) throws IOException
    {
        Files.write(new File(filePath).toPath(), toLines(eventIdToCount));
    }

    @NotNull
    public static List<String> toLines(Map<String, Integer> eventIdToCount)
    {
        List<String> lines = new ArrayList<>();
        lines.add(header());
        eventIdToCount.entrySet().stream().map(e -> toLine(e.getKey(), e.getValue())).sorted().forEach(lines::add);
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIMITER)
                .add("event")
                .add("count")
                .toString();
    }

    private static String toLine(String eventId, int count)
    {
        return new StringJoiner(TSV_DELIMITER)
                .add(eventId)
                .add(Integer.toString(count))
                .toString();
    }
}
