package com.hartwig.hmftools.serve.sources.iclusion.filter;

import static com.hartwig.hmftools.serve.actionability.util.ActionableFileFunctions.FIELD_DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public final class IclusionFilterFile {

    private IclusionFilterFile(){
    }

    @NotNull
    public static List<IclusionFilterEntry> read(@NotNull String iClusionFilterTsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(iClusionFilterTsv).toPath());
        // Skip header
        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    private static List<IclusionFilterEntry> fromLines(@NotNull List<String> lines) {
        List<IclusionFilterEntry> filterEntries = Lists.newArrayList();
        for (String line : lines) {
            filterEntries.add(fromLine(line));
        }
        return filterEntries;
    }

    @NotNull
    private static IclusionFilterEntry fromLine(@NotNull String line) {
        String[] values = line.split(FIELD_DELIMITER);

        return ImmutableIclusionFilterEntry.builder()
                .type(IclusionFilterType.valueOf(values[0]))
                .value(values[1])
                .build();
    }
}
