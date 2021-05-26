package com.hartwig.hmftools.serve.sources.ckb.filter;

import static com.hartwig.hmftools.serve.actionability.util.ActionableFileFunctions.FIELD_DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public final class CkbFilterFile {

    private CkbFilterFile() {
    }

    @NotNull
    public static Set<CkbFilterEntry> read(@NotNull String ckbFilterTsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(ckbFilterTsv).toPath());
        // Skip header
        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    private static Set<CkbFilterEntry> fromLines(@NotNull List<String> lines) {
        Set<CkbFilterEntry> filterEntries = Sets.newHashSet();
        for (String line : lines) {
            filterEntries.add(fromLine(line));
        }
        return filterEntries;
    }

    @NotNull
    private static CkbFilterEntry fromLine(@NotNull String line) {
        String[] values = line.split(FIELD_DELIMITER);

        return ImmutableCkbFilterEntry.builder()
                .type(CkbFilterType.valueOf(values[0]))
                .value(values[1])
                .build();
    }
}
