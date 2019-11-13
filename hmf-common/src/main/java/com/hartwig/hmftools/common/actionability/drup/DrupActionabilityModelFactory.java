package com.hartwig.hmftools.common.actionability.drup;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.io.reader.LineReader;

import org.jetbrains.annotations.NotNull;

public final class DrupActionabilityModelFactory {

    private DrupActionabilityModelFactory() {
    }

    @NotNull
    public static DrupActionabilityModel buildFromCsv(@NotNull String drupGenesCsv) throws IOException {
        List<String> lines = LineReader.build().readLines(new File(drupGenesCsv).toPath(), line -> line.length() > 0);

        Set<String> genes = Sets.newHashSet();

        for (String line : lines) {
            genes.add(line.trim());
        }

        return ImmutableDrupActionabilityModel.of(genes);
    }
}
