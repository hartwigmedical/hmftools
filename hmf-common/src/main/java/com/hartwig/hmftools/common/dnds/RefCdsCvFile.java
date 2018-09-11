package com.hartwig.hmftools.common.dnds;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.amber.AmberBAFFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class RefCdsCvFile {
    private static final Logger LOGGER = LogManager.getLogger(AmberBAFFile.class);
    private static final String DELIMITER = "\t";
    private static final String HEADER_PREFIX = "gene";

    @NotNull
    public static Map<String, RefCdsCv> read(@NotNull final String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    @NotNull
    static Map<String, RefCdsCv> fromLines(@NotNull List<String> lines) {
        Map<String, RefCdsCv> result = Maps.newHashMap();
        int i = 0;
        for (String line : lines) {
            i++;
            try {
                if (!line.startsWith(HEADER_PREFIX)) {
                    final RefCdsCv entry = fromString(line);
                    result.put(entry.gene(), entry);
                }
            } catch (RuntimeException e) {
                LOGGER.info("Unable to parse line {}: {}", i, line);
                throw e;
            }
        }

        return result;
    }

    @NotNull
    private static RefCdsCv fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        final ImmutableRefCdsCv.Builder builder = ImmutableRefCdsCv.builder()
                .gene(values[0])
                .synonymousN(Integer.valueOf(values[1]))
                .missenseN(Integer.valueOf(values[2]))
                .nonsenseN(Integer.valueOf(values[3]))
                .spliceN(Integer.valueOf(values[4]))
                .indelN(Integer.valueOf(values[5]))
                .missenseW(Double.valueOf(values[6]))
                .nonsenseW(Double.valueOf(values[7]))
                .spliceW(Double.valueOf(values[8]))
                .indelW(Double.valueOf(values[9]));

        return builder.build();
    }

}
