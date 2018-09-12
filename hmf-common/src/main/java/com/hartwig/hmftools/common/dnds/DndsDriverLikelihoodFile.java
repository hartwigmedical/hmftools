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

class DndsDriverLikelihoodFile {
    private static final Logger LOGGER = LogManager.getLogger(AmberBAFFile.class);
    private static final String DELIMITER = "\t";
    private static final String HEADER_PREFIX = "gene";

    @NotNull
    public static Map<String, DndsDriverLikelihood> read(@NotNull final String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    @NotNull
    static Map<String, DndsDriverLikelihood> fromLines(@NotNull List<String> lines) {
        Map<String, DndsDriverLikelihood> result = Maps.newHashMap();
        int i = 0;
        for (String line : lines) {
            i++;
            try {
                if (!line.startsWith(HEADER_PREFIX)) {
                    final DndsDriverLikelihood entry = fromString(line);
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
    private static DndsDriverLikelihood fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        final ImmutableDndsDriverLikelihood.Builder builder = ImmutableDndsDriverLikelihood.builder()
                .gene(values[0])
                .indel(Double.valueOf(values[1]))
                .missense(Double.valueOf(values[2]))
                .nonsense(Double.valueOf(values[3]))
                .splice(Double.valueOf(values[4]));

        return builder.build();
    }

}
