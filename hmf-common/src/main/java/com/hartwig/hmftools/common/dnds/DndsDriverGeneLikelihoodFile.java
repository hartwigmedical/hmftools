package com.hartwig.hmftools.common.dnds;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class DndsDriverGeneLikelihoodFile {
    private static final Logger LOGGER = LogManager.getLogger(DndsDriverGeneLikelihoodFile.class);
    private static final String DELIMITER = "\t";
    private static final String HEADER_PREFIX = "gene";

    @NotNull
    public static Map<String, DndsDriverImpactLikelihood> fromOncoInputStream(@NotNull final InputStream genomeInputStream) {
        return fromOncoLines(new BufferedReader(new InputStreamReader(genomeInputStream)).lines().collect(Collectors.toList()));
    }

    @NotNull
    public static Map<String, DndsDriverGeneLikelihood> fromTsgInputStream(@NotNull final InputStream genomeInputStream) {
        return fromTsgLines(new BufferedReader(new InputStreamReader(genomeInputStream)).lines().collect(Collectors.toList()));
    }

    @NotNull
    private static Map<String, DndsDriverGeneLikelihood> fromTsgLines(@NotNull final List<String> lines) {
        Map<String, DndsDriverGeneLikelihood> result = Maps.newHashMap();
        int i = 0;
        for (String line : lines) {
            i++;
            try {
                if (!line.startsWith(HEADER_PREFIX)) {
                    final DndsDriverGeneLikelihood entry = fromString(line);
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
    private static Map<String, DndsDriverImpactLikelihood> fromOncoLines(@NotNull final List<String> lines) {
        Map<String, DndsDriverImpactLikelihood> result = Maps.newHashMap();
        int i = 0;
        for (String line : lines) {
            i++;
            try {
                if (!line.startsWith(HEADER_PREFIX)) {
                    String[] values = line.split(DELIMITER);
                    final DndsDriverImpactLikelihood entry = fromString(1, values);
                    result.put(values[0], entry);
                }
            } catch (RuntimeException e) {
                LOGGER.info("Unable to parse line {}: {}", i, line);
                throw e;
            }
        }

        return result;
    }

    @NotNull
    private static DndsDriverImpactLikelihood fromString(int offset, @NotNull final String[] values) {
        return ImmutableDndsDriverImpactLikelihood.builder()
                .dndsLikelihood(Double.valueOf(values[offset++]))
                .pDriver(Double.valueOf(values[offset++]))
                .pVariantNonDriverFactor(Double.valueOf(values[offset]))
                .build();
    }

    @NotNull
    private static DndsDriverGeneLikelihood fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        final ImmutableDndsDriverGeneLikelihood.Builder builder = ImmutableDndsDriverGeneLikelihood.builder()
                .gene(values[0])
                .missense(fromString(1, values))
                .nonsense(fromString(4, values))
                .splice(fromString(7, values))
                .indel(fromString(10, values));

        return builder.build();
    }
}
