package com.hartwig.hmftools.common.dnds;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class DndsDriverGeneLikelihoodFile {
    private static final Logger LOGGER = LogManager.getLogger(DndsDriverGeneLikelihoodFile.class);
    private static final String DELIMITER = "\t";
    private static final String HEADER_PREFIX = "gene";

    @NotNull
    public static Map<String, DndsDriverGeneLikelihood> read(@NotNull final String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    @NotNull
    static Map<String, DndsDriverGeneLikelihood> fromLines(@NotNull List<String> lines) {
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
    private static DndsDriverGeneLikelihood fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        final DndsDriverImpactLikelihood missense = ImmutableDndsDriverImpactLikelihood.builder()
                .dndsLikelihood(Double.valueOf(values[1]))
                .pDriver(Double.valueOf(values[2]))
                .pVariantNonDriverFactor(Double.valueOf(values[3]))
                .build();

        final DndsDriverImpactLikelihood nonsense = ImmutableDndsDriverImpactLikelihood.builder()
                .dndsLikelihood(Double.valueOf(values[4]))
                .pDriver(Double.valueOf(values[5]))
                .pVariantNonDriverFactor(Double.valueOf(values[6]))
                .build();


        final DndsDriverImpactLikelihood splice = ImmutableDndsDriverImpactLikelihood.builder()
                .dndsLikelihood(Double.valueOf(values[7]))
                .pDriver(Double.valueOf(values[8]))
                .pVariantNonDriverFactor(Double.valueOf(values[9]))
                .build();

        final DndsDriverImpactLikelihood indel = ImmutableDndsDriverImpactLikelihood.builder()
                .dndsLikelihood(Double.valueOf(values[10]))
                .pDriver(Double.valueOf(values[11]))
                .pVariantNonDriverFactor(Double.valueOf(values[12]))
                .build();

        final ImmutableDndsDriverGeneLikelihood.Builder builder = ImmutableDndsDriverGeneLikelihood.builder()
                .gene(values[0])
                .missense(missense)
                .nonsense(nonsense)
                .indel(indel)
                .splice(splice);

        return builder.build();
    }

}
