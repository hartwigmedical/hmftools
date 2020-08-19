package com.hartwig.hmftools.common.drivercatalog.dnds;

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

final class DndsDriverGeneLikelihoodFile {

    private static final Logger LOGGER = LogManager.getLogger(DndsDriverGeneLikelihoodFile.class);
    private static final String DELIMITER = "\t";
    private static final String HEADER_PREFIX = "gene";

    @NotNull
    static Map<String, DndsDriverImpactLikelihood> fromSingleImpactInputStream(@NotNull final InputStream genomeInputStream) {
        return fromSingleImpactLine(new BufferedReader(new InputStreamReader(genomeInputStream)).lines().collect(Collectors.toList()));
    }

    @NotNull
    static Map<String, ModifiableDndsDriverGeneLikelihood> fromMultiImpactInputStream(@NotNull final InputStream genomeInputStream) {
        return fromMultiImpactLine(new BufferedReader(new InputStreamReader(genomeInputStream)).lines().collect(Collectors.toList()));
    }

    @NotNull
    private static Map<String, ModifiableDndsDriverGeneLikelihood> fromMultiImpactLine(@NotNull final List<String> lines) {
        Map<String, ModifiableDndsDriverGeneLikelihood> result = Maps.newHashMap();
        int i = 0;
        for (String line : lines) {
            i++;
            try {
                if (!line.startsWith(HEADER_PREFIX)) {
                    final ModifiableDndsDriverGeneLikelihood entry = fromString(line);
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
    private static Map<String, DndsDriverImpactLikelihood> fromSingleImpactLine(@NotNull final List<String> lines) {
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
                .dndsLikelihood(0.0)
                .driversPerSample(Double.parseDouble(values[offset++]))
                .passengersPerMutation(Double.parseDouble(values[offset]))
                .build();
    }

    @NotNull
    private static ModifiableDndsDriverGeneLikelihood fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        DndsDriverImpactLikelihood missense = fromString(1, values);

        return ModifiableDndsDriverGeneLikelihood.create()
                .setGene(values[0])
                .setMissense(missense)
                .setMissenseBiallelic(missense)
                .setMissenseNonBiallelic(missense)
                .setNonsense(fromString(3, values))
                .setSplice(fromString(5, values))
                .setIndel(fromString(7, values));
    }
}
