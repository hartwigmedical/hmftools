package com.hartwig.hmftools.common.drivercatalog.dnds;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;

public final class DndsDriverGeneLikelihoodFile {

    private static final String DELIMITER = "\t";
    private static final String HEADER_PREFIX = "gene";

    private DndsDriverGeneLikelihoodFile() {
    }

    @NotNull
    public static List<DndsDriverGeneLikelihood> fromInputStream(@NotNull final InputStream genomeInputStream) {
        return fromLines(new BufferedReader(new InputStreamReader(genomeInputStream)).lines().collect(Collectors.toList()));
    }

    @NotNull
    private static List<DndsDriverGeneLikelihood> fromLines(@NotNull final List<String> lines) {
        return lines.stream()
                .filter(line -> !line.startsWith(HEADER_PREFIX))
                .map(DndsDriverGeneLikelihoodFile::fromString)
                .collect(Collectors.toList());
    }

    @NotNull
    private static DndsDriverImpactLikelihood fromString(int offset, @NotNull final String[] values) {
        int current = offset;
        return ImmutableDndsDriverImpactLikelihood.builder()
                .driversPerSample(Double.parseDouble(values[current++]))
                .passengersPerMutation(Double.parseDouble(values[current]))
                .build();
    }

    @NotNull
    private static ModifiableDndsDriverGeneLikelihood fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        DndsDriverImpactLikelihood missense = fromString(1, values);

        return ModifiableDndsDriverGeneLikelihood.create()
                .setGene(values[0])
                .setMissense(missense)
                .setNonsense(fromString(3, values))
                .setSplice(fromString(5, values))
                .setIndel(fromString(7, values));
    }
}
