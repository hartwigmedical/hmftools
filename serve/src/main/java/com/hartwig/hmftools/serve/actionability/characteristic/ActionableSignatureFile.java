package com.hartwig.hmftools.serve.actionability.characteristic;

import static com.hartwig.hmftools.serve.actionability.util.ActionableFileFunctions.FIELD_DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.serve.actionability.util.ActionableFileFunctions;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristic;

import org.jetbrains.annotations.NotNull;

// TODO: Can be removed once PROTECT reads actionable characteristics
public final class ActionableSignatureFile {

    private static final String ACTIONABLE_SIGNATURE_TSV = "ActionableSignatures.tsv";

    private static final Set<TumorCharacteristic> SUPPORTED_CHARACTERISTICS =
            Sets.newHashSet(TumorCharacteristic.HIGH_TUMOR_MUTATIONAL_LOAD,
                    TumorCharacteristic.HOMOLOGOUS_RECOMBINATION_DEFICIENT,
                    TumorCharacteristic.MICROSATELLITE_UNSTABLE);

    private ActionableSignatureFile() {
    }

    @NotNull
    public static String actionableSignatureTsvPath(@NotNull String serveActionabilityDir, @NotNull RefGenomeVersion refGenomeVersion) {
        return refGenomeVersion.addVersionToFilePath(serveActionabilityDir + File.separator + ACTIONABLE_SIGNATURE_TSV);
    }

    public static void write(@NotNull String actionableSignatureTsv, @NotNull Iterable<ActionableCharacteristic> actionableCharacteristics)
            throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(toLines(filter(actionableCharacteristics)));
        Files.write(new File(actionableSignatureTsv).toPath(), lines);
    }

    @NotNull
    private static Iterable<ActionableCharacteristic> filter(@NotNull Iterable<ActionableCharacteristic> actionableCharacteristics) {
        List<ActionableCharacteristic> filtered = Lists.newArrayList();
        for (ActionableCharacteristic characteristic : actionableCharacteristics) {
            if (SUPPORTED_CHARACTERISTICS.contains(characteristic.name())) {
                filtered.add(characteristic);
            }
        }
        return filtered;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(FIELD_DELIMITER).add("name").add(ActionableFileFunctions.header()).toString();
    }

    @NotNull
    private static List<String> toLines(@NotNull Iterable<ActionableCharacteristic> actionableCharacteristics) {
        List<String> lines = Lists.newArrayList();
        for (ActionableCharacteristic actionableCharacteristic : sort(actionableCharacteristics)) {
            lines.add(toLine(actionableCharacteristic));
        }
        return lines;
    }

    @NotNull
    private static List<ActionableCharacteristic> sort(@NotNull Iterable<ActionableCharacteristic> actionableCharacteristics) {
        // Need to make a copy since the input may be immutable and cannot be sorted!
        List<ActionableCharacteristic> sorted = Lists.newArrayList(actionableCharacteristics);
        sorted.sort(new ActionableCharacteristicComparator());

        return sorted;
    }

    @NotNull
    private static String toLine(@NotNull ActionableCharacteristic characteristic) {
        return new StringJoiner(FIELD_DELIMITER).add(characteristic.name().toString())
                .add(ActionableFileFunctions.toLine(characteristic))
                .toString();
    }
}