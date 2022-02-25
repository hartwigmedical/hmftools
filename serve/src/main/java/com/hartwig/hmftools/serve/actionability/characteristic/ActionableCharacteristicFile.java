package com.hartwig.hmftools.serve.actionability.characteristic;

import static com.hartwig.hmftools.serve.actionability.util.ActionableFileFunctions.FIELD_DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.serve.actionability.util.ActionableFileFunctions;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicAnnotation;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicsAtLeast;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ActionableCharacteristicFile {

    private static final String ACTIONABLE_CHARACTERISTIC_TSV = "ActionableCharacteristics.tsv";

    private ActionableCharacteristicFile() {
    }

    @NotNull
    public static String actionableCharacteristicTsvPath(@NotNull String serveActionabilityDir,
            @NotNull RefGenomeVersion refGenomeVersion) {
        return refGenomeVersion.addVersionToFilePath(serveActionabilityDir + File.separator + ACTIONABLE_CHARACTERISTIC_TSV);
    }

    public static void write(@NotNull String actionableCharacteristicTsv,
            @NotNull Iterable<ActionableCharacteristic> actionableCharacteristics) throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(toLines(actionableCharacteristics));
        Files.write(new File(actionableCharacteristicTsv).toPath(), lines);
    }

    @NotNull
    public static List<ActionableCharacteristic> read(@NotNull String actionableCharacteristicTsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(actionableCharacteristicTsv).toPath());
        // Skip header
        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    private static String header() {
        return new StringJoiner(FIELD_DELIMITER).add("name").add("atLeast").add("cutOff").add(ActionableFileFunctions.header()).toString();
    }

    @NotNull
    @VisibleForTesting
    static List<ActionableCharacteristic> fromLines(@NotNull List<String> lines) {
        List<ActionableCharacteristic> actionableCharacteristics = Lists.newArrayList();
        for (String line : lines) {
            actionableCharacteristics.add(fromLine(line));
        }
        return actionableCharacteristics;
    }

    @NotNull
    private static ActionableCharacteristic fromLine(@NotNull String line) {
        String[] values = line.split(FIELD_DELIMITER);

        return ImmutableActionableCharacteristic.builder()
                .from(ActionableFileFunctions.fromLine(values, 3))
                .name(TumorCharacteristicAnnotation.valueOf(values[0]))
                .atLeast(!values[1].equals(Strings.EMPTY) ? TumorCharacteristicsAtLeast.valueOf(values[1]) : null)
                .cutOff(!values[2].equals(Strings.EMPTY) ? Double.valueOf(values[2]) : null)
                .build();
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull Iterable<ActionableCharacteristic> actionableCharacteristics) {
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
                .add(characteristic.atLeast() != null ? characteristic.atLeast().toString() : Strings.EMPTY)
                .add(characteristic.cutOff() != null ? Double.toString(characteristic.cutOff()) : Strings.EMPTY)
                .add(ActionableFileFunctions.toLine(characteristic))
                .toString();
    }
}
