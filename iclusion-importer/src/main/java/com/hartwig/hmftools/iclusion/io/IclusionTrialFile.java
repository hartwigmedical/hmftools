package com.hartwig.hmftools.iclusion.io;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.iclusion.data.IclusionMutation;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.iclusion.data.IclusionTumorLocation;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionMutation;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionTrial;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionTumorLocation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class IclusionTrialFile {

    private static final Logger LOGGER = LogManager.getLogger(IclusionTrialFile.class);

    public static final String MAIN_FIELD_DELIMITER = "\t";
    private static final String SUB_FIELD_DELIMITER = "|";
    private static final String SUB_FIELD_SEPARATOR = " - ";
    private static final String DOID_DELIMITER = ",";

    private static final String NO_LIST_ENTRIES = "-";

    private IclusionTrialFile() {
    }

    public static void write(@NotNull String iClusionTrialsTsv, @NotNull List<IclusionTrial> trials) throws IOException {
        Files.write(new File(iClusionTrialsTsv).toPath(), toLines(trials));
    }

    public static List<IclusionTrial> read(@NotNull String iClusionTrialsTsv) throws IOException {
        return fromLines(Files.readAllLines(new File(iClusionTrialsTsv).toPath()));
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull List<IclusionTrial> trials) {
        List<String> lines = Lists.newArrayList(header());
        for (IclusionTrial trial : trials) {
            lines.add(toString(trial));
        }
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(MAIN_FIELD_DELIMITER).add("id")
                .add("acronym")
                .add("title")
                .add("eudra")
                .add("nct")
                .add("ipn")
                .add("ccmo")
                .add("tumorLocations")
                .add("mutations")
                .add("type")
                .add("age")
                .add("phase")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull IclusionTrial trial) {
        return new StringJoiner(MAIN_FIELD_DELIMITER).add(trial.id())
                .add(trial.acronym())
                .add(trial.title())
                .add(trial.eudra())
                .add(trial.nct())
                .add(trial.ipn())
                .add(trial.ccmo())
                .add(tumorLocationsToString(trial.tumorLocations()))
                .add(mutationsToString(trial.mutations()))
                .add(trial.type())
                .add(trial.age())
                .add(trial.phase())
                .toString();
    }

    @NotNull
    private static String tumorLocationsToString(@NotNull List<IclusionTumorLocation> tumorLocations) {
        StringJoiner tumorLocationString = new StringJoiner(SUB_FIELD_DELIMITER);

        for (IclusionTumorLocation tumorLocation : tumorLocations) {
            StringJoiner doidString = new StringJoiner(DOID_DELIMITER);
            for (String doid : tumorLocation.doids()) {
                if (!containsInvalidString(doid, DOID_DELIMITER)) {
                    doidString.add(doid);
                }
            }

            if (!containsInvalidString(tumorLocation.primaryTumorLocation(), SUB_FIELD_DELIMITER)
                    && !containsInvalidString(tumorLocation.primaryTumorLocation(), SUB_FIELD_SEPARATOR)) {
                tumorLocationString.add(tumorLocation.primaryTumorLocation() + SUB_FIELD_SEPARATOR + doidString.toString());
            }
        }

        String finalString = tumorLocationString.toString();

        return !finalString.isEmpty() ? finalString : NO_LIST_ENTRIES;
    }

    @NotNull
    private static String mutationsToString(@NotNull List<IclusionMutation> mutations) {
        StringJoiner mutationString = new StringJoiner(SUB_FIELD_DELIMITER);

        for (IclusionMutation mutation : mutations) {
            if (!containsInvalidString(mutation.gene(), SUB_FIELD_SEPARATOR) && !containsInvalidString(mutation.name(),
                    SUB_FIELD_SEPARATOR)) {
                mutationString.add(mutation.gene() + SUB_FIELD_SEPARATOR + mutation.name());
            }
        }

        String finalString = mutationString.toString();

        return !finalString.isEmpty() ? finalString : NO_LIST_ENTRIES;
    }

    private static boolean containsInvalidString(@NotNull String stringToCheck, @NotNull String invalidString) {
        if (stringToCheck.contains(invalidString)) {
            LOGGER.error("'{}' found in field '{}'. Cannot serialize!", invalidString, stringToCheck);
            return true;
        }

        return false;
    }

    @NotNull
    @VisibleForTesting
    static List<IclusionTrial> fromLines(@NotNull List<String> lines) {
        List<IclusionTrial> trials = Lists.newArrayList();

        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            trials.add(fromString(line));
        }

        return trials;
    }

    @NotNull
    private static IclusionTrial fromString(@NotNull String line) {
        String[] values = line.split(MAIN_FIELD_DELIMITER);

        return ImmutableIclusionTrial.builder()
                .id(values[0])
                .acronym(values[1])
                .title(values[2])
                .eudra(values[3])
                .nct(values[4])
                .ipn(values[5])
                .ccmo(values[6])
                .tumorLocations(tumorLocationsFromString(values[7]))
                .mutations(mutationsFromString(values[8]))
                .type(values[9])
                .age(values[10])
                .phase(values[11])
                .build();
    }

    @NotNull
    private static List<IclusionTumorLocation> tumorLocationsFromString(@NotNull String tumorLocationsString) {
        List<IclusionTumorLocation> tumorLocations = Lists.newArrayList();

        if (tumorLocationsString.equals(NO_LIST_ENTRIES)) {
            return tumorLocations;
        }

        String[] values = tumorLocationsString.split("\\" + SUB_FIELD_DELIMITER);
        for (String value : values) {
            String[] fields = value.split(SUB_FIELD_SEPARATOR);
            tumorLocations.add(ImmutableIclusionTumorLocation.builder()
                    .primaryTumorLocation(fields[0])
                    .doids(fields.length > 1 ? Lists.newArrayList(fields[1].split(DOID_DELIMITER)) : Lists.newArrayList())
                    .build());
        }

        return tumorLocations;
    }

    @NotNull
    private static List<IclusionMutation> mutationsFromString(@NotNull String mutationsString) {
        List<IclusionMutation> mutations = Lists.newArrayList();

        if (mutationsString.equals(NO_LIST_ENTRIES)) {
            return mutations;
        }

        String[] values = mutationsString.split("\\" + SUB_FIELD_DELIMITER);
        for (String value : values) {
            String[] fields = value.split(SUB_FIELD_SEPARATOR);
            mutations.add(ImmutableIclusionMutation.builder().gene(fields[0]).name(fields[1]).build());
        }

        return mutations;
    }
}
