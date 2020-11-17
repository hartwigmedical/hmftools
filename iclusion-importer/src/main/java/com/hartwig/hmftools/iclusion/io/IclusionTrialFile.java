package com.hartwig.hmftools.iclusion.io;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.iclusion.data.IclusionMutation;
import com.hartwig.hmftools.iclusion.data.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.data.IclusionMutationLogicType;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.iclusion.data.IclusionTumorLocation;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionMutation;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionMutationCondition;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionTrial;
import com.hartwig.hmftools.iclusion.data.ImmutableIclusionTumorLocation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class IclusionTrialFile {

    private static final Logger LOGGER = LogManager.getLogger(IclusionTrialFile.class);

    public static final String MAIN_FIELD_DELIMITER = "\t";
    private static final String SUB_FIELD_DELIMITER = "|";
    private static final String SUB_SUB_FIELD_DELIMITER = "#";
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
                .add("blacklistedTumorLocations")
                .add("mutationConditions")
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
                .add(tumorLocationsToString(trial.blacklistedTumorLocations()))
                .add(mutationConditionsToString(trial.mutationConditions()))
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
    private static String mutationConditionsToString(@NotNull List<IclusionMutationCondition> mutationConditions) {
        StringJoiner mutationConditionString = new StringJoiner(SUB_FIELD_DELIMITER);

        for (IclusionMutationCondition mutationCondition : mutationConditions) {
            StringJoiner mutationString = new StringJoiner(SUB_SUB_FIELD_DELIMITER);
            for (IclusionMutation mutation : mutationCondition.mutations()) {
                if (!containsInvalidString(mutation.name(), SUB_FIELD_SEPARATOR) && !containsInvalidString(mutation.gene(),
                        SUB_FIELD_SEPARATOR)) {
                    mutationString.add(mutation.gene() + SUB_FIELD_SEPARATOR + mutation.name() + SUB_FIELD_SEPARATOR + mutation.negation());
                }
            }
            String finalMutationString = mutationString.toString();
            if (!finalMutationString.isEmpty()) {
                mutationConditionString.add(mutationCondition.logicType().toString() + "(" + mutationString.toString() + ")");
            }
        }

        String finalString = mutationConditionString.toString();

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
                .blacklistedTumorLocations(tumorLocationsFromString(values[8]))
                .mutationConditions(mutationConditionsFromString(values[9]))
                .build();
    }

    @NotNull
    private static List<IclusionTumorLocation> tumorLocationsFromString(@NotNull String tumorLocationsString) {
        if (tumorLocationsString.equals(NO_LIST_ENTRIES)) {
            return Lists.newArrayList();
        }

        List<IclusionTumorLocation> tumorLocations = Lists.newArrayList();

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
    private static List<IclusionMutationCondition> mutationConditionsFromString(@NotNull String mutationConditionString) {
        if (mutationConditionString.equals(NO_LIST_ENTRIES)) {
            return Lists.newArrayList();
        }

        List<IclusionMutationCondition> mutationConditions = Lists.newArrayList();

        String[] values = mutationConditionString.split("\\" + SUB_FIELD_DELIMITER);
        for (String value : values) {
            // Format is "logicType(mutations)"
            int separatorPos = value.indexOf("(");
            IclusionMutationLogicType logicType = IclusionMutationLogicType.valueOf(value.substring(0, separatorPos));
            String mutationEntries = value.substring(separatorPos + 1, value.length() - 1);
            List<IclusionMutation> mutations = Lists.newArrayList();
            for (String mutationEntry : mutationEntries.split(SUB_SUB_FIELD_DELIMITER)) {
                String[] mutationFields = mutationEntry.split(SUB_FIELD_SEPARATOR);
                mutations.add(ImmutableIclusionMutation.builder()
                        .gene(mutationFields[0])
                        .name(mutationFields[1])
                        .negation(Boolean.parseBoolean(mutationFields[2]))
                        .build());
            }
            mutationConditions.add(ImmutableIclusionMutationCondition.builder().logicType(logicType).mutations(mutations).build());
        }

        return mutationConditions;
    }
}
