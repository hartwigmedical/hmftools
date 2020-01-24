package com.hartwig.hmftools.iclusion.io;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.iclusion.data.IclusionMutation;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.iclusion.data.IclusionTumorLocation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class IclusionTrialFile {

    private static final Logger LOGGER = LogManager.getLogger(IclusionTrialFile.class);

    private static final String MAIN_FIELD_DELIMITER = "\t";
    private static final String SUB_FIELD_DELIMITER = "|";
    private static final String SUB_FIELD_SEPARATOR = " - ";
    private static final String DOID_DELIMITER = ",";

    private IclusionTrialFile() {
    }

    public static void write(@NotNull String iClusionTrialsTsv, @NotNull List<IclusionTrial> trials) throws IOException {
        Files.write(new File(iClusionTrialsTsv).toPath(), toLines(trials));
    }

    @NotNull
    private static List<String> toLines(@NotNull List<IclusionTrial> trials) {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        trials.stream().map(IclusionTrialFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(MAIN_FIELD_DELIMITER, "", "").add("id")
                .add("acronym")
                .add("title")
                .add("eudra")
                .add("nct")
                .add("ipn")
                .add("ccmo")
                .add("tumorLocations")
                .add("mutations")
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
        return tumorLocationString.toString();
    }

    @NotNull
    private static String mutationsToString(@NotNull List<IclusionMutation> mutations) {
        StringJoiner mutationString = new StringJoiner(SUB_FIELD_DELIMITER);

        for (IclusionMutation mutation : mutations) {
            if (!containsInvalidString(mutation.gene(), SUB_FIELD_SEPARATOR) && !containsInvalidString(mutation.name(), SUB_FIELD_SEPARATOR)) {
                mutationString.add(mutation.gene() + SUB_FIELD_SEPARATOR + mutation.name());
            }
        }
        return mutationString.toString();
    }

    private static boolean containsInvalidString(@NotNull String stringToCheck, @NotNull String invalidString) {
        if (stringToCheck.contains(invalidString)) {
            LOGGER.error("{} found in primary tumor location '{}'. Cannot serialize", invalidString, stringToCheck);
            return true;
        }

        return false;
    }
}
