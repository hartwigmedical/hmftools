package com.hartwig.hmftools.serve.tumorlocation;

import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Sets;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TumorLocationFactory {

    private static final String DELIMITER_TUM_DOID = ",";
    private static final String DELIMITER_DIFFERENT_TUMOR = ";";

    private TumorLocationFactory() {
    }

    @NotNull
    public static String extractTumorLocationBlacklisting(@NotNull Set<TumorLocation> tumorLocationBlacklistings) {
        StringJoiner joiner = new StringJoiner(DELIMITER_DIFFERENT_TUMOR);
        for (TumorLocation blackListTumorLocation : tumorLocationBlacklistings) {
            if (!blackListTumorLocation.cancerType().equals(Strings.EMPTY) && !blackListTumorLocation.doid().equals(Strings.EMPTY)) {
                joiner.add(blackListTumorLocation.cancerType() + DELIMITER_TUM_DOID + blackListTumorLocation.doid());
            }
        }
        return joiner.toString();
    }

    @NotNull
    public static Set<TumorLocation> readTumorLocationBlacklistingString(@NotNull String tumorLocationBlacklistings) {
        Set<TumorLocation> tumorLocationBlacklistingsList = Sets.newConcurrentHashSet();
        if (!tumorLocationBlacklistings.isEmpty()) {
            String[] splitTumorLocationString = tumorLocationBlacklistings.split(DELIMITER_DIFFERENT_TUMOR);

            for (String tumorLocationDoidArray : splitTumorLocationString) {
                String[] tumorLocationDoid = tumorLocationDoidArray.split(DELIMITER_TUM_DOID);
                tumorLocationBlacklistingsList.add(ImmutableTumorLocation.builder()
                        .cancerType(tumorLocationDoid[0])
                        .doid(tumorLocationDoid[1])
                        .build());
            }
        }

        return tumorLocationBlacklistingsList;
    }
}