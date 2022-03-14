package com.hartwig.hmftools.serve.cancertype;

import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Sets;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class CancerTypeFactory {

    private static final String DELIMITER_TUM_DOID = ",";
    private static final String DELIMITER_DIFFERENT_TUMOR = ";";

    private CancerTypeFactory() {
    }

    @NotNull
    public static String extractCancerTypeBlacklist(@NotNull Set<CancerType> tumorLocationBlacklistings) {
        StringJoiner joiner = new StringJoiner(DELIMITER_DIFFERENT_TUMOR);
        for (CancerType blackListTumorLocation : tumorLocationBlacklistings) {
            if (!blackListTumorLocation.cancerType().equals(Strings.EMPTY) && !blackListTumorLocation.doid().equals(Strings.EMPTY)) {
                joiner.add(blackListTumorLocation.cancerType() + DELIMITER_TUM_DOID + blackListTumorLocation.doid());
            }
        }
        return joiner.toString();
    }

    @NotNull
    public static Set<CancerType> readCancerTypeBlacklistString(@NotNull String tumorLocationBlacklistings) {
        Set<CancerType> tumorLocationBlacklistingsList = Sets.newConcurrentHashSet();
        if (!tumorLocationBlacklistings.isEmpty()) {
            String[] splitTumorLocationString = tumorLocationBlacklistings.split(DELIMITER_DIFFERENT_TUMOR);

            for (String tumorLocationDoidArray : splitTumorLocationString) {
                String[] tumorLocationDoid = tumorLocationDoidArray.split(DELIMITER_TUM_DOID);
                tumorLocationBlacklistingsList.add(ImmutableCancerType.builder()
                        .cancerType(tumorLocationDoid[0])
                        .doid(tumorLocationDoid[1])
                        .build());
            }
        }

        return tumorLocationBlacklistingsList;
    }
}