package com.hartwig.hmftools.serve.blacklisting;

import java.util.List;
import java.util.StringJoiner;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TumorLocationBlacklist {

    private static final String DELIMITER_TUM_DOID = ",";
    private static final String DELIMITER_DIFFERENT_TUMOR = ";";

    private TumorLocationBlacklist() {
    }

    @NotNull
    public static String extractTumorLocationBlacklisting(@NotNull List<TumorLocationBlacklisting> tumorLocationBlacklistings) {
        StringJoiner joiner = new StringJoiner(DELIMITER_DIFFERENT_TUMOR);
        for (TumorLocationBlacklisting tumorLocation : tumorLocationBlacklistings) {
            if (tumorLocation.blacklistCancerType() != null || tumorLocation.blacklistedDoid() != null) {
                if (!tumorLocation.blacklistCancerType().equals(Strings.EMPTY) && !tumorLocation.blacklistedDoid().equals(Strings.EMPTY)) {
                    joiner.add(tumorLocation.blacklistCancerType() + DELIMITER_TUM_DOID + tumorLocation.blacklistedDoid());
                }
            }
        }
        return joiner.toString();
    }

    @NotNull
    public static List<TumorLocationBlacklisting> readTumorLocationBlacklistingString(@NotNull String tumorLocationBlacklistings) {
        List<TumorLocationBlacklisting> tumorLocationBlacklistingsList = Lists.newArrayList();
        String[] splitTumorLocationString = tumorLocationBlacklistings.split(DELIMITER_DIFFERENT_TUMOR);

        for (String tumorLocationDoidArray : splitTumorLocationString) {
            String[] tumorLocationDoid = tumorLocationDoidArray.split(DELIMITER_TUM_DOID);
            tumorLocationBlacklistingsList.add(ImmutableTumorLocationBlacklisting.builder()
                    .blacklistCancerType(tumorLocationDoid[0])
                    .blacklistedDoid(tumorLocationDoid[1])
                    .build());
        }
        return tumorLocationBlacklistingsList;
    }
}