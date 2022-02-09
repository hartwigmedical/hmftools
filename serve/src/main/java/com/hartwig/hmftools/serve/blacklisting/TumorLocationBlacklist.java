package com.hartwig.hmftools.serve.blacklisting;

import java.util.Set;
import java.util.StringJoiner;

import org.jetbrains.annotations.NotNull;

public class TumorLocationBlacklist {

    private static final String DELIMITER = ",";

    @NotNull
    public static String extractTumorLocationBlacklisting(@NotNull Set<TumorLocationBlacklisting> tumorLocationBlacklistings) {
        StringJoiner joiner = new StringJoiner(DELIMITER);
        for (TumorLocationBlacklisting tumorLocation: tumorLocationBlacklistings) {
            joiner.add(tumorLocation.blacklistCancerType());
        }
        return joiner.toString();
    }

    @NotNull
    public static String extractTumorLocationDoid(@NotNull Set<TumorLocationBlacklisting> tumorLocationBlacklistings) {
        StringJoiner joiner = new StringJoiner(DELIMITER);
        for (TumorLocationBlacklisting tumorLocation: tumorLocationBlacklistings) {
            joiner.add(tumorLocation.blacklistedDoid());
        }
        return joiner.toString();
    }
}
