package com.hartwig.hmftools.patientreporter.virusbreakend;

import java.util.Map;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class VirusBlacklistModel {

    @NotNull
    private final Map<Integer, String> virusBlacklistMap;

    VirusBlacklistModel(@NotNull final Map<Integer, String> virusBlacklistMap) {
        this.virusBlacklistMap = virusBlacklistMap;
    }

    public boolean checkVirusForBlacklisting(int id) {
        return virusBlacklistMap.containsKey(id);
    }

    @Nullable
    public String checkTaxusForId(int id) {
        return virusBlacklistMap.get(id);
    }

    @VisibleForTesting
    int count() {
        return virusBlacklistMap.keySet().size();
    }

}
