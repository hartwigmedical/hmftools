package com.hartwig.hmftools.patientreporter.virusbreakend;

import java.util.Map;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class VirusBlackListModel {
    private static final Logger LOGGER = LogManager.getLogger(VirusDbModel.class);

    @NotNull
    private final Map<Integer, String> virusBlacklistMap;

    public VirusBlackListModel(@NotNull final Map<Integer, String> virusBlacklistMap) {
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
    int virusBlacklistcount() {
        return virusBlacklistMap.keySet().size();
    }

}
