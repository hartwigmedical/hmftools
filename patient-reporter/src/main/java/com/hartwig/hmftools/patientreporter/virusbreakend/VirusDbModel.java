package com.hartwig.hmftools.patientreporter.virusbreakend;

import java.util.Map;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class VirusDbModel {

    private static final Logger LOGGER = LogManager.getLogger(VirusDbModel.class);

    @NotNull
    private final Map<Integer, String> virusDbMap;

    public VirusDbModel(@NotNull final Map<Integer, String> virusDbMap) {
        this.virusDbMap = virusDbMap;
    }

    @NotNull
    public String findVirus(int id) {
        boolean mappedVirus = mapIdtoVirusName(id);

        if (!mappedVirus) {
            LOGGER.warn("Could not match id to virusName");
        }
        return mappedVirus ? virusDbMap.get(id) : Strings.EMPTY;
    }

    @VisibleForTesting
    boolean mapIdtoVirusName(int id) {
        return virusDbMap.containsKey(id);
    }
}
