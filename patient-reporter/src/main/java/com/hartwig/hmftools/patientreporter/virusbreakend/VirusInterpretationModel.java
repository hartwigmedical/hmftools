package com.hartwig.hmftools.patientreporter.virusbreakend;

import java.util.Collection;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class VirusInterpretationModel {

    private static final Logger LOGGER = LogManager.getLogger(VirusInterpretationModel.class);

    @NotNull
    private final Map<Integer, String> speciesToInterpretationMap;

    VirusInterpretationModel(@NotNull final Map<Integer, String> speciesToInterpretationMap) {
        this.speciesToInterpretationMap = speciesToInterpretationMap;
    }

    @Nullable
    public String interpretVirusSpecies(int speciesTaxid) {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);

        if (!speciesHasInterpretation) {
            LOGGER.debug("No interpretation found for virus with species taxid {}", speciesTaxid);
        }
        return speciesHasInterpretation ? speciesToInterpretationMap.get(speciesTaxid) : null;
    }

    boolean hasInterpretation(int speciesTaxid) {
        return speciesToInterpretationMap.containsKey(speciesTaxid);
    }

    @VisibleForTesting
    int count() {
        return speciesToInterpretationMap.keySet().size();
    }
}
