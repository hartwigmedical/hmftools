package com.hartwig.hmftools.virusinterpreter.algo;

import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.virus.VirusInterpretation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class VirusWhitelistModel {

    private static final Logger LOGGER = LogManager.getLogger(VirusWhitelistModel.class);

    @NotNull
    private final Map<Integer, VirusWhitelist> speciesToInterpretationMap;

    public VirusWhitelistModel(@NotNull final Map<Integer, VirusWhitelist> speciesToInterpretationMap) {
        this.speciesToInterpretationMap = speciesToInterpretationMap;
    }

    @Nullable
    public VirusInterpretation interpretVirusSpecies(int speciesTaxid) {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);

        VirusWhitelist virusWhitelist = speciesToInterpretationMap.get(speciesTaxid);
        if (!speciesHasInterpretation) {
            LOGGER.error("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        return speciesHasInterpretation ? virusWhitelist.virusInterpretation() : null;
    }

    @Nullable
    public Boolean displayVirusOnreport(int speciesTaxid) {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);

        VirusWhitelist virusWhitelist = speciesToInterpretationMap.get(speciesTaxid);
        if (!speciesHasInterpretation) {
            LOGGER.error("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        return speciesHasInterpretation ? virusWhitelist.reportOnSummary() : null;
    }

    @Nullable
    public Integer integratedMinimalCoverage(int speciesTaxid) {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);

        VirusWhitelist virusWhitelist = speciesToInterpretationMap.get(speciesTaxid);
        if (!speciesHasInterpretation) {
            LOGGER.error("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        return speciesHasInterpretation ? virusWhitelist.integratedMinimalCoverage() : null;
    }

    @Nullable
    public Integer nonintegratedMinimalCoverage(int speciesTaxid) {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);

        VirusWhitelist virusWhitelist = speciesToInterpretationMap.get(speciesTaxid);
        if (!speciesHasInterpretation) {
            LOGGER.error("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        return speciesHasInterpretation ? virusWhitelist.nonintegratedMinimalCoverage() : null;
    }

    public boolean hasInterpretation(int speciesTaxid) {
        return speciesToInterpretationMap.containsKey(speciesTaxid);
    }

    @VisibleForTesting
    int count() {
        return speciesToInterpretationMap.keySet().size();
    }
}
