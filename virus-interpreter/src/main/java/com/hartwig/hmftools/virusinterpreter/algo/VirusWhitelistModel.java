package com.hartwig.hmftools.virusinterpreter.algo;

import java.util.Map;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class VirusWhitelistModel {

    private static final Logger LOGGER = LogManager.getLogger(VirusWhitelistModel.class);

    @NotNull
    private final Map<Integer, VirusWhitelist> speciesToInterpretationMap;

    public VirusWhitelistModel(@NotNull final Map<Integer, VirusWhitelist> speciesToInterpretationMap) {
        this.speciesToInterpretationMap = speciesToInterpretationMap;
    }

    @NotNull
    public String interpretVirusSpecies(int speciesTaxid) {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);

        VirusWhitelist virusWhitelist = speciesToInterpretationMap.get(speciesTaxid);
        if (!speciesHasInterpretation) {
            LOGGER.warn("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        return speciesHasInterpretation ? virusWhitelist.virusInterpretation() : Strings.EMPTY;
    }

    public boolean displayVirusOnSummaryReport(int speciesTaxid) {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);

        VirusWhitelist virusWhitelist = speciesToInterpretationMap.get(speciesTaxid);
        if (!speciesHasInterpretation) {
            LOGGER.warn("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        return speciesHasInterpretation && virusWhitelist.reportOnSummary();
    }

    @Nullable
    public Integer integratedMinimalCoverage(int speciesTaxid) {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);

        VirusWhitelist virusWhitelist = speciesToInterpretationMap.get(speciesTaxid);
        if (!speciesHasInterpretation) {
            LOGGER.warn("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        return speciesHasInterpretation ? virusWhitelist.integratedMinimalCoverage() : null;
    }

    @Nullable
    public Integer nonintegratedMinimalCoverage(int speciesTaxid) {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);

        VirusWhitelist virusWhitelist = speciesToInterpretationMap.get(speciesTaxid);
        if (!speciesHasInterpretation) {
            LOGGER.warn("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        return speciesHasInterpretation ? virusWhitelist.nonintegratedMinimalCoverage() : null;
    }

    public boolean hasInterpretation(int speciesTaxid) {
        return speciesToInterpretationMap.containsKey(speciesTaxid);
    }

    @VisibleForTesting
    public int count() {
        return speciesToInterpretationMap.keySet().size();
    }
}
