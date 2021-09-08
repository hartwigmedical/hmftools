package com.hartwig.hmftools.virusinterpreter.algo;

import java.util.Map;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class VirusReportingModel {

    private static final Logger LOGGER = LogManager.getLogger(VirusReportingModel.class);

    @NotNull
    private final Map<Integer, VirusReporting> speciesToInterpretationMap;

    public VirusReportingModel(@NotNull final Map<Integer, VirusReporting> speciesToInterpretationMap) {
        this.speciesToInterpretationMap = speciesToInterpretationMap;
    }

    @NotNull
    public String interpretVirusSpecies(int speciesTaxid) {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);

        VirusReporting virusReporting = speciesToInterpretationMap.get(speciesTaxid);
        if (!speciesHasInterpretation) {
            LOGGER.warn("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        return speciesHasInterpretation ? virusReporting.virusInterpretation() : Strings.EMPTY;
    }

    public boolean displayVirusOnSummaryReport(int speciesTaxid) {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);

        VirusReporting virusReporting = speciesToInterpretationMap.get(speciesTaxid);
        if (!speciesHasInterpretation) {
            LOGGER.warn("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        return speciesHasInterpretation && virusReporting.reportOnSummary();
    }

    @Nullable
    public Integer integratedMinimalCoverage(int speciesTaxid) {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);

        VirusReporting virusReporting = speciesToInterpretationMap.get(speciesTaxid);
        if (!speciesHasInterpretation) {
            LOGGER.warn("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        return speciesHasInterpretation ? virusReporting.integratedMinimalCoverage() : null;
    }

    @Nullable
    public Integer nonIntegratedMinimalCoverage(int speciesTaxid) {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);

        VirusReporting virusReporting = speciesToInterpretationMap.get(speciesTaxid);
        if (!speciesHasInterpretation) {
            LOGGER.warn("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        return speciesHasInterpretation ? virusReporting.nonIntegratedMinimalCoverage() : null;
    }

    public boolean hasInterpretation(int speciesTaxid) {
        return speciesToInterpretationMap.containsKey(speciesTaxid);
    }

    @VisibleForTesting
    public int count() {
        return speciesToInterpretationMap.keySet().size();
    }
}
