package com.hartwig.hmftools.virusinterpreter.algo;

import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class VirusReportingDbModel {

    private static final Logger LOGGER = LogManager.getLogger(VirusReportingDbModel.class);

    @NotNull
    private final Map<Integer, VirusReportingDb> speciesToInterpretationMap;

    public VirusReportingDbModel(@NotNull final Map<Integer, VirusReportingDb> speciesToInterpretationMap) {
        this.speciesToInterpretationMap = speciesToInterpretationMap;
    }

    @Nullable
    public String interpretVirusSpecies(int speciesTaxid) {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);
        if (!speciesHasInterpretation) {
            LOGGER.debug("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        VirusReportingDb virusReportingDb = speciesToInterpretationMap.get(speciesTaxid);

        return speciesHasInterpretation ? virusReportingDb.virusInterpretation() : null;
    }

    @NotNull
    public VirusLikelihoodType virusLikelihoodType(int speciesTaxid) {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);
        if (!speciesHasInterpretation) {
            LOGGER.debug("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        VirusReportingDb virusReportingDb = speciesToInterpretationMap.get(speciesTaxid);

        return speciesHasInterpretation ? virusReportingDb.virusDriverLikelihoodType() : VirusLikelihoodType.UNKNOWN;
    }

    @Nullable
    public Integer integratedMinimalCoverage(int speciesTaxid) {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);
        if (!speciesHasInterpretation) {
            LOGGER.debug("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        VirusReportingDb virusReportingDb = speciesToInterpretationMap.get(speciesTaxid);

        return speciesHasInterpretation ? virusReportingDb.integratedMinimalCoverage() : null;
    }

    @Nullable
    public Integer nonIntegratedMinimalCoverage(int speciesTaxid) {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);
        if (!speciesHasInterpretation) {
            LOGGER.debug("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        VirusReportingDb virusReportingDb = speciesToInterpretationMap.get(speciesTaxid);

        return speciesHasInterpretation ? virusReportingDb.nonIntegratedMinimalCoverage() : null;
    }

    public boolean hasInterpretation(int speciesTaxid) {
        return speciesToInterpretationMap.containsKey(speciesTaxid);
    }

    @VisibleForTesting
    public int count() {
        return speciesToInterpretationMap.keySet().size();
    }
}
