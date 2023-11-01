package com.hartwig.hmftools.virusinterpreter.algo;

import static com.hartwig.hmftools.virusinterpreter.VirusInterpreterApplication.VI_LOGGER;

import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;

import org.jetbrains.annotations.Nullable;

public class VirusReportingDbModel
{
    private final Map<Integer, VirusReportingDb> speciesToInterpretationMap;

    public VirusReportingDbModel(final Map<Integer, VirusReportingDb> speciesToInterpretationMap)
    {
        this.speciesToInterpretationMap = speciesToInterpretationMap;
    }

    @Nullable
    public VirusType interpretVirusSpecies(int speciesTaxid)
    {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);
        if(!speciesHasInterpretation)
        {
            VI_LOGGER.debug("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        VirusReportingDb virusReportingDb = speciesToInterpretationMap.get(speciesTaxid);

        return speciesHasInterpretation ? virusReportingDb.virusInterpretation() : null;
    }

    public VirusLikelihoodType virusLikelihoodType(int speciesTaxid)
    {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);
        if(!speciesHasInterpretation)
        {
            VI_LOGGER.debug("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        VirusReportingDb virusReportingDb = speciesToInterpretationMap.get(speciesTaxid);

        return speciesHasInterpretation ? virusReportingDb.virusDriverLikelihoodType() : VirusLikelihoodType.UNKNOWN;
    }

    @Nullable
    public Integer integratedMinimalCoverage(int speciesTaxid)
    {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);
        if(!speciesHasInterpretation)
        {
            VI_LOGGER.debug("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        VirusReportingDb virusReportingDb = speciesToInterpretationMap.get(speciesTaxid);

        return speciesHasInterpretation ? virusReportingDb.integratedMinimalCoverage() : null;
    }

    @Nullable
    public Integer nonIntegratedMinimalCoverage(int speciesTaxid)
    {
        boolean speciesHasInterpretation = hasInterpretation(speciesTaxid);
        if(!speciesHasInterpretation)
        {
            VI_LOGGER.debug("No interpretation found for virus with species taxid {}", speciesTaxid);
        }

        VirusReportingDb virusReportingDb = speciesToInterpretationMap.get(speciesTaxid);

        return speciesHasInterpretation ? virusReportingDb.nonIntegratedMinimalCoverage() : null;
    }

    public boolean hasInterpretation(int speciesTaxid)
    {
        return speciesToInterpretationMap.containsKey(speciesTaxid);
    }

    @VisibleForTesting
    public int count()
    {
        return speciesToInterpretationMap.keySet().size();
    }
}
