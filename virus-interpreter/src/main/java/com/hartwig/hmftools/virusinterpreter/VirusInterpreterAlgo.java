package com.hartwig.hmftools.virusinterpreter;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.ImmutableAnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusBreakend;
import com.hartwig.hmftools.common.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.virusinterpreter.algo.VirusReportingDbModel;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;
import com.hartwig.hmftools.virusinterpreter.coverages.CoveragesAnalysis;
import com.hartwig.hmftools.virusinterpreter.taxonomy.TaxonomyDb;

import org.jetbrains.annotations.Nullable;

public class VirusInterpreterAlgo
{
    private final TaxonomyDb taxonomyDb;
    private final VirusReportingDbModel virusReportingDbModel;
    private final CoveragesAnalysis coveragesAnalysis;

    public VirusInterpreterAlgo(
            final TaxonomyDb taxonomyDb, final VirusReportingDbModel virusReportingDbModel, final CoveragesAnalysis coveragesAnalysis)
    {
        this.taxonomyDb = taxonomyDb;
        this.virusReportingDbModel = virusReportingDbModel;
        this.coveragesAnalysis = coveragesAnalysis;
    }

    public List<AnnotatedVirus> analyze(final List<VirusBreakend> virusBreakends, final PurityContext purityContext)
    {
        List<AnnotatedVirus> annotatedViruses = Lists.newArrayList();
        for(VirusBreakend virusBreakend : virusBreakends)
        {
            VirusType interpretation = virusReportingDbModel.interpretVirusSpecies(virusBreakend.taxidSpecies());

            int taxid = virusBreakend.referenceTaxid();
            boolean reported = report(virusBreakend, coveragesAnalysis.ExpectedClonalCoverage, purityContext.qc().status());
            annotatedViruses.add(ImmutableAnnotatedVirus.builder()
                    .taxid(taxid)
                    .name(taxonomyDb.lookupName(taxid))
                    .qcStatus(virusBreakend.qcStatus())
                    .integrations(virusBreakend.integrations())
                    .interpretation(interpretation)
                    .percentageCovered(virusBreakend.coverage())
                    .meanCoverage(virusBreakend.meanDepth())
                    .expectedClonalCoverage(hasAcceptablePurpleQuality(purityContext.qc().status())
                            ? coveragesAnalysis.ExpectedClonalCoverage
                            : null)
                    .reported(reported)
                    .virusDriverLikelihoodType(virusLikelihoodType(virusBreakend, reported))
                    .build());
        }

        return annotatedViruses;
    }

    @VisibleForTesting
    VirusLikelihoodType virusLikelihoodType(VirusBreakend virusBreakend, Boolean reported)
    {
        return reported ? virusReportingDbModel.virusLikelihoodType(virusBreakend.taxidSpecies()) : VirusLikelihoodType.UNKNOWN;
    }

    @VisibleForTesting
    boolean report(VirusBreakend virusBreakend, double expectedClonalCoverage, final Set<PurpleQCStatus> purpleQCStatuses)
    {
        boolean reported = false;
        if(virusReportingDbModel.hasInterpretation(virusBreakend.taxidSpecies()))
        {
            Integer minimalCoveragePercentage = determineMinimalCoverageVirus(virusBreakend.integrations(), virusBreakend.taxidSpecies());
            double viralPercentageCovered = virusBreakend.coverage();
            double viralCoverage = virusBreakend.meanDepth();

            if(hasAcceptablePurpleQuality(purpleQCStatuses))
            {
                if(minimalCoveragePercentage != null)
                {
                    if(viralPercentageCovered > minimalCoveragePercentage && viralCoverage > expectedClonalCoverage)
                    {
                        reported = true;
                    }
                }
                else
                {
                    reported = true;
                }
            }
            else
            {
                if(minimalCoveragePercentage == null)
                {
                    reported = true;
                }
            }
        }

        boolean virusQCStatus = virusBreakend.qcStatus() != VirusBreakendQCStatus.LOW_VIRAL_COVERAGE;
        return reported && virusQCStatus;
    }

    @Nullable
    @VisibleForTesting
    Integer determineMinimalCoverageVirus(int integrations, int taxidSpecies)
    {
        if(integrations >= 1)
        {
            return virusReportingDbModel.integratedMinimalCoverage(taxidSpecies);
        }
        else
        {
            return virusReportingDbModel.nonIntegratedMinimalCoverage(taxidSpecies);
        }
    }

    @VisibleForTesting
    static boolean hasAcceptablePurpleQuality(final Set<PurpleQCStatus> purpleQCStatus)
    {
        return !purpleQCStatus.contains(PurpleQCStatus.FAIL_NO_TUMOR) && !purpleQCStatus.contains(PurpleQCStatus.FAIL_CONTAMINATION);
    }
}
