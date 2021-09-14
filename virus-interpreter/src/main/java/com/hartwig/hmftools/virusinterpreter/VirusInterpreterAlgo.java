package com.hartwig.hmftools.virusinterpreter;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.ImmutableAnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusBreakend;
import com.hartwig.hmftools.common.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.virusinterpreter.algo.VirusReportingDbModel;
import com.hartwig.hmftools.virusinterpreter.coverages.CoveragesAnalysis;
import com.hartwig.hmftools.virusinterpreter.taxonomy.TaxonomyDb;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class VirusInterpreterAlgo {

    @NotNull
    private final TaxonomyDb taxonomyDb;
    @NotNull
    private final VirusReportingDbModel virusReportingDbModel;
    @NotNull
    private final CoveragesAnalysis coveragesAnalysis;

    public VirusInterpreterAlgo(@NotNull final TaxonomyDb taxonomyDb, @NotNull final VirusReportingDbModel virusReportingDbModel,
            @NotNull final CoveragesAnalysis coveragesAnalysis) {
        this.taxonomyDb = taxonomyDb;
        this.virusReportingDbModel = virusReportingDbModel;
        this.coveragesAnalysis = coveragesAnalysis;
    }

    public boolean isPurpleQcPass(@NotNull Set<PurpleQCStatus> purpleQCStatus) {
        return !purpleQCStatus.contains(PurpleQCStatus.FAIL_NO_TUMOR) && !purpleQCStatus.contains(PurpleQCStatus.FAIL_CONTAMINATION);
    }

    @NotNull
    public List<AnnotatedVirus> analyze(@NotNull List<VirusBreakend> virusBreakends, @NotNull PurityContext purityContext) {
        List<AnnotatedVirus> annotatedViruses = Lists.newArrayList();
        for (VirusBreakend virusBreakend : virusBreakends) {
            String interpretation = virusReportingDbModel.interpretVirusSpecies(virusBreakend.taxidSpecies());

            int taxid = virusBreakend.referenceTaxid();
            annotatedViruses.add(ImmutableAnnotatedVirus.builder()
                    .taxid(taxid)
                    .name(taxonomyDb.lookupName(taxid))
                    .qcStatus(virusBreakend.qcStatus())
                    .integrations(virusBreakend.integrations())
                    .interpretation(interpretation)
                    .percentageCovered(virusBreakend.coverage())
                    .meanCoverage(virusBreakend.meanDepth())
                    .expectedClonalCoverage(isPurpleQcPass(purityContext.qc().status()) ? coveragesAnalysis.expectedClonalCoverage() : null)
                    .reported(report(virusBreakend, coveragesAnalysis.expectedClonalCoverage(), purityContext.qc().status())).build());
        }

        return annotatedViruses;
    }

    @VisibleForTesting
    public boolean report(@NotNull VirusBreakend virusBreakend, double expectedClonalCoverage, @NotNull Set<PurpleQCStatus> purpleQCStatuses) {
        boolean reported = false;
        if (virusReportingDbModel.hasInterpretation(virusBreakend.taxidSpecies())) {
            Integer minimalCoveragePercentage = determineMinimalCoverageVirus(virusBreakend.integrations(), virusBreakend.taxidSpecies());
            double viralPercentageCovered = virusBreakend.coverage();
            double viralCoverage = virusBreakend.meanDepth();

            if (isPurpleQcPass(purpleQCStatuses)) {
                if (minimalCoveragePercentage != null) {
                    if (viralPercentageCovered > minimalCoveragePercentage && viralCoverage > expectedClonalCoverage) {
                        reported = true;
                    }
                } else {
                    reported = true;
                }
            } else {
                if (minimalCoveragePercentage == null) {
                    reported = true;
                }
            }
        }

        boolean virusQCStatus = virusBreakend.qcStatus() != VirusBreakendQCStatus.LOW_VIRAL_COVERAGE;
        return reported && virusQCStatus;
    }

    @Nullable
    public Integer determineMinimalCoverageVirus(int integrations, int taxidSpecies) {
        if (integrations >= 1) {
            return virusReportingDbModel.integratedMinimalCoverage(taxidSpecies);
        } else {
            return virusReportingDbModel.nonIntegratedMinimalCoverage(taxidSpecies);
        }
    }
}
