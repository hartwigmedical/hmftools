package com.hartwig.hmftools.virusinterpreter;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.ImmutableAnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusBreakend;
import com.hartwig.hmftools.common.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.virusinterpreter.algo.VirusReportingModel;
import com.hartwig.hmftools.virusinterpreter.coverages.CoveragesAnalysis;
import com.hartwig.hmftools.virusinterpreter.taxonomy.TaxonomyDb;

import org.jetbrains.annotations.NotNull;

public class VirusInterpreterAlgo {

    @NotNull
    private final TaxonomyDb taxonomyDb;
    @NotNull
    private final VirusReportingModel virusReportingModel;

    @NotNull
    private final CoveragesAnalysis coveragesAnalysis;

    public VirusInterpreterAlgo(@NotNull final TaxonomyDb taxonomyDb, @NotNull final VirusReportingModel virusReportingModel,
            @NotNull CoveragesAnalysis coveragesAnalysis) {
        this.taxonomyDb = taxonomyDb;
        this.virusReportingModel = virusReportingModel;
        this.coveragesAnalysis = coveragesAnalysis;
    }

    @NotNull
    public List<AnnotatedVirus> analyze(@NotNull List<VirusBreakend> virusBreakends) throws IOException {

        List<AnnotatedVirus> annotatedViruses = Lists.newArrayList();
        for (VirusBreakend virusBreakend : virusBreakends) {
            String interpretation = virusReportingModel.interpretVirusSpecies(virusBreakend.taxidSpecies());

            int taxid = virusBreakend.referenceTaxid();
            annotatedViruses.add(ImmutableAnnotatedVirus.builder()
                    .taxid(taxid)
                    .name(taxonomyDb.lookupName(taxid))
                    .qcStatus(virusBreakend.qcStatus())
                    .integrations(virusBreakend.integrations())
                    .interpretation(interpretation)
                    .percentageCovered(virusBreakend.coverage())
                    .coverage(virusBreakend.meanDepth())
                    .expectedClonalCoverage(coveragesAnalysis.expectedClonalCoverage())
                    .reported(report(virusBreakend, coveragesAnalysis.expectedClonalCoverage()))
                    .reportedSummary(virusReportingModel.displayVirusOnSummaryReport(taxid))
                    .build());
        }

        return annotatedViruses;
    }

    private boolean report(@NotNull VirusBreakend virusBreakend, double expectedClonalCoverage) {
        double viralPercentageCovered = virusBreakend.coverage();
        double viralCoverage = virusBreakend.meanDepth();

        if (virusReportingModel.hasInterpretation(virusBreakend.taxidSpecies())) {
            if (virusBreakend.qcStatus() == VirusBreakendQCStatus.LOW_VIRAL_COVERAGE) {
                return false;
            }

            if (virusBreakend.integrations() == 0) {
                Integer minimalCoverage = virusReportingModel.nonIntegratedMinimalCoverage(virusBreakend.taxidSpecies());
                if (minimalCoverage != null) {
                    if (viralPercentageCovered <= minimalCoverage && viralCoverage <= expectedClonalCoverage) {
                        return false;
                    }
                }
            }

            if (virusBreakend.integrations() >= 1) {
                Integer minimalCoverage = virusReportingModel.integratedMinimalCoverage(virusBreakend.taxidSpecies());
                if (minimalCoverage != null) {
                    if (viralPercentageCovered <= minimalCoverage && viralCoverage <= expectedClonalCoverage) {
                        return false;
                    }
                }
            }
        }

        return true;
    }
}
