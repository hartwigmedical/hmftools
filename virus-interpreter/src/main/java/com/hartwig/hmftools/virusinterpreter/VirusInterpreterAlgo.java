package com.hartwig.hmftools.virusinterpreter;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.ImmutableAnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusBreakend;
import com.hartwig.hmftools.common.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.virusinterpreter.algo.VirusReportingModel;
import com.hartwig.hmftools.virusinterpreter.coverages.CoveragesAnalysis;
import com.hartwig.hmftools.virusinterpreter.taxonomy.TaxonomyDb;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class VirusInterpreterAlgo {
    private static final Logger LOGGER = LogManager.getLogger(VirusInterpreterAlgo.class);

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
                    .build());
        }

        return annotatedViruses;
    }

    private boolean report(@NotNull VirusBreakend virusBreakend, double expectedClonalCoverage) {
        double viralPercentageCovered = virusBreakend.coverage();
        double viralCoverage = virusBreakend.meanDepth();

        Integer coverage = determineMinimalCoverageVirus(virusBreakend.integrations(), virusBreakend.taxidSpecies());

        boolean reported = false;
        boolean virusQCStatus = false;
        if (virusReportingModel.hasInterpretation(virusBreakend.taxidSpecies())) {
            if (virusBreakend.qcStatus() != VirusBreakendQCStatus.LOW_VIRAL_COVERAGE) {
                virusQCStatus = true;
            }

            if (coverage != null) {
                if (viralPercentageCovered > coverage && viralCoverage > expectedClonalCoverage) {
                    reported = true;
                }
            } else {
                reported = true;
            }
        }
        return reported && virusQCStatus;
    }

    @Nullable
    public Integer determineMinimalCoverageVirus(int integrations, int taxidSpecies) {
        if (integrations >= 1) {
            return virusReportingModel.integratedMinimalCoverage(taxidSpecies);
        } else {
            return virusReportingModel.nonIntegratedMinimalCoverage(taxidSpecies);
        }
    }
}
