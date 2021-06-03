package com.hartwig.hmftools.virusinterpreter.algo;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.virus.VirusBreakend;
import com.hartwig.hmftools.common.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.common.virus.VirusInterpretation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ReportableVirusBreakendFactory {

    private static final Logger LOGGER = LogManager.getLogger(ReportableVirusBreakendFactory.class);

    @NotNull
    private final TaxonomyDb taxonomyDb;
    @NotNull
    private final VirusInterpretationModel virusInterpretationModel;
    @NotNull
    private final VirusBlacklistModel virusBlacklistModel;

    public ReportableVirusBreakendFactory(@NotNull final TaxonomyDb taxonomyDb,
            @NotNull final VirusInterpretationModel virusInterpretationModel, @NotNull final VirusBlacklistModel virusBlacklistModel) {
        this.taxonomyDb = taxonomyDb;
        this.virusInterpretationModel = virusInterpretationModel;
        this.virusBlacklistModel = virusBlacklistModel;
    }

    @NotNull
    public List<ReportableVirusBreakend> analyze(@NotNull List<VirusBreakend> virusBreakends) {
        List<VirusBreakend> virusBreakendsFiltered = Lists.newArrayList();
        for (VirusBreakend virusBreakend : virusBreakends) {
            if (include(virusBreakend)) {
                virusBreakendsFiltered.add(virusBreakend);
            }
        }

        List<ReportableVirusBreakend> reportableVirusBreakends = Lists.newArrayList();
        for (VirusBreakend virusBreakend : virusBreakendsFiltered) {
            String virusName = taxonomyDb.lookupName(virusBreakend.referenceTaxid());

            VirusInterpretation interpretation = null;
            if (virusInterpretationModel.hasInterpretation(virusBreakend.taxidSpecies())) {
                interpretation = virusInterpretationModel.interpretVirusSpecies(virusBreakend.taxidSpecies());
            } else {
                LOGGER.info(" VIRUS breakend has called a non-interpreted virus: {}", virusName);
            }

            reportableVirusBreakends.add(ImmutableReportableVirusBreakend.builder()
                    .virusName(virusName)
                    .integrations(virusBreakend.integrations())
                    .interpretation(interpretation)
                    .build());
        }

        return reportableVirusBreakends;
    }

    private boolean include(@NotNull VirusBreakend virusBreakend) {
        if (virusBreakend.qcStatus() == VirusBreakendQCStatus.LOW_VIRAL_COVERAGE) {
            return false;
        }

        if (virusBreakend.integrations() == 0) {
            return false;
        }

        return !virusBlacklistModel.isBlacklisted(virusBreakend);
    }
}
