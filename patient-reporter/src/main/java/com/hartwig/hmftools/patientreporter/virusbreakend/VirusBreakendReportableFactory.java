package com.hartwig.hmftools.patientreporter.virusbreakend;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.virusbreakend.VirusBreakend;
import com.hartwig.hmftools.common.virusbreakend.VirusBreakendQCStatus;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class VirusBreakendReportableFactory {

    private static final Logger LOGGER = LogManager.getLogger(VirusBreakendReportableFactory.class);

    @NotNull
    private final TaxonomyDb taxonomyDb;
    @NotNull
    private final VirusInterpretationModel virusInterpretationModel;
    @NotNull
    private final VirusBlacklistModel virusBlacklistModel;

    public VirusBreakendReportableFactory(@NotNull final TaxonomyDb taxonomyDb,
            @NotNull final VirusInterpretationModel virusInterpretationModel, @NotNull final VirusBlacklistModel virusBlacklistModel) {
        this.taxonomyDb = taxonomyDb;
        this.virusInterpretationModel = virusInterpretationModel;
        this.virusBlacklistModel = virusBlacklistModel;
    }

    @NotNull
    public List<ReportableVirusBreakend> analyzeVirusBreakend(@NotNull List<VirusBreakend> virusBreakends) {
        List<VirusBreakend> virusBreakendsFiltered = Lists.newArrayList();
        for (VirusBreakend virusBreakend : virusBreakends) {
            if (include(virusBreakend)) {
                virusBreakendsFiltered.add(virusBreakend);
            }
        }

        List<ReportableVirusBreakend> virusBreakendsReportable = Lists.newArrayList();
        String interpretations;

        for (VirusBreakend virusBreakend : virusBreakendsFiltered) {
            String virusName = taxonomyDb.lookupName(virusBreakend.referenceTaxid());


            if (virusInterpretationModel.hasInterpretation(virusBreakend.taxidSpecies())) {
                interpretations = virusInterpretationModel.interpretVirusSpecies(virusBreakend.taxidSpecies());
            } else {
                LOGGER.info(" VIRUS breakend has called a non-interpreted virus: {}", virusName);
                interpretations = null;
            }

            virusBreakendsReportable.add(ImmutableReportableVirusBreakend.builder()
                    .virusName(virusName)
                    .integrations(virusBreakend.integrations())
                    .interpretations(interpretations)
                    .build());
        }

        return virusBreakendsReportable;
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
