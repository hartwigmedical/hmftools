package com.hartwig.hmftools.virusinterpreter;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.ImmutableAnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusBreakend;
import com.hartwig.hmftools.common.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.common.virus.VirusInterpretation;
import com.hartwig.hmftools.virusinterpreter.algo.VirusBlacklistModel;
import com.hartwig.hmftools.virusinterpreter.algo.VirusInterpretationModel;
import com.hartwig.hmftools.virusinterpreter.taxonomy.TaxonomyDb;

import org.jetbrains.annotations.NotNull;

public class VirusInterpreterAlgo {

    @NotNull
    private final TaxonomyDb taxonomyDb;
    @NotNull
    private final VirusInterpretationModel virusInterpretationModel;
    @NotNull
    private final VirusBlacklistModel virusBlacklistModel;

    public VirusInterpreterAlgo(@NotNull final TaxonomyDb taxonomyDb,
            @NotNull final VirusInterpretationModel virusInterpretationModel, @NotNull final VirusBlacklistModel virusBlacklistModel) {
        this.taxonomyDb = taxonomyDb;
        this.virusInterpretationModel = virusInterpretationModel;
        this.virusBlacklistModel = virusBlacklistModel;
    }

    @NotNull
    public List<AnnotatedVirus> analyze(@NotNull List<VirusBreakend> virusBreakends) {
        List<AnnotatedVirus> annotatedViruses = Lists.newArrayList();
        for (VirusBreakend virusBreakend : virusBreakends) {
            VirusInterpretation interpretation = null;
            if (virusInterpretationModel.hasInterpretation(virusBreakend.taxidSpecies())) {
                interpretation = virusInterpretationModel.interpretVirusSpecies(virusBreakend.taxidSpecies());
            }

            int taxid = virusBreakend.referenceTaxid();
            annotatedViruses.add(ImmutableAnnotatedVirus.builder()
                    .taxid(taxid)
                    .name(taxonomyDb.lookupName(taxid))
                    .qcStatus(virusBreakend.qcStatus())
                    .integrations(virusBreakend.integrations())
                    .interpretation(interpretation)
                    .reported(report(virusBreakend))
                    .build());
        }

        return annotatedViruses;
    }

    private boolean report(@NotNull VirusBreakend virusBreakend) {
        if (virusBreakend.qcStatus() == VirusBreakendQCStatus.LOW_VIRAL_COVERAGE) {
            return false;
        }

        if (virusBreakend.integrations() == 0) {
            return false;
        }

        return !virusBlacklistModel.isBlacklisted(virusBreakend);
    }
}
