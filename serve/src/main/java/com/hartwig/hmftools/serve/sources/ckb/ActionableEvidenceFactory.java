package com.hartwig.hmftools.serve.sources.ckb;

import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.sources.vicc.curation.DrugCurator;
import com.hartwig.hmftools.serve.sources.vicc.curation.EvidenceLevelCurator;

import org.jetbrains.annotations.NotNull;

class ActionableEvidenceFactory {

    @NotNull
    private final DoidLookup missingDoidLookup;
    @NotNull
    private final DrugCurator drugCurator;
    @NotNull
    private final EvidenceLevelCurator evidenceLevelCurator;

    public ActionableEvidenceFactory(@NotNull final DoidLookup missingDoidLookup, @NotNull final DrugCurator drugCurator,
            @NotNull final EvidenceLevelCurator evidenceLevelCurator) {
        this.missingDoidLookup = missingDoidLookup;
        this.drugCurator = drugCurator;
        this.evidenceLevelCurator = evidenceLevelCurator;
    }
}
