package com.hartwig.hmftools.serve.sources.ckb;

import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.extraction.EventExtractorFactory;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;
import com.hartwig.hmftools.serve.sources.ckb.curation.DrugCurator;
import com.hartwig.hmftools.serve.sources.ckb.curation.EvidenceLevelCurator;

import org.jetbrains.annotations.NotNull;

public class CkbExtractorFactory {

    private CkbExtractorFactory() {
    }

    @NotNull
    public static CkbExtractor buildCkbExtractor(@NotNull EventClassifierConfig config, @NotNull RefGenomeResource refGenomeResource,
            @NotNull DoidLookup missingDoidLookup) {
        return new CkbExtractor(EventExtractorFactory.create(config, refGenomeResource),
                new ActionableEvidenceFactory(missingDoidLookup, new DrugCurator(), new EvidenceLevelCurator()));
    }
}
