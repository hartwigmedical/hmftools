package com.hartwig.hmftools.serve.sources.vicc;

import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.extraction.EventExtractorFactory;
import com.hartwig.hmftools.serve.extraction.catalog.DriverInconsistencyMode;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;
import com.hartwig.hmftools.serve.sources.vicc.curation.DrugCurator;
import com.hartwig.hmftools.serve.sources.vicc.curation.EvidenceLevelCurator;

import org.jetbrains.annotations.NotNull;

public final class ViccExtractorFactory {
    
    private ViccExtractorFactory() {
    }

    @NotNull
    public static ViccExtractor buildViccExtractor(@NotNull EventClassifierConfig config, @NotNull RefGenomeResource refGenomeResource,
            @NotNull DoidLookup missingDoidLookup) {
        return new ViccExtractor(EventExtractorFactory.create(config, refGenomeResource, DriverInconsistencyMode.FILTER),
                new ActionableEvidenceFactory(missingDoidLookup, new DrugCurator(), new EvidenceLevelCurator()));
    }
}
