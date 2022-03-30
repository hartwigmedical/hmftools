package com.hartwig.hmftools.serve.sources.actin;

import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.extraction.EventExtractorFactory;
import com.hartwig.hmftools.serve.extraction.util.DriverInconsistencyMode;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;

import org.jetbrains.annotations.NotNull;

public final class ActinExtractorFactory {

    private ActinExtractorFactory() {
    }

    @NotNull
    public static ActinExtractor buildActinExtractor(@NotNull EventClassifierConfig config, @NotNull RefGenomeResource refGenomeResource) {
        return new ActinExtractor(EventExtractorFactory.create(config, refGenomeResource, DriverInconsistencyMode.WARN_ONLY));
    }
}
