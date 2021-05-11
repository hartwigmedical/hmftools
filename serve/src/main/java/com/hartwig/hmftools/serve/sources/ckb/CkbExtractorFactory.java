package com.hartwig.hmftools.serve.sources.ckb;

import java.util.List;

import com.hartwig.hmftools.common.refseq.RefSeq;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.extraction.EventExtractorFactory;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;

import org.jetbrains.annotations.NotNull;

public final class CkbExtractorFactory {

    // For CKB we want to ignore driver inconsistencies for now (their gene panel is larger than hmf driver panel to start with)
    private static final boolean REPORT_DRIVER_INCONSISTENCIES = false;

    private CkbExtractorFactory() {
    }

    @NotNull
    public static CkbExtractor buildCkbExtractor(@NotNull EventClassifierConfig config, @NotNull RefGenomeResource refGenomeResource,
            @NotNull List<RefSeq> refSeqMappings) {
        return new CkbExtractor(EventExtractorFactory.create(config, refGenomeResource, REPORT_DRIVER_INCONSISTENCIES),
                new ActionableEntryFactory(),
                refSeqMappings);
    }
}
