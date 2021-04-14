package com.hartwig.hmftools.serve.sources.ckb;

import java.util.List;

import com.hartwig.hmftools.common.refseq.RefSeq;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.extraction.EventExtractorFactory;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;

import org.jetbrains.annotations.NotNull;

public class CkbExtractorFactory {

    private CkbExtractorFactory() {
    }

    @NotNull
    public static CkbExtractor buildCkbExtractor(@NotNull EventClassifierConfig config, @NotNull RefGenomeResource refGenomeResource,
            @NotNull List<RefSeq> refSeqMappings) {
        return new CkbExtractor(EventExtractorFactory.create(config, refGenomeResource), new ActionableEvidenceFactory(), refSeqMappings);
    }
}
