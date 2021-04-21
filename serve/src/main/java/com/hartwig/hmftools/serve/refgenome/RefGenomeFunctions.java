package com.hartwig.hmftools.serve.refgenome;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ImmutableExtractionResult;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public final class RefGenomeFunctions {

    private RefGenomeFunctions() {
    }

    @NotNull
    public static List<ExtractionResult> makeRefDependent(@NotNull List<ExtractionResult> results,
            @NotNull RefGenomeVersion refGenomeVersion, @NotNull RefGenomeResource resource) {
        // TODO Implement properly
        List<ExtractionResult> refDependentResults = Lists.newArrayList();
        for (ExtractionResult result : results) {
            refDependentResults.add(ImmutableExtractionResult.builder().from(result).refGenomeVersion(refGenomeVersion).build());
        }
        return refDependentResults;
    }
}
