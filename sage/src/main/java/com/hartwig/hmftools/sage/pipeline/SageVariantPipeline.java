package com.hartwig.hmftools.sage.pipeline;

import java.util.List;
import java.util.concurrent.CompletableFuture;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

public interface SageVariantPipeline
{
    @NotNull
    CompletableFuture<List<SageVariant>> variants(@NotNull final BaseRegion region);
}
