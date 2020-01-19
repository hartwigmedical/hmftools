package com.hartwig.hmftools.sage.phase;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.pipeline.MnvPipeline;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

public class Phase implements Consumer<SageVariant> {
    private final SnvSnvMerge mnvMerge;
    private final LocalPhaseSet localPhaseSet;
    private final DedupIndel dedupIndel;

    public Phase(@NotNull final SageConfig config, @NotNull final List<VariantHotspot> hotspots,
            @NotNull final List<GenomeRegion> panelRegions, @NotNull final MnvPipeline mnvPipeline,
            @NotNull final Consumer<SageVariant> consumer) {
        dedupIndel = new DedupIndel(consumer);
        mnvMerge = new SnvSnvMerge(config, dedupIndel, panelRegions, hotspots, mnvPipeline);
        localPhaseSet = new LocalPhaseSet(config.germlineOnly(), mnvMerge);

    }

    @Override
    public void accept(final SageVariant sageVariant) {
        localPhaseSet.accept(sageVariant);
    }

    public void flush() {
        localPhaseSet.flush();
        mnvMerge.flush();
        dedupIndel.flush();
    }
}
