package com.hartwig.hmftools.sage.phase;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class Phase implements Consumer<SageVariant> {
    private final SnvSnvMerge mnvMerge;
    private final DedupSnv snvIndelMerge;
    private final LocalPhaseSet localPhaseSet;
    private final DedupIndel dedupIndel;

    public Phase(
            @NotNull final SageConfig config,
            @NotNull final IndexedFastaSequenceFile reference,
            @NotNull final SageVariantFactory sageVariantFactory,
            @NotNull final List<VariantHotspot> hotspots,
            @NotNull final List<GenomeRegion> panelRegions,

            @NotNull final Consumer<SageVariant> consumer) {
        final MnvFactory mnvFactory = new MnvFactory(config, sageVariantFactory, reference, hotspots, panelRegions);
        dedupIndel = new DedupIndel(consumer);
        snvIndelMerge = new DedupSnv(dedupIndel);
        mnvMerge = new SnvSnvMerge(config, snvIndelMerge, mnvFactory, panelRegions, hotspots);
        localPhaseSet = new LocalPhaseSet(config.germlineOnly(), mnvMerge);

    }

    @Override
    public void accept(final SageVariant sageVariant) {
        localPhaseSet.accept(sageVariant);
    }

    public void flush() {
        localPhaseSet.flush();
        mnvMerge.flush();
        snvIndelMerge.flush();
        dedupIndel.flush();
    }
}
