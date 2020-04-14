package com.hartwig.hmftools.sage.phase;

import java.util.function.Consumer;

import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

public class Phase implements Consumer<SageVariant> {
    private final DedupMnv mnvMerge;
    private final LocalPhaseSet localPhaseSet;
    private final DedupIndel dedupIndel;
    private final MixedGermlineMnv mixedGermlineMnv;
    private final MixedGermlineImpact mixedGermlineImpact;

    public Phase(@NotNull final SageConfig config, @NotNull final Consumer<SageVariant> consumer) {
        dedupIndel = new DedupIndel(consumer);
        mnvMerge = new DedupMnv(dedupIndel);
        mixedGermlineImpact = new MixedGermlineImpact(mnvMerge);
        mixedGermlineMnv = new MixedGermlineMnv(mixedGermlineImpact);
        localPhaseSet = new LocalPhaseSet(config.readContextFlankSize(), mixedGermlineMnv);
    }

    @Override
    public void accept(final SageVariant sageVariant) {
        localPhaseSet.accept(sageVariant);
    }

    public void flush() {
        localPhaseSet.flush();
        mixedGermlineMnv.flush();
        mixedGermlineImpact.flush();
        mnvMerge.flush();
        dedupIndel.flush();
    }
}
