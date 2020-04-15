package com.hartwig.hmftools.sage.phase;

import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

public class Phase implements Consumer<SageVariant> {
    private final DedupMnv mnvMerge;
    private final LocalPhaseSet localPhaseSet;
    private final DedupIndel dedupIndel;
    private final MixedGermlineMnv mixedGermlineMnv;
    private final MixedGermlineImpact mixedGermlineImpact;
    private final PhasedInframeIndel phasedInframeIndel;
    private final RightAlignMicrohomology rightAlignMicrohomology;

    public Phase(@NotNull final SageConfig config, @NotNull final String chromosome, @NotNull final Consumer<SageVariant> consumer) {

        final List<HmfTranscriptRegion> transcripts =
                config.transcriptRegions().stream().filter(x -> x.chromosome().equals(chromosome)).collect(Collectors.toList());

        dedupIndel = new DedupIndel(consumer);
        mnvMerge = new DedupMnv(dedupIndel);
        mixedGermlineImpact = new MixedGermlineImpact(mnvMerge);
        mixedGermlineMnv = new MixedGermlineMnv(mixedGermlineImpact);
        phasedInframeIndel = new PhasedInframeIndel(mixedGermlineImpact, transcripts);
        rightAlignMicrohomology = new RightAlignMicrohomology(phasedInframeIndel, transcripts);
        localPhaseSet = new LocalPhaseSet(config.readContextFlankSize(), rightAlignMicrohomology);
    }

    @Override
    public void accept(final SageVariant sageVariant) {
        localPhaseSet.accept(sageVariant);
    }

    public void flush() {
        localPhaseSet.flush();
        rightAlignMicrohomology.flush();
        phasedInframeIndel.flush();
        mixedGermlineMnv.flush();
        mixedGermlineImpact.flush();
        mnvMerge.flush();
        dedupIndel.flush();
    }
}
