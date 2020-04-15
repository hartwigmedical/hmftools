package com.hartwig.hmftools.sage.phase;

import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

public class Phase implements Consumer<SageVariant> {
    private final DedupMnv dedupMnv;
    private final LocalPhaseSet localPhaseSet;
    private final DedupIndel dedupIndel;
    private final MixedGermlineMnv mixedGermlineMnv;
    private final MixedGermlineDedup mixedGermlineDedup;
    private final PhasedInframeIndel phasedInframeIndel;
    private final RightAlignMicrohomology rightAlignMicrohomology;

    public Phase(@NotNull final SageConfig config, @NotNull final String chromosome, @NotNull final Consumer<SageVariant> consumer) {

        final List<HmfTranscriptRegion> transcripts =
                config.transcriptRegions().stream().filter(x -> x.chromosome().equals(chromosome)).collect(Collectors.toList());

        dedupIndel = new DedupIndel(consumer);
        dedupMnv = new DedupMnv(dedupIndel);
        mixedGermlineDedup = new MixedGermlineDedup(dedupMnv, transcripts);
        mixedGermlineMnv = new MixedGermlineMnv(mixedGermlineDedup);
        phasedInframeIndel = new PhasedInframeIndel(mixedGermlineMnv, transcripts);
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
        mixedGermlineDedup.flush();
        dedupMnv.flush();
        dedupIndel.flush();
    }
}
