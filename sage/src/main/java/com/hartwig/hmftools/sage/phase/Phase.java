package com.hartwig.hmftools.sage.phase;

import java.util.List;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

public class Phase implements Consumer<SageVariant>
{

    static final int PHASE_BUFFER = 150;

    private final DedupRealign dedupRealign;
    private final DedupMnv dedupMnv;
    private final LocalPhaseSet localPhaseSet;
    private final LocalRealignSet localRealignSet;
    private final DedupIndel dedupIndel;
    private final MixedSomaticGermlineIdentifier mixedSomaticGermlineIdentifier;
    private final MixedSomaticGermlineDedup mixedSomaticGermlineDedup;
    private final PhasedInframeIndel phasedInframeIndel;
    private final RightAlignMicrohomology rightAlignMicrohomology;

    public Phase(@NotNull final SageConfig config, @NotNull final String chromosome, @NotNull final Consumer<SageVariant> consumer)
    {
        final List<HmfTranscriptRegion> transcripts =
                config.transcriptRegions().stream().filter(x -> x.chromosome().equals(chromosome)).collect(Collectors.toList());

        dedupRealign = new DedupRealign(consumer);
        dedupIndel = new DedupIndel(dedupRealign);
        dedupMnv = new DedupMnv(dedupIndel);
        mixedSomaticGermlineDedup = new MixedSomaticGermlineDedup(dedupMnv, transcripts);
        mixedSomaticGermlineIdentifier = new MixedSomaticGermlineIdentifier(mixedSomaticGermlineDedup);
        phasedInframeIndel = new PhasedInframeIndel(mixedSomaticGermlineIdentifier, transcripts);
        rightAlignMicrohomology = new RightAlignMicrohomology(phasedInframeIndel, transcripts);
        localRealignSet = new LocalRealignSet(rightAlignMicrohomology);
        localPhaseSet = new LocalPhaseSet(localRealignSet);
    }

    @NotNull
    public Set<Integer> passingPhaseSets()
    {
        return localPhaseSet.passingPhaseSets();
    }

    @Override
    public void accept(final SageVariant sageVariant)
    {
        localPhaseSet.accept(sageVariant);
    }

    public void flush()
    {
        localPhaseSet.flush();
        localRealignSet.flush();
        rightAlignMicrohomology.flush();
        phasedInframeIndel.flush();
        mixedSomaticGermlineIdentifier.flush();
        mixedSomaticGermlineDedup.flush();
        dedupMnv.flush();
        dedupIndel.flush();
        dedupRealign.flush();
    }
}
