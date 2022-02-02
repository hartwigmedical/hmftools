package com.hartwig.hmftools.sage.phase;

import java.util.List;
import java.util.Set;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.sage.common.SageVariant;

public class VariantDeduper implements Consumer<SageVariant>
{
    public static final int PHASE_BUFFER = 150;

    private final DedupRealign mDedupRealign;
    private final DedupMnv mDedupMnv;
    private final LocalPhaseSet mLocalPhaseSet;
    private final LocalRealignSet mLocalRealignSet;
    private final DedupIndel mDedupIndel;
    private final MixedSomaticGermlineIdentifier mMixedSomaticGermlineIdentifier;
    private final MixedSomaticGermlineDedup mMixedSomaticGermlineDedup;

    public VariantDeduper(final List<TranscriptData> transcripts, final PhaseSetCounter phaseSetCounter, final Consumer<SageVariant> consumer)
    {
        mDedupRealign = new DedupRealign(consumer);
        mDedupIndel = new DedupIndel(mDedupRealign);
        mDedupMnv = new DedupMnv(mDedupIndel);
        mMixedSomaticGermlineDedup = new MixedSomaticGermlineDedup(mDedupMnv, transcripts);
        mMixedSomaticGermlineIdentifier = new MixedSomaticGermlineIdentifier(mMixedSomaticGermlineDedup);
        mLocalRealignSet = new LocalRealignSet(mMixedSomaticGermlineIdentifier);
        mLocalPhaseSet = new LocalPhaseSet(phaseSetCounter, mLocalRealignSet);
    }

    public Set<Integer> passingPhaseSets()
    {
        return mLocalPhaseSet.passingPhaseSets();
    }

    @Override
    public void accept(final SageVariant sageVariant)
    {
        mLocalPhaseSet.accept(sageVariant);
    }

    public void flush()
    {
        mLocalPhaseSet.flush();
        mLocalRealignSet.flush();
        mMixedSomaticGermlineIdentifier.flush();
        mMixedSomaticGermlineDedup.flush();
        mDedupMnv.flush();
        mDedupIndel.flush();
        mDedupRealign.flush();
    }
}
