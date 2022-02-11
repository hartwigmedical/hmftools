package com.hartwig.hmftools.sage.phase;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.sage.common.SageVariant;

public class VariantDeduperOld implements Consumer<SageVariant>
{
    public static final int PHASE_BUFFER = 150;

    private final DedupRealign mDedupRealign;
    private final DedupMnvOld mDedupMnv;
    private final LocalRealignSet mLocalRealignSet;
    private final DedupIndelOld mDedupIndel;
    private final MixedSomaticGermlineIdentifier mMixedSomaticGermlineIdentifier;
    private final MixedSomaticGermlineDedup mMixedSomaticGermlineDedup;

    public VariantDeduperOld(final List<TranscriptData> transcripts, final Consumer<SageVariant> consumer)
    {
        mDedupRealign = new DedupRealign(consumer);
        mDedupIndel = new DedupIndelOld(mDedupRealign);
        mDedupMnv = new DedupMnvOld(mDedupIndel);
        mMixedSomaticGermlineDedup = new MixedSomaticGermlineDedup(mDedupMnv, transcripts);
        mMixedSomaticGermlineIdentifier = new MixedSomaticGermlineIdentifier(mMixedSomaticGermlineDedup);
        mLocalRealignSet = new LocalRealignSet(mMixedSomaticGermlineIdentifier);
    }

    @Override
    public void accept(final SageVariant sageVariant)
    {
        mLocalRealignSet.accept(sageVariant);
    }

    public void flush()
    {
        mLocalRealignSet.flush();
        mMixedSomaticGermlineIdentifier.flush();
        mMixedSomaticGermlineDedup.flush();
        mDedupMnv.flush();
        mDedupIndel.flush();
        mDedupRealign.flush();
    }
}
