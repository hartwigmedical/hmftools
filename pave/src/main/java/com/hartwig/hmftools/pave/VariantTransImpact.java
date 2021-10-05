package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.variant.impact.VariantEffect.effectsToString;
import static com.hartwig.hmftools.pave.PaveConstants.ITEM_DELIM;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;

public class VariantTransImpact
{
    public final TranscriptData TransData;

    private final List<VariantEffect> mEffects;

    private CodingContext mCodingContext;
    private ProteinContext mProteinContext;

    private boolean mInSpliceRegion;
    private int mExonRank;

    private int mLocalPhaseSetId;

    public VariantTransImpact(final TranscriptData transData)
    {
        TransData = transData;

        mEffects = Lists.newArrayList();

        mCodingContext = null;
        mExonRank = 0;
        mLocalPhaseSetId = -1;
        mInSpliceRegion = false;
    }

    public void addEffect(final VariantEffect effect)
    {
        if(effect == null)
            return;

        if(mEffects.contains(effect))
            return;

        // add in order of significance
        int index = 0;
        while(index < mEffects.size())
        {
            if(effect.rank() > mEffects.get(index).rank())
                break;

            ++index;
        }

        mEffects.add(index, effect);
    }

    public List<VariantEffect> effects() { return mEffects; }

    public VariantEffect topEffect() { return mEffects.get(0); };
    public int topRank() { return mEffects.stream().mapToInt(x -> x.rank()).max().orElse(-1); }

    public CodingContext getCodingContext() { return mCodingContext; }
    public void setCodingContext(final CodingContext context) { mCodingContext = context; }

    public ProteinContext getProteinContext() { return mProteinContext; }
    public void setProteinContext(final ProteinContext context) { mProteinContext = context; }
    public boolean hasCodingData() { return mProteinContext != null && mProteinContext.hasCodingBases(); }

    public void setExonRank(int rank) { mExonRank = rank; }
    public int exonRank() { return mExonRank; }

    public void markSpliceRegion() { mInSpliceRegion = true; }
    public boolean inSpliceRegion() { return mInSpliceRegion; }

    public String hgvsCodingChange() { return ""; }
    public String hgvsProteinChange() { return ""; }

    public String effectsStr()
    {
        return effectsToString(mEffects);
    }
    public String effectsToCsv() { return effectsToString(mEffects, ITEM_DELIM); }

    public String toString()
    {
        return String.format("trans(%s) effects(%s) inSplice(%s)",
                TransData.TransName, effectsStr(), mInSpliceRegion);
    }
}
