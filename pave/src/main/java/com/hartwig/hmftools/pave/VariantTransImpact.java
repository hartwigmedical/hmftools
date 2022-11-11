package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.effectsToString;
import static com.hartwig.hmftools.pave.PaveConstants.DELIM;
import static com.hartwig.hmftools.pave.PaveConstants.ITEM_DELIM;

import java.util.List;
import java.util.StringJoiner;

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
    private boolean mRealigned;
    private SpliceImpactType mSpliceImpactType;
    public boolean mPhasedFrameshift;

    public VariantTransImpact(final TranscriptData transData)
    {
        TransData = transData;

        mEffects = Lists.newArrayList();

        mCodingContext = new CodingContext();
        mProteinContext = null;
        mInSpliceRegion = false;
        mPhasedFrameshift = false;
        mRealigned = false;
        mSpliceImpactType = SpliceImpactType.UNKNOWN;
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

    public VariantEffect topEffect() { return !mEffects.isEmpty() ? mEffects.get(0) : null; }; // now only used in testing
    public boolean hasEffect(final VariantEffect effect) { return mEffects.contains(effect); };
    public int topRank() { return mEffects.stream().mapToInt(x -> x.rank()).max().orElse(-1); }

    public CodingContext codingContext() { return mCodingContext; }

    public boolean isExonic()
    {
        return mCodingContext != null && mCodingContext.RegionType == EXONIC;
    }

    public ProteinContext proteinContext() { return mProteinContext; }
    public void setProteinContext(final ProteinContext context) { mProteinContext = context; }

    public boolean hasCodingBases() { return mProteinContext != null && mProteinContext.hasCodingBases(); }
    public boolean hasAminoAcids() { return mProteinContext != null && mProteinContext.hasAminoAcids(); }
    public boolean hasProteinContext() { return mProteinContext != null; }

    public void markSpliceRegion() { mInSpliceRegion = true; }
    public boolean inSpliceRegion() { return mInSpliceRegion; }

    public void markPhasedFrameshift() { mPhasedFrameshift = true; }
    public boolean phasedFrameshift() { return mPhasedFrameshift; }

    public void markRealigned() { mRealigned = true; }
    public boolean realigned() { return mRealigned; }

    public void setSpliceImpactType(SpliceImpactType type) { mSpliceImpactType = type; }
    public SpliceImpactType spliceImpactType() { return mSpliceImpactType; }

    public String hgvsCoding() { return mCodingContext.Hgvs; }
    public String hgvsProtein() { return mProteinContext != null ? mProteinContext.Hgvs  :  ""; }

    public String effectsStr()
    {
        return effectsToString(mEffects);
    }
    public String effectsToCsv() { return effectsToString(mEffects, ITEM_DELIM); }

    public static String csvHeader()
    {
        return "TransId,Canonical,IsCoding,Strand,SpliceRegion,SpliceImpact,PhasedInframe,Realigned,Effects";
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DELIM);

        sj.add(TransData.TransName);
        sj.add(String.valueOf(TransData.IsCanonical));
        sj.add(String.valueOf(!TransData.nonCoding()));
        sj.add(String.valueOf(TransData.Strand));
        sj.add(String.valueOf(mInSpliceRegion));
        sj.add(String.valueOf(mSpliceImpactType));
        sj.add(String.valueOf(mPhasedFrameshift));
        sj.add(String.valueOf(mRealigned));
        sj.add(effectsToCsv());

        return sj.toString();
    }

    public String toString()
    {
        return String.format("trans(%s) region(%s - %s) effects(%s) inSplice(%s)",
                TransData.TransName,
                mCodingContext != null ? mCodingContext.RegionType : "",
                mCodingContext != null ? mCodingContext.CodingType : "",
                effectsStr(), mInSpliceRegion);
    }
}
