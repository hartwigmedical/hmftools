package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
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

    public boolean mPhasedFrameshift;

    public VariantTransImpact(final TranscriptData transData)
    {
        TransData = transData;

        mEffects = Lists.newArrayList();

        mCodingContext = null;
        mInSpliceRegion = false;
        mPhasedFrameshift = false;
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

    public CodingContext codingContext() { return mCodingContext; }
    public void setCodingContext(final CodingContext context) { mCodingContext = context; }

    public boolean isExonic()
    {
        return mCodingContext != null && mCodingContext.RegionType == EXONIC;
    }

    public ProteinContext proteinContext() { return mProteinContext; }
    public void setProteinContext(final ProteinContext context) { mProteinContext = context; }

    public boolean hasAminoAcids() { return mProteinContext != null && mProteinContext.hasCodingBases(); }
    public boolean hasProteinContext() { return mProteinContext != null; }

    public void markSpliceRegion() { mInSpliceRegion = true; }
    public boolean inSpliceRegion() { return mInSpliceRegion; }

    public void markPhasedFrameshift() { mPhasedFrameshift = true; }
    public boolean phasedFrameshift() { return mPhasedFrameshift; }

    public String hgvsCoding() { return mCodingContext.hgvsStr(); }
    public String hgvsProtein() { return mProteinContext != null ? mProteinContext.hgvsStr() :  ""; }

    public String effectsStr()
    {
        return effectsToString(mEffects);
    }
    public String effectsToCsv() { return effectsToString(mEffects, ITEM_DELIM); }

    public static String csvHeader()
    {
        return "TransId,Canonical,IsCoding,Strand,SpliceRegion,Effects";
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DELIM);

        sj.add(TransData.TransName);
        sj.add(String.valueOf(TransData.IsCanonical));
        sj.add(String.valueOf(!TransData.nonCoding()));
        sj.add(String.valueOf(TransData.Strand));
        sj.add(String.valueOf(mInSpliceRegion));
        sj.add(effectsToCsv());

        return sj.toString();
    }

    public String toString()
    {
        return String.format("trans(%s) region(%s - %s) effects(%s) inSplice(%s)",
                TransData.TransName, mCodingContext.RegionType, mCodingContext.CodingType, effectsStr(), mInSpliceRegion);
    }
}
