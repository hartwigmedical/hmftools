package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.variant.VariantConsequence.consequencesToString;
import static com.hartwig.hmftools.pave.PaveConstants.ITEM_DELIM;

import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.variant.VariantConsequence;

public class VariantTransImpact
{
    public final TranscriptData TransData;

    private final List<String> mConsequenceEffects;
    private final List<VariantConsequence> mConsequences;

    private CodingContext mCodingContext;
    private ProteinContext mProteinContext;

    private boolean mInSpliceRegion;
    private int mExonRank;

    private int mLocalPhaseSetId;

    @Deprecated
    public VariantTransImpact(final TranscriptData transData, final String consequenceEffect)
    {
        TransData = transData;

        mConsequenceEffects = Lists.newArrayList(consequenceEffect);
        mConsequences = Lists.newArrayList();

        mCodingContext = null;
        mExonRank = 0;
        mLocalPhaseSetId = -1;
        mInSpliceRegion = false;
    }

    public VariantTransImpact(final TranscriptData transData)
    {
        TransData = transData;

        mConsequenceEffects = Lists.newArrayList();
        mConsequences = Lists.newArrayList();

        mCodingContext = null;
        mExonRank = 0;
        mLocalPhaseSetId = -1;
        mInSpliceRegion = false;
    }

    public void addConsequence(final String consequenceEffect)
    {
        if(consequenceEffect == null)
            return;

        VariantConsequence consequence = VariantConsequence.fromEffect(consequenceEffect);
        addConsequence(consequence, consequenceEffect);
    }

    public void addConsequence(final VariantConsequence consequence)
    {
        addConsequence(consequence, consequence.parentTerm());
    }

    public void addConsequence(final VariantConsequence consequence, final String consequenceEffect)
    {
        if(consequence == null || consequenceEffect == null)
            return;

        if(mConsequences.contains(consequence))
            return;

        // add in order of significance
        int index = 0;
        while(index < mConsequences.size())
        {
            if(consequence.rank() > mConsequences.get(index).rank())
                break;

            ++index;
        }

        mConsequences.add(index, consequence);
        mConsequenceEffects.add(index, consequenceEffect);
    }

    public VariantConsequence topConsequence() { return mConsequences.get(0); };

    public List<String> consequenceEffects() { return mConsequenceEffects; }

    public List<VariantConsequence> consequences()
    {
        return mConsequences;
    }

    public int topRank() { return consequences().stream().mapToInt(x -> x.rank()).max().orElse(-1); }

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
        return consequencesToString(consequences());
    }

    public String rawConsequencesStr()
    {
        StringJoiner sj = new StringJoiner(ITEM_DELIM);
        mConsequenceEffects.forEach(x -> sj.add(x));
        return sj.toString();
    }

    public String toString()
    {
        return String.format("trans(%s) conseq(%s) effects(%s) inSplice(%s)",
                TransData.TransName, rawConsequencesStr(), effectsStr(), mInSpliceRegion);
    }
}
