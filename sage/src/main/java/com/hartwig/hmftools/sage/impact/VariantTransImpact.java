package com.hartwig.hmftools.sage.impact;

import static com.hartwig.hmftools.common.variant.VariantConsequence.STOP_GAINED;
import static com.hartwig.hmftools.common.variant.VariantConsequence.STOP_LOST;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.ITEM_DELIM;

import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;

public class VariantTransImpact
{
    public final TranscriptData TransData;

    private final List<String> mConsequenceEffects;

    private String mUpstreamAA;
    private String mWildtypeAA;
    private String mNovelAA;
    private String mDownstreamAA;

    private int mCodingBases;
    private int mBasesToLastExonJunction;
    private boolean mInSpliceRegion;

    private int mLocalPhaseSetId;

    public VariantTransImpact(final TranscriptData transData, final VariantConsequence consequence)
    {
        this(transData, consequence.parentTerm());
    }

    public VariantTransImpact(final TranscriptData transData, final String consequenceEffect)
    {
        TransData = transData;

        mConsequenceEffects = Lists.newArrayList(consequenceEffect);

        mUpstreamAA = "";
        mWildtypeAA = "";
        mNovelAA = "";
        mDownstreamAA = "";
        mCodingBases = 0;
        mBasesToLastExonJunction = 0;
        mLocalPhaseSetId = -1;
        mInSpliceRegion = false;
    }

    public void addConsequence(final String consequenceEffect)
    {
        VariantConsequence consequence = VariantConsequence.fromEffect(consequenceEffect);

        List<VariantConsequence> consequences = consequences();

        if(consequences.contains(consequence))
            return;

        int index = 0;
        while(index < consequences.size())
        {
            if(consequence.rank() > consequences.get(index).rank())
                break;

            ++index;
        }

        mConsequenceEffects.add(index, consequenceEffect);
    }

    public VariantConsequence consequence() { return VariantConsequence.fromEffect(mConsequenceEffects.get(0)); };

    public List<String> consequenceEffects() { return mConsequenceEffects; }

    public List<VariantConsequence> consequences()
    {
        return mConsequenceEffects.stream().map(x -> VariantConsequence.fromEffect(x)).collect(Collectors.toList());
    }

    public int topRank() { return consequences().stream().mapToInt(x -> x.rank()).max().orElse(-1); }

    public String upstreamAA() { return mUpstreamAA; }
    public String wildtypeAA() { return mWildtypeAA; }
    public String novelAA() { return mNovelAA; }
    public String downstreamAA() { return mDownstreamAA; }

    public void setAminoAcids(final String up, final String wildtype, final String novel, final String down)
    {
        mUpstreamAA = up;
        mWildtypeAA = wildtype;
        mNovelAA = novel;
        mDownstreamAA = down;
    }

    public void setCodingBases(int codingBase, int toLastExonJunc)
    {
        mCodingBases = codingBase;
        mBasesToLastExonJunction = toLastExonJunc;
    }

    public int codingBases() { return mCodingBases; }
    public int basesToLastExonJunction() { return mBasesToLastExonJunction; }

    public void markSpliceRegion() { mInSpliceRegion = true; }
    public boolean inSpliceRegion() { return mInSpliceRegion; }

    public String hgvsCodingChange() { return ""; }
    public String hgvsProteinChange() { return ""; }

    public SnpEffAnnotation findMatchingAnnotation(final List<SnpEffAnnotation> annotations)
    {
        return annotations.stream()
                .filter(x -> x.featureID() != null && x.featureID().equals(TransData.TransName))
                .findFirst().orElse(null);
    }

    public String consequencesStr()
    {
        StringJoiner sj = new StringJoiner(ITEM_DELIM);
        consequences().forEach(x -> sj.add(x.toString()));
        return sj.toString();
    }

    public String effectsStr()
    {
        StringJoiner sj = new StringJoiner(ITEM_DELIM);
        mConsequenceEffects.forEach(x -> sj.add(x.toString()));
        return sj.toString();
    }

    public String toString()
    {
        return String.format("trans(%s) conseq(%s) effects(%s) inSplice(%s)",
                TransData.TransName, consequencesStr(), effectsStr(), mInSpliceRegion);
    }
}
