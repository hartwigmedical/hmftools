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
        this(transData, consequence.description());
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

        if(consequences().contains(consequence))
            return;

        // apply any ordering / prioritisation
        if(consequence == STOP_GAINED || consequence == STOP_LOST)
            mConsequenceEffects.add(0, consequenceEffect);
        else
            mConsequenceEffects.add(consequenceEffect);
    }

    public VariantConsequence consequence() { return VariantConsequence.fromEffect(mConsequenceEffects.get(0)); };

    public List<VariantConsequence> consequences()
    {
        return mConsequenceEffects.stream().map(x -> VariantConsequence.fromEffect(x)).collect(Collectors.toList());
    }

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

    public String codingChange() { return ""; }
    public String proteinChange() { return ""; }

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
}
