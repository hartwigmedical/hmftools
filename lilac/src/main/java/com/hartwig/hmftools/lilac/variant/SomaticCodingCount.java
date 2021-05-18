package com.hartwig.hmftools.lilac.variant;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import java.util.List;
import java.util.stream.Collectors;

public class SomaticCodingCount
{
    public final HlaAllele Allele;

    private double mInframeIndel;
    private double mMissense;
    private double mNonsense;
    private double mSplice;
    private double mSynonymous;

    public SomaticCodingCount(final HlaAllele allele, double inframeIndel, double missense, double nonsense, double splice, double synonymous)
    {
        Allele = allele;
        mInframeIndel = inframeIndel;
        mMissense = missense;
        mNonsense = nonsense;
        mSplice = splice;
        mSynonymous = synonymous;
    }

    public double inframeIndel() { return mInframeIndel; }
    public double missense() { return mMissense; }
    public double nonsense() { return mNonsense; }
    public double splice() { return mSplice; }
    public double synonymous() { return mSynonymous; }

    public String toString()
    {
        return "SomaticCodingCount(allele=" + Allele + ", inframeIndel=" + mInframeIndel + ", missense=" + mMissense
                + ", nonsense=" + mNonsense + ", splice=" + mSplice + ", synonymous=" + mSynonymous + ")";
    }

    public double total() { return mInframeIndel + mMissense + mNonsense + mSplice + mSynonymous; }

    public static List<SomaticCodingCount> create(final List<HlaAllele> winners)
    {
        return winners.stream()
                .map(x -> new SomaticCodingCount(x, 0, 0, 0, 0, 0)).collect(Collectors.toList());
    }

    public static void addVariant(
            final List<SomaticCodingCount> codingCounts, final VariantContextDecorator variant, final List<HlaAllele> variantAlleles)
    {
        boolean isIndel = variant.alt().length() != variant.ref().length();
        CodingEffect codingEffect = variant.canonicalCodingEffect();
        addVariant(codingCounts, isIndel, codingEffect, variantAlleles);
    }

    public static void addVariant(
            final List<SomaticCodingCount> codingCounts, boolean isIndel, final CodingEffect effect, final List<HlaAllele> variantAlleles)
    {
        double contribution = 1.0 / variantAlleles.size();

        for(HlaAllele allele : variantAlleles)
        {
            codingCounts.stream().filter(x -> x.Allele.equals(allele)).forEach(x -> x.addVariant(isIndel, effect, contribution));
        }
        /*
            // Alternative, give only to first to first
            if (counts.isNotEmpty()) {
                result.add(counts[0].addVariant(indel, effect, contribution))
                result.addAll(counts.takeLast(counts.size - 1))
            }
        */
    }

    private void addVariant(boolean indel, CodingEffect effect, double contribution)
    {
        if(indel && effect == CodingEffect.MISSENSE)
        {
            mInframeIndel += contribution;
        }

        switch(effect)
        {
            case MISSENSE:
                mMissense += contribution;
                break;

            case NONSENSE_OR_FRAMESHIFT:
                mNonsense += contribution;
                break;

            case SPLICE:
                mSplice += contribution;
                break;

            case SYNONYMOUS:
                mSynonymous += contribution;
                break;

            default:
                break;
        }
    }

}
