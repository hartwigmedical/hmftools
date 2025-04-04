package com.hartwig.hmftools.sage.tinc;

import static java.lang.Math.max;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.sage.SageConstants.LONG_INSERT_LENGTH;
import static com.hartwig.hmftools.sage.tinc.TincConstants.RECOVERY_FILTERS;

import java.util.Arrays;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.sage.filter.SoftFilter;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class VariantData
{
    public final VariantContext Context;
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;

    public final int ReferenceDepth;
    public final int ReferenceAltFrags;
    public final int TumorDepth;
    public final int TumorAltFrags;

    public final Genotype RefGenotype;
    public final Genotype TumorGenotype;

    public final boolean OnlyGermlineFiltered;

    private boolean mPonFiltered;
    public int mPonSampleCount;
    public int mPonMaxReadCount;
    public int mPonMeanReadCount;
    private Double mGnomadFrequency;

    public double mReferenceAltFragsReduction;

    private Set<SoftFilter> mNewFilters;

    public VariantData(
            final String chromosome, final int position, final String ref, final String alt,
            final int referenceDepth, final int referenceAltFrags, final int tumorDepth, final int tumorAltFrags)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        ReferenceDepth = referenceDepth;
        ReferenceAltFrags = referenceAltFrags;
        TumorDepth = tumorDepth;
        TumorAltFrags = tumorAltFrags;

        RefGenotype = null;
        TumorGenotype = null;

        Context = null;
        OnlyGermlineFiltered = false;

        mPonFiltered = false;
        mPonSampleCount = 0;
        mPonMaxReadCount = 0;
        mPonMeanReadCount = 0;
        mGnomadFrequency = null;

        mReferenceAltFragsReduction = 0;
        mNewFilters = null;
    }

    public VariantData(final VariantContext variantContext, final GenotypeIds genotypeIds)
    {
        Context = variantContext;

        Position = variantContext.getStart();
        Chromosome = variantContext.getContig();

        Ref = variantContext.getReference().getBaseString();
        Alt = !variantContext.getAlternateAlleles().isEmpty() ? variantContext.getAlternateAlleles().get(0).toString() : Ref;

        RefGenotype = Context.getGenotype(genotypeIds.ReferenceOrdinal);
        TumorGenotype = Context.getGenotype(genotypeIds.TumorOrdinal);

        ReferenceDepth = Context.getGenotype(genotypeIds.ReferenceOrdinal).getDP();
        ReferenceAltFrags = Context.getGenotype(genotypeIds.ReferenceOrdinal).getAD()[1];
        TumorDepth = Context.getGenotype(genotypeIds.TumorOrdinal).getDP();
        TumorAltFrags = Context.getGenotype(genotypeIds.TumorOrdinal).getAD()[1];

        OnlyGermlineFiltered = variantContext.getFilters().stream().allMatch(x -> RECOVERY_FILTERS.contains(x));

        mPonFiltered = false;
        mPonSampleCount = 0;
        mPonMaxReadCount = 0;
        mPonMeanReadCount = 0;
        mGnomadFrequency = null;

        mReferenceAltFragsReduction = 0;
        mNewFilters = null;
    }

    public boolean isSnv() { return Ref.length() == 1 && Alt.length() == 1; }
    public boolean isMnv() { return Ref.length() > 1 && Alt.equals(Ref); }
    public boolean isIndel() { return Alt.length() != Ref.length(); }
    public boolean isInsert() { return Alt.length() > Ref.length(); }

    public boolean isLongIndel() { return max(Ref.length(), Alt.length()) >= LONG_INSERT_LENGTH; }

    public VariantTier tier() { return VariantTier.fromContext(Context); }

    public boolean isPassing()
    {
        if(Context == null)
            return true; // loaded from file, ignore existing filter status

        return !Context.isFiltered() || Context.getFilters().size() == 1 && Context.getFilters().contains(PASS);
    }

    public void setPonFiltered() { mPonFiltered = true;}
    public boolean ponFiltered() { return mPonFiltered;}

    public int ponSampleCount() { return mPonSampleCount; }
    public int ponMaxReadCount() { return mPonMaxReadCount; }
    public int ponMeanReadCount() { return mPonMeanReadCount; }

    public void setPonFrequency(int sampleCount, int maxReadCount, int meanReadCount)
    {
        mPonSampleCount = sampleCount;
        mPonMaxReadCount = maxReadCount;
        mPonMeanReadCount = meanReadCount;
    }

    public Double gnomadFrequency() { return mGnomadFrequency; }
    public void setGnomadFrequency(Double frequency) { mGnomadFrequency = frequency; }

    public boolean recovered() { return mNewFilters != null && mNewFilters.isEmpty(); }

    public double tumorAf() { return TumorDepth > 0 ? TumorAltFrags / (double)TumorDepth : 0; }
    public double referenceAf() { return ReferenceDepth > 0 ? ReferenceAltFrags / (double)ReferenceDepth : 0; }

    public void setReferenceAltFragReduction(double perc) { mReferenceAltFragsReduction = perc; }

    public int calcReducedAltCount(int altCount) { return (int)round(calcReducedAltValue(altCount)); }

    public double calcReducedAltValue(double altValue)
    {
        double referenceAf = referenceAf();
        if(referenceAf <= 0)
            return 0;

        double reducedAf = max(referenceAf - mReferenceAltFragsReduction, 0.0);
        double reducedRefAfPercent = reducedAf / referenceAf;
        return altValue * reducedRefAfPercent;
    }

    public void buildNewFilters()
    {
        mNewFilters = Sets.newHashSet();

        for(String filterStr : Context.getFilters())
        {
            SoftFilter filter = Arrays.stream(SoftFilter.values()).filter(x -> x.filterName().equals(filterStr)).findFirst().orElse(null);
            if(filter != null)
                mNewFilters.add(filter);
        }
    }

    public Set<SoftFilter> newFilters() { return mNewFilters; }

    public String toString()
    {
        return format("%s:%d %s>%s ref(%d/%d) tumor(%d/%d)",
                Chromosome, Position, Ref, Alt, ReferenceAltFrags, ReferenceDepth, TumorAltFrags, TumorDepth);
    }

}
