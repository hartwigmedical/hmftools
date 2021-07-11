package com.hartwig.hmftools.sage.impact;

import java.util.List;

import com.google.common.collect.Lists;

public class VariantData
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;

    private final List<VariantTransImpact> mImpacts;

    public VariantData(final String chromosome, final int position, final String ref, final String alt)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;

        mImpacts = Lists.newArrayList();
    }

    public List<VariantTransImpact> getImpacts() { return mImpacts; }

    public void addImpacts(final VariantTransImpact impact)
    {
        if(mImpacts.stream().anyMatch(x -> x.TransData.TransId == impact.TransData.TransId))
            return;

        mImpacts.add(impact);
    }

    public VariantTransImpact getWorstImpact()
    {
        return null;
    }

    public VariantTransImpact getCanonicalImpact()
    {
        return mImpacts.stream().filter(x -> x.TransData.IsCanonical).findFirst().orElse(null);
    }

    public String toString()
    {
        return String.format("pos(%s:%d) variant(%s>%s)", Chromosome, Position, Ref, Alt);
    }
}
