package com.hartwig.hmftools.sage.impact;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.VariantConsequence;

public class VariantData
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;

    private final List<VariantImpact> mImpacts;

    public VariantData(final String chromosome, final int position, final String ref, final String alt)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;

        mImpacts = Lists.newArrayList();
    }

    public List<VariantImpact> getImpacts() { return mImpacts; }

    public void addImpacts(final VariantImpact impact)
    {
        if(mImpacts.stream().anyMatch(x -> x.TransData.TransId == impact.TransData.TransId))
            return;

        mImpacts.add(impact);
    }

    public VariantImpact getWorstImpact()
    {

    }

    public VariantImpact getCanonicalImpact()
    {
        return mImpacts.stream().filter(x -> x.TransData.IsCanonical).findFirst().orElse(null);
    }

    public String toString()
    {
        return String.format("pos(%s:%d) variant(%s>%s)", Chromosome, Position, Ref, Alt);
    }
}
