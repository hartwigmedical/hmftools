package com.hartwig.hmftools.ctdna.purity;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;

public class SomaticVariant
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;
    public final VariantTier Tier;
    public final VariantType Type;
    public final int RepeatCount;
    public final double Mappability;
    public final double SubclonalPerc;

    public final List<GenotypeFragments> Samples;
    public final boolean PassFilters;

    public SomaticVariant(
            final String chromosome, final int position, final String ref, final String alt, final VariantTier tier,
            final VariantType type, final int repeatCount, final double mappability, final double subclonalPerc, final boolean passFilters)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Tier = tier;
        Type = type;
        RepeatCount = repeatCount;
        Mappability = mappability;
        SubclonalPerc = subclonalPerc;
        Samples = Lists.newArrayList();
        PassFilters = passFilters;
    }

    public GenotypeFragments findGenotypeData(final String sampleId)
    {
        return Samples.stream().filter(x -> x.SampleName.equals(sampleId)).findFirst().orElse(null);
    }

    public String toString() { return format("%s:%d %s>%s", Chromosome, Position, Ref, Alt); }
}
