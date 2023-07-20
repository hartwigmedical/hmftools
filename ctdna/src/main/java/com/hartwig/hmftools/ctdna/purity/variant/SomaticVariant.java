package com.hartwig.hmftools.ctdna.purity.variant;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;

public class SomaticVariant
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;
    public final VariantType Type;
    public final double SubclonalPerc;

    public final List<GenotypeFragments> Samples;
    public final boolean PassFilters;

    private final VariantContextDecorator mVariant;
    private double mSequenceVcRatio;

    public SomaticVariant(final VariantContextDecorator variant, final double subclonalPerc, final boolean passFilters)
    {
        Chromosome = variant.chromosome();
        Position = variant.position();
        Ref = variant.ref();
        Alt = variant.alt();
        Type = variant.type();
        SubclonalPerc = subclonalPerc;
        Samples = Lists.newArrayList();
        PassFilters = passFilters;
        mVariant = variant;
        mSequenceVcRatio = 0;
    }

    public GenotypeFragments findGenotypeData(final String sampleId)
    {
        return Samples.stream().filter(x -> x.SampleName.equals(sampleId)).findFirst().orElse(null);
    }

    public void setSequenceGcRatio(double ratio) { mSequenceVcRatio = ratio; }
    public double sequenceGcRatio() { return mSequenceVcRatio; }

    public VariantContextDecorator decorator() { return mVariant; }

    public String toString() { return format("%s:%d %s>%s", Chromosome, Position, Ref, Alt); }
}
