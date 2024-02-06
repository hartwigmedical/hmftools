package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.String.format;

import static com.hartwig.hmftools.wisp.purity.variant.FilterReason.NO_FILTER;

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

    private List<FilterReason> mFilterReasons;

    private final VariantContextDecorator mVariant;
    private double mSequenceGcRatio;
    private boolean mIsProbeVariant;

    public SomaticVariant(final VariantContextDecorator variant, final double subclonalPerc, final List<FilterReason> filterReasons)
    {
        Chromosome = variant.chromosome();
        Position = variant.position();
        Ref = variant.ref();
        Alt = variant.alt();
        Type = variant.type();
        SubclonalPerc = subclonalPerc;
        Samples = Lists.newArrayList();
        mFilterReasons = filterReasons;
        mVariant = variant;
        mSequenceGcRatio = 0;
        mIsProbeVariant = false;
    }

    public void addFilterReason(final FilterReason filterReason) { mFilterReasons.add(filterReason); }

    public List<FilterReason> filterReasons() { return mFilterReasons; }

    public boolean isFiltered() { return !mFilterReasons.isEmpty(); }

    public boolean isProbeVariant() { return mIsProbeVariant; }
    public void markProbeVariant() { mIsProbeVariant = true; }

    public GenotypeFragments findGenotypeData(final String sampleId)
    {
        return Samples.stream().filter(x -> x.SampleName.equals(sampleId)).findFirst().orElse(null);
    }

    public void setSequenceGcRatio(double ratio) { mSequenceGcRatio = ratio; }
    public double sequenceGcRatio() { return mSequenceGcRatio; }

    public double copyNumber() { return mVariant.adjustedCopyNumber(); }
    public double variantCopyNumber() { return mVariant.variantCopyNumber(); }

    public VariantContextDecorator decorator() { return mVariant; }

    public String toString() { return format("%s:%d %s>%s", Chromosome, Position, Ref, Alt); }
}
