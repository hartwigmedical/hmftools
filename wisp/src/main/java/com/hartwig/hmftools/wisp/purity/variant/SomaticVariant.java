package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

public class SomaticVariant
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;
    public final VariantType Type;
    public final double SubclonalPerc;
    public final double CopyNumber;
    public final double VariantCopyNumber;
    public final double Mappability;
    public final VariantTier Tier;
    public final boolean Reported;
    public final int RepeatCount;
    public final boolean Hotspot;
    public final String TriNucContext;
    public final String CanonicalGeneName;
    public final String CanonicalCodingEffect;

    public final List<GenotypeFragments> Samples;

    private List<FilterReason> mFilterReasons;

    private double mSequenceGcRatio;
    private boolean mIsProbeVariant;

    public SomaticVariant(
            final VariantContextDecorator variantDecorator, final double subclonalPerc, final List<FilterReason> filterReasons,
            boolean hasSyntheticTumor)
    {
        Chromosome = variantDecorator.chromosome();
        Position = variantDecorator.position();
        Ref = variantDecorator.ref();
        Alt = variantDecorator.alt();
        Type = variantDecorator.type();
        SubclonalPerc = subclonalPerc;
        CopyNumber = !hasSyntheticTumor ? variantDecorator.adjustedCopyNumber() : 2;
        VariantCopyNumber = !hasSyntheticTumor ? variantDecorator.variantCopyNumber() : 1;
        Tier = variantDecorator.tier();
        Reported = variantDecorator.reported();
        RepeatCount = variantDecorator.repeatCount();
        Hotspot = variantDecorator.isHotspot();
        TriNucContext = variantDecorator.trinucleotideContext();
        Mappability = variantDecorator.mappability();
        VariantImpact variantImpact = variantDecorator.variantImpact();
        CanonicalGeneName = variantImpact.GeneName;
        CanonicalCodingEffect = variantImpact.CanonicalCodingEffect.toString();

        Samples = Lists.newArrayList();
        mFilterReasons = filterReasons;
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

    public String toString() { return format("%s:%d %s>%s", Chromosome, Position, Ref, Alt); }
}
