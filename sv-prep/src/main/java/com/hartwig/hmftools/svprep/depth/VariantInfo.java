package com.hartwig.hmftools.svprep.depth;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.CIPOS;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.VARIANT_FRAGMENT_BREAKEND_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.VARIANT_FRAGMENT_BREAKPOINT_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.isSingleBreakend;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.util.List;

import htsjdk.variant.variantcontext.VariantContext;

class VariantInfo
{
    public final int Position;
    public final int PositionMin;
    public final int PositionMax;
    public final byte Orientation;
    public final int FragmentCount;
    public final int RefFragsCap;

    public final RefSupportCounts[] SampleSupportCounts;

    public VariantInfo(final VariantContext variant, int sampleCount, double vafCap)
    {
        Position = variant.getStart();
        Orientation = getOrientation(variant);

        final int[] homology = { 0, 0 };

        if(variant.hasAttribute(CIPOS))
        {
            final List<Integer> ihompos = variant.getAttributeAsIntList(CIPOS, 0);
            homology[0] = ihompos.get(0);
            homology[1] = ihompos.get(1);
        }

        PositionMin = Position + homology[0];
        PositionMax = Position + homology[1];

        FragmentCount = max(
                variant.getAttributeAsInt(VARIANT_FRAGMENT_BREAKPOINT_COVERAGE, 0),
                variant.getAttributeAsInt(VARIANT_FRAGMENT_BREAKEND_COVERAGE, 0));

        RefFragsCap = vafCap > 0 ? (int) (FragmentCount / vafCap) : 0;

        SampleSupportCounts = new RefSupportCounts[sampleCount];

        for(int i = 0; i < sampleCount; ++i)
        {
            SampleSupportCounts[i] = new RefSupportCounts();
        }
    }

    public RefSupportCounts totalSupport()
    {
        if(SampleSupportCounts.length == 1)
            return SampleSupportCounts[0];

        RefSupportCounts totalCounts = new RefSupportCounts();
        for(int i = 0; i < SampleSupportCounts.length; ++i)
        {
            totalCounts.RefSupport += SampleSupportCounts[i].RefSupport;
            totalCounts.RefPairSupport += SampleSupportCounts[i].RefPairSupport;
        }

        return totalCounts;
    }

    private static byte getOrientation(final VariantContext variant)
    {
        String alt = variant.getAlternateAllele(0).getDisplayString();

        if(isSingleBreakend(variant))
        {
            return alt.startsWith(".") ? NEG_ORIENT : POS_ORIENT;
        }
        else
        {
            return alt.startsWith("]") || alt.startsWith("[") ? NEG_ORIENT : POS_ORIENT;
        }
    }

    public String toString()
    {
        return format("pos(%d %d-%d) Frags(%d cap=%d)", Position, PositionMin, PositionMax, FragmentCount, RefFragsCap);
    }
}
