package com.hartwig.hmftools.svprep.append;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.parseSingleOrientation;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.parseSvOrientation;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.svprep.reads.JunctionData;
import com.hartwig.hmftools.svprep.reads.ReadGroup;
import com.hartwig.hmftools.svprep.reads.ReadRecord;
import com.hartwig.hmftools.svprep.reads.ReadType;

import htsjdk.variant.variantcontext.VariantContext;

public class BreakendData
{
    public final String Chromosome;
    public final int Position;
    public final byte Orientation;
    public final boolean IsSingle;

    private final int[] mReadTypeSupport;
    private int mDepth;

    private final VariantContext mVariant;

    public BreakendData(
            final VariantContext variant, final String chromosome, final int position, final byte orientation, final boolean isSingle)
    {
        Chromosome = chromosome;
        Position = position;
        Orientation = orientation;
        IsSingle = isSingle;
        mVariant = variant;

        mReadTypeSupport = new int[ReadType.values().length];
        mDepth = 0;
    }

    public static BreakendData fromVariant(final VariantContext variant)
    {
        boolean isSingle = StructuralVariantFactory.isSingleBreakend(variant);
        byte orientation = isSingle ? parseSingleOrientation(variant) : parseSvOrientation(variant);
        return new BreakendData(variant, variant.getContig(), variant.getStart(), orientation, isSingle);
    }

    public VariantContext variant() { return mVariant; }

    public void setJunctionData(final JunctionData junctionData)
    {
        for(Map.Entry<ReadType, List<ReadRecord>> entry : junctionData.ReadTypeReads.entrySet())
        {
            mReadTypeSupport[entry.getKey().ordinal()] = entry.getValue().size();
        }

        mDepth = junctionData.depth();
    }

    public final int[] readTypeSupport() { return mReadTypeSupport; }

    public int totalSupport()
    {
        return mReadTypeSupport[ReadType.JUNCTION.ordinal()]
                + mReadTypeSupport[ReadType.EXACT_SUPPORT.ordinal()]
                + mReadTypeSupport[ReadType.SUPPORT.ordinal()];
    }

    public int depth() { return mDepth; }

}
