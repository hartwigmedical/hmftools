package com.hartwig.hmftools.svprep.append;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.UMI_TYPE_ATTRIBUTE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.parseSingleOrientation;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.parseSvOrientation;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.samtools.UmiReadType;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.svprep.reads.JunctionData;
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
    private final int[] mUmiTypeCounts;
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
        mUmiTypeCounts = new int[UmiReadType.values().length];
        mDepth = 0;
    }

    public static BreakendData fromVariant(final VariantContext variant)
    {
        boolean isSingle = StructuralVariantFactory.isSingleBreakend(variant);
        byte orientation = isSingle ? parseSingleOrientation(variant) : parseSvOrientation(variant);
        return new BreakendData(variant, variant.getContig(), variant.getStart(), orientation, isSingle);
    }

    public VariantContext variant() { return mVariant; }

    public void addJunctionData(final List<JunctionData> junctions)
    {
        Set<String> processedReads = Sets.newHashSet();

        for(JunctionData junctionData : junctions)
        {
            for(Map.Entry<ReadType, List<ReadRecord>> entry : junctionData.ReadTypeReads.entrySet())
            {
                ReadType readType = entry.getKey();
                List<ReadRecord> reads = entry.getValue();

                if(!supportsJunction(readType))
                    continue;

                for(ReadRecord read : reads)
                {
                    if(processedReads.contains(read.id()))
                        continue;

                    processedReads.add(read.id());

                    ++mReadTypeSupport[readType.ordinal()];

                    /* rework or move into Esvee
                    String umiType = read.record().getStringAttribute(UMI_TYPE_ATTRIBUTE);

                    UmiReadType umiReadType = umiType != null ? UmiReadType.valueOf(umiType) : UmiReadType.NONE;
                    ++mUmiTypeCounts[umiReadType.ordinal()];
                    */
                }
            }

            mDepth = max(mDepth, junctionData.depth());
        }
    }

    private static boolean supportsJunction(final ReadType readType)
    {
        return readType == ReadType.JUNCTION || readType == ReadType.EXACT_SUPPORT; //  || readType == ReadType.SUPPORT
    }

    public final int[] readTypeSupport() { return mReadTypeSupport; }

    public int totalSupport()
    {
        int total = 0;

        for(ReadType readType : ReadType.values())
        {
            if(supportsJunction(readType))
                total += mReadTypeSupport[readType.ordinal()];
        }

        return total;
    }

    public int depth() { return mDepth; }

    public int[] umiTypeCounts() { return mUmiTypeCounts; }

    public String toString() { return format("%s:%d:%d", mVariant.getContig(), Position, Orientation); }

}
