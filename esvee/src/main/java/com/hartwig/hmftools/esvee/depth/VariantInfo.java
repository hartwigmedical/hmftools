package com.hartwig.hmftools.esvee.depth;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.isSingleBreakend;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.parseSingleOrientation;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.parseSvOrientation;
import static com.hartwig.hmftools.common.sv.SvUtils.isIndel;
import static com.hartwig.hmftools.common.sv.SvUtils.isShortLocalDelDupIns;
import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.fromRefAlt;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.sv.VariantAltInsertCoords;

import htsjdk.variant.variantcontext.VariantContext;

public class VariantInfo
{
    public final int Position;
    public final int PositionMin;
    public final int PositionMax;
    public final boolean IsSgl;
    public final boolean IsShortIndel;
    public final Orientation Orient;

    public final RefSupportCounts[] SampleSupportCounts;

    public VariantInfo(final VariantContext variant, final List<Integer> genotypeIds, double vafCap)
    {
        Position = variant.getStart();
        IsSgl = isSingleBreakend(variant);
        Orient = Orientation.fromByte(IsSgl ? parseSingleOrientation(variant) : parseSvOrientation(variant));

        if(!IsSgl)
        {
            StructuralVariantType svType = StructuralVariantType.fromContext(variant);

            if(isIndel(svType))
            {
                String ref = variant.getAlleles().get(0).getDisplayString();
                VariantAltInsertCoords altInsertCoords = fromRefAlt(variant.getAlleles().get(1).getDisplayString(), ref);
                int svLength = abs(Position - altInsertCoords.OtherPosition);
                IsShortIndel = isShortLocalDelDupIns(svType, svLength);
            }
            else
            {
                IsShortIndel = false;
            }
        }
        else
        {
            IsShortIndel = false;
        }

        final int[] homology = { 0, 0 };

        if(variant.hasAttribute(CIPOS))
        {
            final List<Integer> ihompos = variant.getAttributeAsIntList(CIPOS, 0);
            homology[0] = ihompos.get(0);
            homology[1] = ihompos.get(1);
        }

        PositionMin = Position + homology[0];
        PositionMax = Position + homology[1];

        SampleSupportCounts = new RefSupportCounts[genotypeIds.size()];

        for(int i = 0; i < genotypeIds.size(); ++i)
        {
            Integer genotypeIndex = genotypeIds.get(i);
            int sampleFragments = getGenotypeAttributeAsInt(variant.getGenotype(genotypeIndex), TOTAL_FRAGS, 0);
            int refFragsCap = vafCap > 0 ? (int) (max(sampleFragments, 1) / vafCap) : 0;
            SampleSupportCounts[i] = new RefSupportCounts(refFragsCap);
        }
    }

    public RefSupportCounts totalSupport()
    {
        if(SampleSupportCounts.length == 1)
            return SampleSupportCounts[0];

        RefSupportCounts totalCounts = new RefSupportCounts(0);

        for(int i = 0; i < SampleSupportCounts.length; ++i)
        {
            totalCounts.RefSupport += SampleSupportCounts[i].RefSupport;
            totalCounts.RefPairSupport += SampleSupportCounts[i].RefPairSupport;
        }

        return totalCounts;
    }

    public String toString()
    {
        return format("pos(%d %d-%d)", Position, PositionMin, PositionMax);
    }
}
