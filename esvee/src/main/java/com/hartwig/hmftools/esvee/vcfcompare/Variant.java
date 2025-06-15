package com.hartwig.hmftools.esvee.vcfcompare;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SV_ID;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.variant.GenotypeIds;

import htsjdk.variant.variantcontext.VariantContext;

public class Variant
{
    private final String mId;
    private final StructuralVariantType mType;

    private final Breakend[] mBreakends;

    // repeatedly used values for filtering are cached
    private final double mQual;

    public Variant(final StructuralVariant sv, final GenotypeIds genotypeIds)
    {
        // is this used for anything
        mId = sv.startContext().getAttributeAsString(SV_ID, sv.id());

        // special handling of DUPs at the same base - need to switch their breakends
        mType = sv.type();

        Breakend breakendStart = Breakend.from(
                this, true, sv.start(), sv.startContext(), genotypeIds.ReferenceOrdinal, genotypeIds.TumorOrdinal);

        Breakend breakendEnd = sv.end() != null ?
                Breakend.from(this, false, sv.end(), sv.endContext(), genotypeIds.ReferenceOrdinal, genotypeIds.TumorOrdinal) : null;

        mBreakends = new Breakend[] { breakendStart, breakendEnd };

        mQual = sv.qualityScore();
    }

    public String id() { return mId; }

    public String chromosomeStart() { return mBreakends[SE_START].Chromosome; }
    public String chromosomeEnd() { return !isSgl() ? mBreakends[SE_END].Chromosome : ""; }

    public int posStart() { return mBreakends[SE_START].Position; }
    public int posEnd() { return !isSgl() ? mBreakends[SE_END].Position : -1; }

    public Orientation orientStart() { return mBreakends[SE_START].Orient; }
    public Orientation orientEnd() { return !isSgl() ? mBreakends[SE_END].Orient : null; }

    public StructuralVariantType type() { return mType; }
    public boolean isSgl() { return mType == SGL; }

    public Breakend[] breakends() { return mBreakends; }
    public Breakend breakendStart() { return mBreakends[SE_START]; }
    public Breakend breakendEnd() { return mBreakends[SE_END]; }

    public VariantContext contextStart() { return mBreakends[SE_START].Context; }
    public VariantContext contextEnd() { return !isSgl() ? mBreakends[SE_END].Context : null; }

    public double qual() { return mQual; }

    public boolean isLineSite()
    {
        for(Breakend breakend : mBreakends)
        {
            if(breakend != null && (breakend.isLine()))
                return true;
        }

        return false;
    }

    public boolean inChainedAssembly()
    {
        for(Breakend breakend : mBreakends)
        {
            if(breakend != null && breakend.inChainedAssembly())
                return true;
        }

        return false;
    }

    public boolean isPass() { return mBreakends[0].isPass(); }
    public boolean isFiltered() { return !isPass(); }

    public String svString()
    {
        if(!isSgl())
        {
            return String.format("%s %s:%d:%s_%s:%d:%s",
                    mType.toString(), chromosomeStart(), posStart(), orientStart(),
                    chromosomeEnd(), posEnd(), orientEnd());
        }
        else
        {
            return String.format("%s %s:%d:%s",
                    mType.toString(), chromosomeStart(), posStart(), orientStart().toString());
        }
    }


    public String toString()
    {
        if(!isSgl())
        {
            return String.format("%s:%s pos(%s:%d - %s:%d)",
                    mId, mType.toString(), chromosomeStart(), posStart(), chromosomeEnd(), posEnd());
        }
        else
        {
            return String.format("%s:%s pos(%s:%d)",
                    mId, mType.toString(), chromosomeStart(), posStart());
        }
    }
}
