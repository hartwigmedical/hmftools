package com.hartwig.hmftools.gripss.common;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.SHORT_CALLING_SIZE;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_EVENT;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_HOMSEQ;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.gripss.filters.FilterType;

import htsjdk.variant.variantcontext.VariantContext;

public class SvData
{
    // VCF data
    private final String mId;
    private final StructuralVariantType mType;

    private final int mReferenceOrdinal;
    private final int mTumorOrdinal;
    private final Breakend[] mBreakends;

    // repeatedly used values for filtering are cached
    private final String mInsertSequence;
    private final boolean mImprecise;
    private final boolean mIsShortLocal;

    private final List<FilterType> mFilters;

    private int mPonCount;

    public SvData(final StructuralVariant sv, final GenotypeIds genotypeIds)
    {
        mId = sv.startContext().getAttributeAsString(VT_EVENT, sv.id());

        mType = sv.type();
        mReferenceOrdinal = genotypeIds.ReferenceOrdinal;
        mTumorOrdinal = genotypeIds.TumorOrdinal;

        Breakend breakendStart = Breakend.from(
                this, mType, true, sv.start(), sv.startContext(), genotypeIds.ReferenceOrdinal, genotypeIds.TumorOrdinal);

        Breakend breakendEnd = sv.end() != null ?
                Breakend.from(this, mType, false, sv.end(), sv.endContext(), genotypeIds.ReferenceOrdinal, genotypeIds.TumorOrdinal) : null;

        mBreakends = new Breakend[] { breakendStart, breakendEnd };

        mIsShortLocal = (mType == DEL || mType == DUP || mType == INS) && length() < SHORT_CALLING_SIZE;

        mFilters = Lists.newArrayList();

        mPonCount = 0;

        mImprecise = sv.imprecise();
        mInsertSequence = sv.insertSequence();
    }

    public String id() { return mId; }

    public String chromosomeStart() { return mBreakends[SE_START].Chromosome; }
    public String chromosomeEnd() { return !isSgl() ? mBreakends[SE_END].Chromosome : ""; }

    public int posStart() { return mBreakends[SE_START].Position; }
    public int posEnd() { return !isSgl() ? mBreakends[SE_END].Position : -1; }

    public byte orientStart() { return mBreakends[SE_START].Orientation; }
    public byte orientEnd() { return !isSgl() ? mBreakends[SE_END].Orientation : 0; }

    public StructuralVariantType type() { return mType; }
    public boolean isSgl() { return mType == SGL; }

    public Breakend[] breakends() { return mBreakends; }
    public Breakend breakendStart() { return mBreakends[SE_START]; }
    public Breakend breakendEnd() { return mBreakends[SE_END]; }

    public VariantContext contextStart() { return mBreakends[SE_START].Context; }
    public VariantContext contextEnd() { return !isSgl() ? mBreakends[SE_END].Context : null; }

    public String insertSequence() { return mInsertSequence; }
    public int insertSequenceLength() { return mInsertSequence.length(); }
    public int duplicationLength() { return mType == DUP ? length() + 1 : 0; }

    public boolean hasReference() { return mReferenceOrdinal >= 0; }

    public boolean isShortLocal() { return mIsShortLocal; }
    public boolean imprecise() { return mImprecise; }

    public List<FilterType> getFilters() { return mFilters; }

    public void addFilter(final FilterType filter)
    {
        if(!mFilters.contains(filter))
            mFilters.add(filter);
    }

    public void setPonCount(int count) { mPonCount = count; }
    public int getPonCount() { return mPonCount; }

    public static boolean hasLength(final StructuralVariantType type)
    {
        return type == INS || type == INV || type == DEL || type == DUP;
    }

    public int length() { return hasLength(mType) ? posEnd() - posStart() : 0; }

    public static int length(final StructuralVariant sv)
    {
        if(hasLength(sv.type()))
            return (int)(sv.position(false) - sv.position(true));

        return 0;
    }

    public String startHomology() { return contextStart().getAttributeAsString(VT_HOMSEQ, ""); }
    public String endHomology() { return contextEnd() != null ? contextEnd().getAttributeAsString(VT_HOMSEQ, "") : ""; }

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
