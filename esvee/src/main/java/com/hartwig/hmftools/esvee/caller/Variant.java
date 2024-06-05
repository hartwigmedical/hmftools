package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.sv.SvVcfTags.AVG_FRAG_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.SHORT_CALLING_SIZE;

import java.util.Arrays;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.esvee.caller.annotation.RepeatMaskAnnotation;
import com.hartwig.hmftools.esvee.common.FilterType;

import htsjdk.variant.variantcontext.VariantContext;

public class Variant
{
    // VCF data
    private final String mId;
    private final StructuralVariantType mType;

    private final int mReferenceOrdinal;
    private final Breakend[] mBreakends;

    // repeatedly used values for filtering are cached
    private final double mQual;
    private final String mInsertSequence;
    private boolean mIsShortLocal;
    private int mPonCount;
    private boolean mIsHotspot;
    private boolean mGermline;
    private RepeatMaskAnnotation mRmAnnotation;
    private final Set<FilterType> mFilters;

    public Variant(final StructuralVariant sv, final GenotypeIds genotypeIds)
    {
        // is this used for anything
        mId = "NONE";// sv.startContext().getAttributeAsString(EVENT, sv.id());

        mReferenceOrdinal = genotypeIds.ReferenceOrdinal;

        mIsShortLocal = (sv.type() == DEL || sv.type() == DUP || sv.type() == INS)
                && (sv.end().position() - sv.start().position()) < SHORT_CALLING_SIZE;

        // special handling of DUPs at the same base - need to switch their breakends
        if(sv.type() == INS && sv.position(true).equals(sv.position(false)))
        {
            mType = sv.type();

            mBreakends = new Breakend[] {
                    Breakend.from(this, true, sv.end(), sv.endContext(), genotypeIds.ReferenceOrdinal, genotypeIds.TumorOrdinal),
                    Breakend.from(this, false, sv.start(), sv.startContext(), genotypeIds.ReferenceOrdinal, genotypeIds.TumorOrdinal) };
        }
        else
        {
            mType = sv.type();

            Breakend breakendStart = Breakend.from(
                    this, true, sv.start(), sv.startContext(), genotypeIds.ReferenceOrdinal, genotypeIds.TumorOrdinal);

            Breakend breakendEnd = sv.end() != null ?
                    Breakend.from(this, false, sv.end(), sv.endContext(), genotypeIds.ReferenceOrdinal, genotypeIds.TumorOrdinal) : null;

            mBreakends = new Breakend[] { breakendStart, breakendEnd };
        }

        mQual = sv.qualityScore();

        mInsertSequence = sv.insertSequence();

        mFilters = Sets.newHashSet();

        // keep any non-pass filters from assembly
        addExistingFilters(sv.startContext());
        addExistingFilters(sv.endContext());

        mPonCount = 0;
        mIsHotspot = false;
        mGermline = false;
        mRmAnnotation = null;
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
    public String insertSequence() { return mInsertSequence; }

    public int ponCount() { return mPonCount; }
    public void setPonCount(int count) { mPonCount = count; }

    public RepeatMaskAnnotation getRmAnnotation() { return mRmAnnotation; }
    public void setRepeatMaskAnnotation(final RepeatMaskAnnotation annotation) { mRmAnnotation = annotation; }

    public boolean hasReference() { return mReferenceOrdinal >= 0; }

    public boolean isShortLocal() { return mIsShortLocal; }

    public static boolean hasLength(final StructuralVariantType type)
    {
        return type == INS || type == INV || type == DEL || type == DUP;
    }

    public int length() { return hasLength(mType) ? posEnd() - posStart() : 0; }

    public static int length(final StructuralVariant sv)
    {
        if(hasLength(sv.type()))
            return (sv.position(false) - sv.position(true));

        return 0;
    }

    public int averageFragmentLength() { return contextStart().getAttributeAsInt(AVG_FRAG_LENGTH, 0); }

    public int splitFragmentCount() { return contextStart().getAttributeAsInt(TOTAL_FRAGS, 0); }

    public void markHotspot() { mIsHotspot = true; }
    public boolean isHotspot() { return mIsHotspot; }

    public void markGermline() { mGermline = true; }
    public boolean isGermline() { return mGermline; }

    public void addFilter(final FilterType filter) { mFilters.add(filter); }
    public Set<FilterType> filters() { return mFilters; }
    public boolean isPass() { return mFilters.isEmpty(); }

    private void addExistingFilters(final VariantContext variantContext)
    {
        if(variantContext == null)
            return;

        // keep any non-pass filters from assembly
        for(String filterStr : variantContext.getFilters())
        {
            if(filterStr.equals(PASS))
                continue;

            FilterType filterType = Arrays.stream(FilterType.values()).filter(x -> x.vcfTag().equals(x)).findFirst().orElse(null);

            if(filterType != null)
                mFilters.add(filterType);
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
