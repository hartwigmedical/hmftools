package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.seIndex;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.svtools.germline.FilterConstants.SHORT_CALLING_SIZE;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.BQ;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.BVF;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.HOMSEQ;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.QUAL;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.REF;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.REFPAIR;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.VF;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class SvData
{
    // VCF data
    private final String mVcfId;
    private final StructuralVariantType mType;

    private final int mReferenceOrdinal;
    private final int mTumorOrdinal;
    private final Breakend[] mBreakends;

    // repeatedly used values for filtering are cached
    private final String mAlt;
    private final String mInsertSequence;
    private final boolean mImprecise;
    private final boolean mIsShortLocal;

    private final List<FilterType> mFilters;

    private int mPonCount;
    private String mHotspotGeneInfo;

    public SvData(final StructuralVariant sv, final GenotypeIds genotypeIds)
    {
        mVcfId = sv.id();
        mType = sv.type();
        mReferenceOrdinal = genotypeIds.ReferenceOrdinal;
        mTumorOrdinal = genotypeIds.TumorOrdinal;

        Breakend breakendStart = Breakend.from(
                mVcfId, mType, sv.start(), sv.startContext(), genotypeIds.ReferenceOrdinal, genotypeIds.TumorOrdinal);

        Breakend breakendEnd = sv.end() != null ?
                Breakend.from(mVcfId, mType, sv.end(), sv.endContext(), genotypeIds.ReferenceOrdinal, genotypeIds.TumorOrdinal) : null;

        mBreakends = new Breakend[] { breakendStart, breakendEnd };

        mIsShortLocal = (mType == DEL || mType == DUP || mType == INS) && length() < SHORT_CALLING_SIZE;

        mFilters = Lists.newArrayList();

        mPonCount = 0;
        mHotspotGeneInfo = "";

        mImprecise = sv.imprecise();
        mInsertSequence = sv.insertSequence();

        // TODO
        mAlt = "";
    }

    public String id() { return mVcfId; }

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

    /*
    public Genotype refGenotype() { return mContexts[SE_START].getGenotype(mReferenceOrdinal); }
    public Genotype tumorGenotype() { return mContexts[SE_START].getGenotype(mTumorOrdinal); }

    public Allele refAllele() { return mReferenceOrdinal >= 0 ? mContexts[SE_START].getAlleles().get(mReferenceOrdinal) : null; }
    public Allele tumorAllele() { return mContexts[SE_START].getAlleles().get(mTumorOrdinal); }
    */

    public String altString()
    {
        // TODO
        if(isSgl())
            return ""; // mOrientation[SE_START] == POS_ORIENT ? mAlt + mInsertSequence : mInsertSequence + mAlt;
        else
            return mAlt;
    }

    public String insertSequence() { return mInsertSequence; }
    public int insertSequenceLength() { return mInsertSequence.length(); }

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

    public void setHotspotGeneInfo(final String info) { mHotspotGeneInfo = info; }
    public String getHotspotGeneInfo() { return mHotspotGeneInfo; }

    public boolean hasLength() { return hasLength(mType); }

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

    public String startHomology() { return contextStart().getAttributeAsString(HOMSEQ, ""); }
    public String endHomology() { return contextEnd() != null ? contextEnd().getAttributeAsString(HOMSEQ, "") : ""; }


    {
        /*
            Paired logic:

            val match = Regex(BREAKPOINT_REGEX).find(alt)!!
            val (initialSequence, bracket, location, finalSequence) = match.destructured
            val (otherChromosome, otherPosition) = LocationString(location)

            val endOrientation: Byte = if (bracket == "]") 1 else -1

            val startOrientation: Byte
            val insertSequence: String
            val altBase: String
            if (initialSequence.isNotEmpty()) {
                startOrientation = 1
                insertSequence = initialSequence.substring(ref.length, initialSequence.length)
                altBase = alt.substring(0, 1)
            } else {
                startOrientation = -1
                insertSequence = finalSequence.substring(0, finalSequence.length - ref.length)
                altBase = alt.substring(alt.length - 1, alt.length)
            }

            if (alt.startsWith(".")) {
                return Single(alt.substring(alt.length - 1, alt.length), alt.substring(ref.length, alt.length - 1), -1)
            }

            if (alt.endsWith(".")) {
                return Single(alt.substring(0, 1), alt.substring(1, alt.length - ref.length), 1)
            }

            return paired(chromosome, position, ref, alt)

         */
    }

    public String toString()
    {
        if(!isSgl())
        {
            return String.format("%s:%s pos(%s:%d - %s:%d)",
                    mVcfId, mType.toString(), chromosomeStart(), posStart(), chromosomeEnd(), posEnd());
        }
        else
        {
            return String.format("%s:%s pos(%s:%d)",
                    mVcfId, mType.toString(), chromosomeStart(), posStart());
        }
    }

}
