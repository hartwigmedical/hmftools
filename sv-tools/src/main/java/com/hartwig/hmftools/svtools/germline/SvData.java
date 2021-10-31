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
    private final String[] mChromosome;
    private final int[] mPosition;
    private final byte[] mOrientation;

    private final int mReferenceOrdinal;
    private final int mTumorOrdinal;
    private VariantContext[] mContexts;

    // repeatedly used values for filtering are cached
    private final int mTumorFragments;
    private final int mReferenceFragments;
    private final int mReferenceReads;
    private final int mReferencePairReads;
    private final double mTumorQuality;
    private final String mAlt;
    private final String mInsertSequence;
    private final boolean mIsShortLocal;

    private final List<FilterType> mFilters;

    private int mPonCount;
    private String mHotspotGeneInfo;

    // classification state
    private String[] mAsmbSvIds;

    public SvData(
            final String id, final StructuralVariantType type,
            final String chrStart, final String chrEnd, int posStart, int posEnd, byte orientStart, byte orientEnd,
            final VariantContext contextStart, final VariantContext contextEnd, final int referenceOrdinal, final int tumorOrdinal)
    {
        mVcfId = id;
        mType = type;
        mChromosome = new String[] { chrStart, chrEnd };
        mPosition = new int[] { posStart, posEnd };
        mOrientation = new byte[] { orientStart, orientEnd };
        mReferenceOrdinal = referenceOrdinal;
        mTumorOrdinal = tumorOrdinal;
        mContexts = new VariantContext[] { contextStart, contextEnd };

        mFilters = Lists.newArrayList();

        mPonCount = 0;
        mHotspotGeneInfo = "";
        mAsmbSvIds = new String[SE_PAIR];

        final Genotype refGenotype = mContexts[SE_START].getGenotype(mReferenceOrdinal);
        final Genotype tumorGenotype = mContexts[SE_START].getGenotype(mTumorOrdinal);

        if(isSgl())
        {
            mReferenceFragments = VcfUtils.getGenotypeAttributeAsInt(refGenotype, BVF, 0);
            mTumorFragments = VcfUtils.getGenotypeAttributeAsInt(tumorGenotype, BVF, 0);
            mTumorQuality = VcfUtils.getGenotypeAttributeAsDouble(tumorGenotype, BQ, 0);
        }
        else
        {
            mReferenceFragments = VcfUtils.getGenotypeAttributeAsInt(refGenotype, VF, 0);
            mTumorFragments = VcfUtils.getGenotypeAttributeAsInt(tumorGenotype, VF, 0);
            mTumorQuality = VcfUtils.getGenotypeAttributeAsDouble(tumorGenotype, QUAL, 0);
        }

        mIsShortLocal = (mType == DEL || mType == DUP || mType == INS) && length() < SHORT_CALLING_SIZE;

        mReferenceReads = VcfUtils.getGenotypeAttributeAsInt(refGenotype, REF, 0);
        mReferencePairReads = VcfUtils.getGenotypeAttributeAsInt(refGenotype, REFPAIR, 0);

        // set alt and insert sequence
        mAlt = "";
        mInsertSequence = "";
    }

    public static SvData from(final StructuralVariant sv, final GenotypeIds genotypeIds)
    {
        boolean sglBreakend = sv.type() == SGL;

        return new SvData(
                sv.id(), sv.type(), sv.chromosome(true), !sglBreakend ? sv.chromosome(false) : "0",
                sv.position(true).intValue(), !sglBreakend ? sv.position(false).intValue() : -1,
                sv.orientation(true), !sglBreakend ? sv.orientation(false) : 0,
                sv.startContext(), sv.endContext(), genotypeIds.ReferenceOrdinal, genotypeIds.TumorOrdinal);
    }

    public String id() { return mVcfId; }

    public String[] chromosomes() { return mChromosome; }
    public String chromosomeStart() { return mChromosome[SE_START]; }
    public String chromosomeEnd() { return mChromosome[SE_END]; }

    public int[] positions() { return mPosition; }
    public int posStart() { return mPosition[SE_START]; }
    public int posEnd() { return mPosition[SE_END]; }

    public byte[] orientations() { return mOrientation; }
    public byte orientStart() { return mOrientation[SE_START]; }
    public byte orientEnd() { return mOrientation[SE_END]; }

    public StructuralVariantType type() { return mType; }
    public boolean isSgl() { return mType == SGL; }

    public VariantContext[] contexts() { return mContexts; }
    public VariantContext contextStart() { return mContexts[SE_START]; }
    public VariantContext contextEnd() { return mContexts[SE_END]; }

    public Genotype refGenotype() { return mContexts[SE_START].getGenotype(mReferenceOrdinal); }
    public Genotype tumorGenotype() { return mContexts[SE_START].getGenotype(mTumorOrdinal); }

    public Allele refAllele() { return mReferenceOrdinal >= 0 ? mContexts[SE_START].getAlleles().get(mReferenceOrdinal) : null; }
    public Allele tumorAllele() { return mContexts[SE_START].getAlleles().get(mTumorOrdinal); }

    public String altString()
    {
        if(isSgl())
            return mOrientation[SE_START] == POS_ORIENT ? mAlt + mInsertSequence : mInsertSequence + mAlt;
        else
            return mAlt;
    }

    public String insertSequence() { return mInsertSequence; }

    public boolean hasReference() { return mReferenceOrdinal >= 0; }

    public int tumorFragments() { return mTumorFragments; }
    public int referenceFragments() { return mReferenceFragments; }
    public int referenceReads() { return mReferenceReads; }
    public int referencePairReads() { return mReferencePairReads; }
    public double tumorQuality() { return mTumorQuality; }
    public boolean isShortLocal() { return mIsShortLocal; }

    public void addFilter(final FilterType filter)
    {
        if(!mFilters.contains(filter))
            mFilters.add(filter);
    }

    public List<FilterType> getFilters() { return mFilters; }

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

    public String assemblySvIds(boolean isStart) { return mAsmbSvIds[seIndex(isStart)]; }
    public void setAssemblySvId(int seIndex, final String svId) { mAsmbSvIds[seIndex] = svId; }

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
        return String.format("%s:%s %s:%d - %s:%d",
                mVcfId, mType.toString(), mChromosome[SE_START], mPosition[SE_START], mChromosome[SE_END], mPosition[SE_END]);
    }

}
