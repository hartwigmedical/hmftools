package com.hartwig.hmftools.common.variant.structural.annotation;

import static com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation.isDownstream;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class Transcript {

    @NotNull
    private final GeneAnnotation mGene;

    private final int mTransId;
    private final String mStableId;
    private final int mExonUpstream;
    private final int mExonUpstreamPhase;
    private final int mExonDownstream;
    private final int mExonDownstreamPhase;
    private final int mExonMax;
    private int mExactCodingBase;

    private final int mCodingBases;
    private final int mTotalCodingBases;

    private final boolean mCanonical;
    private String mBioType;

    private String mCodingType;
    private final String mRegionType;

    @Nullable
    private final Long mCodingStart;

    @Nullable
    private final Long mCodingEnd;

    private final long mTranscriptStart;
    private final long mTranscriptEnd;

    private int mExonDistanceUp;
    private int mExonDistanceDown;

    private boolean mIsDisruptive;

    private String mProteinFeaturesKept;
    private String mProteinFeaturesLost;

    public static String TRANS_REGION_TYPE_UPSTREAM = "Upstream";
    public static String TRANS_REGION_TYPE_EXONIC = "Exonic";
    public static String TRANS_REGION_TYPE_INTRONIC = "Intronic";

    public static String TRANS_CODING_TYPE_CODING = "Coding";
    public static String TRANS_CODING_TYPE_UPSTREAM = "5P_UTR";
    public static String TRANS_CODING_TYPE_DOWNSTREAM = "3P_UTR";
    public static String TRANS_CODING_TYPE_NON_CODING = "NonCoding";

    private static int STOP_CODON_LENGTH = 3;

    public Transcript(@NotNull final GeneAnnotation parent, int transId, final String stableId,
            final int exonUpstream, final int exonUpstreamPhase, final int exonDownstream, final int exonDownstreamPhase,
            final int codingBases, final int totalCodingBases,
            final int exonMax, final boolean canonical, final long transcriptStart, final long transcriptEnd,
            final Long codingStart, final Long codingEnd)
    {
        mGene = parent;
        mTransId = transId;
        mStableId = stableId;

        mExonUpstream = exonUpstream;
        mExonDownstream = exonDownstream;

        mExactCodingBase = -1;

        mExonMax = exonMax;
        mBioType = "";
        mExonDistanceUp = 0;
        mExonDistanceDown = 0;

        if(totalCodingBases > STOP_CODON_LENGTH)
        {
            // remove the stop codon from what is consider coding
            mTotalCodingBases = totalCodingBases - STOP_CODON_LENGTH;
        }
        else
        {
            mTotalCodingBases = 0;
        }

        if(codingBases > mTotalCodingBases)
        {
            // limit as well
            mCodingBases = mTotalCodingBases;
        }
        else
        {
            mCodingBases = codingBases;
        }

        mCanonical = canonical;
        mCodingStart = codingStart;
        mCodingEnd = codingEnd;
        mTranscriptStart = transcriptStart;
        mTranscriptEnd = transcriptEnd;

        mCodingType = calcCodingType();
        mRegionType = calcRegionType();

        mExonDownstreamPhase = exonDownstreamPhase;
        mExonUpstreamPhase = exonUpstreamPhase;

        if(isDownstream(mGene) && mRegionType == TRANS_REGION_TYPE_UPSTREAM)
            mIsDisruptive = false;
        else
            mIsDisruptive = true;

        mProteinFeaturesKept = "";
        mProteinFeaturesLost = "";
    }

    public int transId() { return mTransId; }
    public final String transcriptId() { return mStableId; }

    public boolean isExonic()
    {
        return mExonUpstream > 0 && mExonUpstream == mExonDownstream;
    }

    public boolean isPromoter()
    {
        return mExonUpstream == 0 && (mExonDownstream == 1 || mExonDownstream == 2);
    }

    public boolean isIntronic()
    {
        return mExonUpstream > 0 && (mExonDownstream - mExonUpstream) == 1;
    }

    public long getDistanceUpstream()
    {
        if(!isPromoter())
            return 0;

        if(mGene.strand() == 1)
            return mTranscriptStart - svPosition();
        else
            return svPosition() - mTranscriptEnd;
    }

    public final String codingType() { return mCodingType; }
    public final String regionType() { return mRegionType; }

    // for convenience
    public boolean isCoding() { return mCodingType.equals(TRANS_CODING_TYPE_CODING); }

    public boolean preCoding()
    {
        return mCodingType.equals(TRANS_CODING_TYPE_UPSTREAM);
    }

    public boolean postCoding()
    {
        return mCodingType.equals(TRANS_CODING_TYPE_DOWNSTREAM);
    }

    public boolean nonCoding()
    {
        return mCodingType.equals(TRANS_CODING_TYPE_NON_CODING);
    }

    private final String calcRegionType()
    {
        if(isIntronic())
            return TRANS_REGION_TYPE_INTRONIC;

        if(isExonic())
            return TRANS_REGION_TYPE_EXONIC;

        if(isPromoter())
            return TRANS_REGION_TYPE_UPSTREAM;

        return "Unknown";
    }

    private final String calcCodingType()
    {
        if(mCodingStart == null || mCodingEnd == null || mTotalCodingBases == 0)
        {
            return TRANS_CODING_TYPE_NON_CODING;
        }
        else if(mCodingBases == 0)
        {
            return TRANS_CODING_TYPE_UPSTREAM;
        }
        else if(mCodingBases == mTotalCodingBases)
        {
            return TRANS_CODING_TYPE_DOWNSTREAM;
        }
        else
        {
            return TRANS_CODING_TYPE_CODING;
        }
    }

    public void setCodingType(final String type)
    {
        mCodingType = type;
    }

    public boolean isCanonical() { return mCanonical; }

    public void setExonDistances(int up, int down)
    {
        mExonDistanceUp = up;
        mExonDistanceDown = down;
    }

    public int exonDistanceUp() { return mExonDistanceUp; }
    public int exonDistanceDown() { return mExonDistanceDown; }

    public void setBioType(final String type) { mBioType = type; }
    public final String bioType() { return mBioType; }

    @NotNull
    public GeneAnnotation parent() { return mGene; }

    public long svPosition() { return mGene.position(); }

    public String geneName() { return mGene.geneName(); }

    public int exonUpstream() { return mExonUpstream; }
    public int exonUpstreamPhase() { return mExonUpstreamPhase; }

    public int codingBases() { return mCodingBases; }
    public int calcCodingBases(boolean isUpstream) { return isUpstream ? mCodingBases : mTotalCodingBases - mCodingBases; }
    public int totalCodingBases() { return mTotalCodingBases; }

    public int exonDownstream() { return mExonDownstream; }
    public int exonDownstreamPhase() { return mExonDownstreamPhase; }

    public int exonMax() { return mExonMax; }

    public int exactCodingBase() { return mExactCodingBase; }
    public void setExactCodingBase(int base) { mExactCodingBase = base; }

    public long transcriptStart() { return mTranscriptStart; }
    public long transcriptEnd() { return mTranscriptEnd; }

    public final long length() { return mTranscriptEnd - mTranscriptStart; }

    public boolean isDisruptive() { return mIsDisruptive; }
    public void setIsDisruptive(boolean toggle) { mIsDisruptive = toggle; }

    public long codingStart() { return mCodingStart != null ? mCodingStart : 0; }
    public long codingEnd() { return mCodingEnd != null ? mCodingEnd : 0; }

    public final String toString()
    {
        return mGene.geneName() + " " + mStableId;
    }

    public void addProteinFeature(final String feature, boolean isPreserved)
    {
        if(isPreserved)
        {
            mProteinFeaturesKept += mProteinFeaturesKept.isEmpty() ? feature : ";" + feature;
        }
        else
        {
            mProteinFeaturesLost += mProteinFeaturesLost.isEmpty() ? feature : ";" + feature;
        }
    }

    public final String getProteinFeaturesKept() { return mProteinFeaturesKept; }
    public final String getProteinFeaturesLost() { return mProteinFeaturesLost; }

}
