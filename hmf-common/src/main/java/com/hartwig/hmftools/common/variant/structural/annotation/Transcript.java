package com.hartwig.hmftools.common.variant.structural.annotation;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class Transcript {

    @NotNull
    private final GeneAnnotation mGene;

    private final String mTranscriptId;
    private final int mExonUpstream;
    private final int mExonUpstreamPhase;
    private final int mExonDownstream;
    private final int mExonDownstreamPhase;
    private final int mExonMax;

    private final long mCodingBases;
    private final long mTotalCodingBases;

    private final boolean mCanonical;

    private final String mCodingType;
    private final String mRegionType;

    @Nullable
    private final Long mCodingStart;

    @Nullable
    private final Long mCodingEnd;

    private final long mTranscriptStart;
    private final long mTranscriptEnd;

    public static String TRANS_REGION_TYPE_PROMOTOR = "Promotor";
    public static String TRANS_REGION_TYPE_EXONIC = "Exonic";
    public static String TRANS_REGION_TYPE_INTRONIC = "Intronic";

    public static String TRANS_CODING_TYPE_CODING = "Coding";
    public static String TRANS_CODING_TYPE_UPSTREAM = "Upstream";
    public static String TRANS_CODING_TYPE_DOWNSTREAM = "Downstream";
    public static String TRANS_CODING_TYPE_NON_CODING = "NonCoding";

    public static int PROMOTOR_REGION_MAX = 20000;

    private static int STOP_CODON_LENGTH = 3;

    public Transcript(@NotNull final GeneAnnotation parent, @NotNull final String transcriptId,
            final int exonUpstream, final int exonUpstreamPhase, final int exonDownstream, final int exonDownstreamPhase,
            final long codingBases, final long totalCodingBases,
            final int exonMax, final boolean canonical, final long transcriptStart, final long transcriptEnd,
            @Nullable final Long codingStart, @Nullable final Long codingEnd)
    {
        mGene = parent;
        mTranscriptId = transcriptId;

        mExonUpstream = exonUpstream;
        mExonUpstreamPhase = exonUpstreamPhase;
        mExonDownstream = exonDownstream;
        mExonDownstreamPhase = exonDownstreamPhase;

        mExonMax = exonMax;

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

    }

    @NotNull
    public String transcriptId() { return mTranscriptId; }

    public boolean isExonic()
    {
        return mExonUpstream > 0 && mExonUpstream == mExonDownstream;
    }

    public boolean isPromoter()
    {
        return mExonUpstream == 0 && mExonDownstream == 1;
    }

    public boolean isIntronic()
    {
        return mExonUpstream > 0 && (mExonDownstream - mExonUpstream) == 1;
    }

    public boolean isPrePromotor()
    {
        if(!isPromoter())
            return false;

        if(mGene.strand() == 1 && mTranscriptStart > 0 && svPosition() < mTranscriptStart - PROMOTOR_REGION_MAX)
            return true;
        else if(mGene.strand() == -1 && mTranscriptEnd > 0 && svPosition() > mTranscriptEnd + PROMOTOR_REGION_MAX)
            return true;
        else
            return false;
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
            return TRANS_REGION_TYPE_PROMOTOR;

        return "Unknown";
    }

    private final String calcCodingType()
    {
        if(mTotalCodingBases == 0)
            return TRANS_CODING_TYPE_NON_CODING;
        else if(mCodingBases == 0)
            return TRANS_CODING_TYPE_UPSTREAM;
        else if(mCodingBases == mTotalCodingBases)
            return TRANS_CODING_TYPE_DOWNSTREAM;
        else
            return TRANS_CODING_TYPE_CODING;
    }

    public boolean isCanonical() { return mCanonical; }

    @NotNull
    public GeneAnnotation parent() { return mGene; }

    public long svPosition() { return mGene.position(); }

    @NotNull
    public String geneName() { return mGene.geneName(); }

    public int exonUpstream() { return mExonUpstream; }
    public int exonUpstreamPhase() { return mExonUpstreamPhase; }

    public long codingBases() { return mCodingBases; }
    public long totalCodingBases() { return mTotalCodingBases; }

    public int exonDownstream() { return mExonDownstream; }
    public int exonDownstreamPhase() { return mExonDownstreamPhase; }

    public int exonMax() { return mExonMax; }

    public long transcriptStart() { return mTranscriptStart; }
    public long transcriptEnd() { return mTranscriptEnd; }

    @Nullable
    public Long codingStart() { return mCodingStart; }

    @Nullable
    public Long codingEnd() { return mCodingEnd; }

    public final String toString()
    {
        return mGene.geneName() + " " + mTranscriptId;
    }
}
