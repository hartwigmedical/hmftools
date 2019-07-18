package com.hartwig.hmftools.common.variant.structural.annotation;

import static com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation.isDownstream;

import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class Transcript {

    public final int TransId;
    public final String StableId;

    @Nullable
    public final Long CodingStart;

    @Nullable
    public final Long CodingEnd;

    public final long TranscriptStart;
    public final long TranscriptEnd;

    public final int ExonUpstream;
    public final int ExonUpstreamPhase;
    public final int ExonDownstream;
    public final int ExonDownstreamPhase;
    public final int ExonMax;

    @NotNull
    private final GeneAnnotation mGene;

    private int mExactCodingBase;
    private final int mCodingBases;
    private final int mTotalCodingBases;

    private final boolean mCanonical;
    private String mBioType;

    private String mCodingType;
    private final String mRegionType;

    private int mExonDistanceUp;
    private int mExonDistanceDown;

    private Map<Integer,Integer> mAlternativePhasing;

    private boolean mIsDisruptive;
    private boolean mReportableDisruption;

    private String mProteinFeaturesKept;
    private String mProteinFeaturesLost;

    public static final String TRANS_REGION_TYPE_UPSTREAM = "Upstream"; // promotor and earlier
    public static final String TRANS_REGION_TYPE_EXONIC = "Exonic";
    public static final String TRANS_REGION_TYPE_INTRONIC = "Intronic";

    public static final String TRANS_CODING_TYPE_CODING = "Coding";
    public static final String TRANS_CODING_TYPE_5P_UTR = "5P_UTR";
    public static final String TRANS_CODING_TYPE_3P_UTR = "3P_UTR";
    public static final String TRANS_CODING_TYPE_NON_CODING = "NonCoding";

    private static final int STOP_CODON_LENGTH = 3;

    public Transcript(@NotNull final GeneAnnotation parent, int transId, final String stableId,
            final int exonUpstream, final int exonUpstreamPhase, final int exonDownstream, final int exonDownstreamPhase,
            final int codingBases, final int totalCodingBases,
            final int exonMax, final boolean canonical, long transcriptStart, long transcriptEnd,
            final Long codingStart, final Long codingEnd)
    {
        TransId = transId;
        StableId = stableId;
        mCanonical = canonical;
        CodingStart = codingStart;
        CodingEnd = codingEnd;
        TranscriptStart = transcriptStart;
        TranscriptEnd = transcriptEnd;
        ExonDownstreamPhase = exonDownstreamPhase;
        ExonUpstreamPhase = exonUpstreamPhase;
        ExonUpstream = exonUpstream;
        ExonDownstream = exonDownstream;
        ExonMax = exonMax;

        mGene = parent;

        mExactCodingBase = -1;

        mBioType = "";
        mExonDistanceUp = 0;
        mExonDistanceDown = 0;
        mAlternativePhasing = Maps.newHashMap();

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

        mCodingType = calcCodingType();
        mRegionType = calcRegionType();

        if(isDownstream(mGene) && mRegionType == TRANS_REGION_TYPE_UPSTREAM)
            mIsDisruptive = false;
        else
            mIsDisruptive = true;

        mReportableDisruption = false;

        mProteinFeaturesKept = "";
        mProteinFeaturesLost = "";
    }

    public boolean isExonic()
    {
        return ExonUpstream > 0 && ExonUpstream == ExonDownstream;
    }

    public boolean isPromoter()
    {
        return ExonUpstream == 0 && (ExonDownstream == 1 || ExonDownstream == 2);
    }

    public boolean isIntronic()
    {
        return ExonUpstream > 0 && (ExonDownstream - ExonUpstream) == 1;
    }

    public long getDistanceUpstream()
    {
        if(!isPromoter())
            return 0;

        if(mGene.Strand == 1)
            return TranscriptStart - svPosition();
        else
            return svPosition() - TranscriptEnd;
    }

    public final String codingType() { return mCodingType; }
    public final String regionType() { return mRegionType; }

    // for convenience
    public boolean isCoding() { return mCodingType.equals(TRANS_CODING_TYPE_CODING); }

    public boolean preCoding()
    {
        return mCodingType.equals(TRANS_CODING_TYPE_5P_UTR);
    }

    public boolean postCoding()
    {
        return mCodingType.equals(TRANS_CODING_TYPE_3P_UTR);
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
        if(CodingStart == null || CodingEnd == null || mTotalCodingBases == 0)
        {
            return TRANS_CODING_TYPE_NON_CODING;
        }
        else if(mCodingBases == 0)
        {
            return TRANS_CODING_TYPE_5P_UTR;
        }
        else if(mCodingBases == mTotalCodingBases)
        {
            return TRANS_CODING_TYPE_3P_UTR;
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

    public final Map<Integer,Integer> getAlternativePhasing() { return mAlternativePhasing; }
    public void setAlternativePhasing(final Map<Integer,Integer> phasings) { mAlternativePhasing = phasings; }

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
    public String geneName() { return mGene.GeneName; }
    public boolean isUpstream() { return mGene.isUpstream(); }

    public int codingBases() { return mCodingBases; }
    public int calcCodingBases(boolean isUpstream) { return isUpstream ? mCodingBases : mTotalCodingBases - mCodingBases; }
    public int totalCodingBases() { return mTotalCodingBases; }

    public void setExonicCodingBase()
    {
        int calcStartPhase = calcPositionPhasing(this, isUpstream());
        setExactCodingBase(calcStartPhase);
    }

    public static int calcPositionPhasing(final Transcript transcript, boolean isUpstream)
    {
        // if upstream then can just use the coding bases
        // if downstream then coding bases are what's remaining
        long codingBases = transcript.calcCodingBases(isUpstream);

        // factor in insert sequence for the upstream partner
        if(isUpstream && !transcript.parent().insertSequence().isEmpty())
        {
            codingBases += transcript.parent().insertSequence().length();
        }

        int adjustedPhase = (int)(codingBases % 3);

        return adjustedPhase;
    }

    public int exactCodingBase() { return mExactCodingBase; }
    public void setExactCodingBase(int base) { mExactCodingBase = base; }

    public final long length() { return TranscriptEnd - TranscriptStart; }

    public boolean isDisruptive() { return mIsDisruptive; }
    public void setIsDisruptive(boolean toggle) { mIsDisruptive = toggle; }

    public boolean reportableDisruption() { return mReportableDisruption; }
    public void setReportableDisruption(boolean toggle) { mReportableDisruption = toggle; }


    public long codingStart() { return CodingStart != null ? CodingStart : 0; }
    public long codingEnd() { return CodingEnd != null ? CodingEnd : 0; }

    public final String toString()
    {
        return mGene.GeneName + " " + StableId;
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
