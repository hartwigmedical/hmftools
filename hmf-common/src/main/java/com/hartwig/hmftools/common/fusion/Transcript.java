package com.hartwig.hmftools.common.fusion;

import static com.hartwig.hmftools.common.fusion.GeneAnnotation.isDownstream;

import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

// class linking an SV breakend to a potentially impacted transcript
public class Transcript {

    public final int TransId;
    public final String StableId;

    @Nullable
    public final Integer CodingStart;

    @Nullable
    public final Integer CodingEnd;

    public final int TranscriptStart;
    public final int TranscriptEnd;

    public final int ExonUpstream;
    public final int ExonDownstream;

    public final int ExonDownstreamPhase;
    public final int ExonUpstreamPhase;

    public final int ExonMax;

    @NotNull
    private final GeneAnnotation mGene;

    private final int mCodingBases; // number of bases into coding where this breakend occurs
    private final int mTotalCodingBases;
    private int mExonicBasePhase; // phase of base if in exon

    private final boolean mCanonical;
    private String mBioType;

    private String mCodingType;
    private final String mRegionType;

    private Integer mPrevSpliceAcceptorDistance;
    private Integer mNextSpliceAcceptorDistance;

    private Map<Integer,Integer> mAlternativePhasing;

    private boolean mIsDisruptive;
    private boolean mReportableDisruption;
    private double mUndisruptedCopyNumber;

    private String mProteinFeaturesKept;
    private String mProteinFeaturesLost;

    public static final String TRANS_REGION_TYPE_UPSTREAM = "Upstream"; // promotor and earlier
    public static final String TRANS_REGION_TYPE_EXONIC = "Exonic";
    public static final String TRANS_REGION_TYPE_INTRONIC = "Intronic";

    public static final String TRANS_CODING_TYPE_CODING = "Coding";
    public static final String TRANS_CODING_TYPE_5P_UTR = "5P_UTR";
    public static final String TRANS_CODING_TYPE_3P_UTR = "3P_UTR";
    public static final String TRANS_CODING_TYPE_NON_CODING = "NonCoding";

    public static final int POST_CODING_PHASE = -2;

    private static final int STOP_CODON_LENGTH = 3;

    public static final int NO_NEXT_SPLICE_ACCEPTOR = -1;

    public Transcript(@NotNull final GeneAnnotation gene, int transId, final String stableId,
            final int exonUpstream, final int exonUpstreamPhase, final int exonDownstream, final int exonDownstreamPhase,
            final int codingBases, final int totalCodingBases,
            final int exonMax, final boolean canonical, int transcriptStart, int transcriptEnd,
            final Integer codingStart, final Integer codingEnd)
    {
        TransId = transId;
        StableId = stableId;
        mCanonical = canonical;
        CodingStart = codingStart;
        CodingEnd = codingEnd;
        TranscriptStart = transcriptStart;
        TranscriptEnd = transcriptEnd;
        ExonUpstream = exonUpstream;
        ExonDownstream = exonDownstream;
        ExonMax = exonMax;

        mGene = gene;

        mExonicBasePhase = -1;

        mBioType = "";
        mPrevSpliceAcceptorDistance = null;
        mNextSpliceAcceptorDistance = null;

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

        if(mCodingType == TRANS_CODING_TYPE_3P_UTR)
        {
            ExonUpstreamPhase = POST_CODING_PHASE;
            ExonDownstreamPhase = POST_CODING_PHASE;
        }
        else
        {
            ExonDownstreamPhase = exonDownstreamPhase;
            ExonUpstreamPhase = exonUpstreamPhase;
        }

        if(isDownstream(mGene) && mRegionType == TRANS_REGION_TYPE_UPSTREAM)
            mIsDisruptive = false;
        else
            mIsDisruptive = true;

        mReportableDisruption = false;
        mUndisruptedCopyNumber = 0;

        mProteinFeaturesKept = "";
        mProteinFeaturesLost = "";
    }

    @NotNull
    public GeneAnnotation gene() { return mGene; }

    public int svPosition() { return mGene.position(); }
    public String geneName() { return mGene.GeneName; }
    public boolean isUpstream() { return mGene.isUpstream(); }

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

    public int getDistanceUpstream()
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

    private String calcRegionType()
    {
        if(isIntronic())
            return TRANS_REGION_TYPE_INTRONIC;

        if(isExonic())
            return TRANS_REGION_TYPE_EXONIC;

        if(isPromoter())
            return TRANS_REGION_TYPE_UPSTREAM;

        return "Unknown";
    }

    public int nextSpliceExonRank()
    {
        if(isExonic())
            return isUpstream() ? ExonUpstream - 1 : ExonDownstream + 1;
        else
            return isUpstream() ? ExonUpstream : ExonDownstream;
    }

    public int nextSpliceExonPhase()
    {
        return isUpstream() ? ExonUpstreamPhase : ExonDownstreamPhase;
    }

    private String calcCodingType()
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

    public void setSpliceAcceptorDistance(boolean isPrevious, Integer distance)
    {
        if(isPrevious)
            mPrevSpliceAcceptorDistance = distance;
        else
            mNextSpliceAcceptorDistance = distance;
    }

    public boolean hasPrevSpliceAcceptorDistance()
    {
        return mPrevSpliceAcceptorDistance != null;
    }

    public boolean hasNegativePrevSpliceAcceptorDistance()
    {
        return mPrevSpliceAcceptorDistance != null && mPrevSpliceAcceptorDistance < 0;
    }

    public int prevSpliceAcceptorDistance()
    {
        return mPrevSpliceAcceptorDistance != null ? mPrevSpliceAcceptorDistance : NO_NEXT_SPLICE_ACCEPTOR;
    }

    public int nextSpliceAcceptorDistance()
    {
        return mNextSpliceAcceptorDistance != null ? mNextSpliceAcceptorDistance : NO_NEXT_SPLICE_ACCEPTOR;
    }

    public void setBioType(final String type) { mBioType = type; }
    public final String bioType() { return mBioType; }

    public int codingBases() { return mCodingBases; }

    public int calcCodingBases()
    {
        // returns number of coding bases preserved in the context of this breakend and whether up or down stream
        return mGene.isUpstream() ? mCodingBases : mTotalCodingBases - mCodingBases;
    }

    public int totalCodingBases() { return mTotalCodingBases; }

    public void setExonicCodingBase()
    {
        mExonicBasePhase = calcPositionPhasing(this, isUpstream());
    }

    public static int calcPositionPhasing(final Transcript transcript, boolean isUpstream)
    {
        // int codingBases = transcript.calcCodingBases(isUpstream);
        int codingBases = transcript.codingBases();

        // subtract 1 to get back to phasing starting at zero, ie start with the first coding base where position == coding start:
        // coding base: 1 2 3 4 5 6 7 8 9 10
        // phase:       0 1 2 0 1 2 0 1 2 0
        codingBases -= 1;

        // factor in insert sequence for the upstream partner
        if(isUpstream && !transcript.gene().insertSequence().isEmpty())
        {
            codingBases += transcript.gene().insertSequence().length();
        }

        int adjustedPhase = (int)(codingBases % 3);

        return adjustedPhase;
    }

    public int exonicBasePhase() { return mExonicBasePhase; }

    public final int length() { return TranscriptEnd - TranscriptStart; }

    public boolean isDisruptive() { return mIsDisruptive; }
    public void setIsDisruptive(boolean toggle) { mIsDisruptive = toggle; }

    public boolean reportableDisruption() { return mReportableDisruption; }
    public void setReportableDisruption(boolean toggle) { mReportableDisruption = toggle; }

    public double undisruptedCopyNumber() { return mUndisruptedCopyNumber; }
    public void setUndisruptedCopyNumber(double copyNumber) { mUndisruptedCopyNumber = copyNumber; }

    public int codingStart() { return CodingStart != null ? CodingStart : 0; }
    public int codingEnd() { return CodingEnd != null ? CodingEnd : 0; }

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
