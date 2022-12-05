package com.hartwig.hmftools.linx.gene;

import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.NON_CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.DOWNSTREAM;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.UNKNOWN;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.UPSTREAM;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.tickPhaseForward;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.TranscriptCodingType;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;

import org.jetbrains.annotations.NotNull;

// class linking an SV breakend to a potentially impacted transcript
public class BreakendTransData
{
    public final TranscriptData TransData;

    public final int ExonUpstream;
    public final int ExonDownstream;

    @NotNull
    private final BreakendGeneData mGene;

    public final int CodingBases; // number of bases into coding where this breakend occurs
    public final int TotalCodingBases;

    public int Phase; // intronic phase or the start phase for upstream, exon end for downstream
    public int ExonicBasePhase; // phase of base if in exon

    private TranscriptCodingType mCodingType;
    private TranscriptRegionType mRegionType;

    private Integer mPrevSpliceAcceptorDistance;
    private Integer mNextSpliceAcceptorDistance;

    private final Map<Integer,Integer> mAlternativePhasing;

    private boolean mIsDisruptive;
    private boolean mReportableDisruption;
    private double mUndisruptedCopyNumber;
    private boolean mUndisruptedCopyNumberSet;

    private String mProteinFeaturesKept;
    private String mProteinFeaturesLost;

    public static final int POST_CODING_PHASE = -2;
    private static final int STOP_CODON_LENGTH = 3;
    public static final int NO_NEXT_SPLICE_ACCEPTOR = -1;

    public BreakendTransData(
            final BreakendGeneData gene, final TranscriptData transData,
            final int exonUpstream, final int exonDownstream, int phase, int exonicBasePhase, int codingBases, int totalCodingBases)
    {
        TransData = transData;

        ExonUpstream = exonUpstream;
        ExonDownstream = exonDownstream;

        mGene = gene;

        Phase = PHASE_NONE;

        mPrevSpliceAcceptorDistance = null;
        mNextSpliceAcceptorDistance = null;

        mAlternativePhasing = Maps.newHashMap();

        if(totalCodingBases > STOP_CODON_LENGTH)
        {
            // remove the stop codon from what is consider coding
            TotalCodingBases = totalCodingBases - STOP_CODON_LENGTH;
        }
        else
        {
            TotalCodingBases = 0;
        }

        // factor in insert sequence for the upstream partner if exonic
        if(gene.isUpstream() && !gene.insertSequence().isEmpty() && ExonUpstream == ExonDownstream
        && codingBases > 0 && codingBases < totalCodingBases)
        {
            int insSeqLength = gene.insertSequence().length();
            codingBases += insSeqLength;
            exonicBasePhase = tickPhaseForward(phase, insSeqLength);
        }

        ExonicBasePhase = exonicBasePhase;

        if(codingBases > TotalCodingBases)
            CodingBases = TotalCodingBases;
        else
            CodingBases = codingBases;

        setCodingType();
        setRegionType();

        if(mCodingType == UTR_3P)
        {
            Phase = POST_CODING_PHASE;
        }
        else
        {
            Phase = phase;
        }

        if(mRegionType == EXONIC || mRegionType == INTRONIC)
            mIsDisruptive = true;
        else
            mIsDisruptive = false;

        mReportableDisruption = false;
        mUndisruptedCopyNumber = 0;
        mUndisruptedCopyNumberSet = false;

        mProteinFeaturesKept = "";
        mProteinFeaturesLost = "";
    }

    @NotNull
    public BreakendGeneData gene() { return mGene; }

    public int svPosition() { return mGene.position(); }
    public String geneName() { return mGene.geneName(); }
    public boolean isUpstream() { return mGene.isUpstream(); }

    // convenience
    public int transId() { return TransData.TransId; }
    public String transName() { return TransData.TransName; }
    public int transStart() { return TransData.TransStart; }
    public int transEnd() { return TransData.TransEnd; }
    public final String bioType() { return TransData.BioType; }
    public int codingStart() { return TransData.CodingStart != null ? TransData.CodingStart : 0; }
    public int codingEnd() { return TransData.CodingEnd != null ? TransData.CodingEnd : 0; }
    public int exonCount() { return TransData.exons().size(); }

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

    public boolean isPostTranscript()
    {
        if(mGene.strand() == POS_STRAND)
            return svPosition() > TransData.TransEnd;
        else
            return svPosition() < TransData.TransStart;
    }

    public int getDistanceUpstream()
    {
        if(!isPromoter())
            return 0;

        if(mGene.strand() == POS_STRAND)
            return TransData.TransStart - svPosition();
        else
            return svPosition() - TransData.TransEnd;
    }

    public final TranscriptCodingType codingType() { return mCodingType; }
    public final TranscriptRegionType regionType() { return mRegionType; }

    // for convenience
    public boolean isCoding() { return mCodingType == CODING; }
    public boolean preCoding() { return mCodingType == UTR_5P; }
    public boolean postCoding() { return mCodingType == UTR_3P; }

    public boolean nonCoding()
    {
        return mCodingType == NON_CODING;
    }

    private void setRegionType()
    {
        if(isPromoter())
        {
            mRegionType = UPSTREAM;
        }
        else if(isPostTranscript())
        {
            mRegionType = DOWNSTREAM;
        }
        else if(isIntronic())
        {
            mRegionType = INTRONIC;
        }
        else if(isExonic())
        {
            mRegionType = EXONIC;
        }
        else
        {
            mRegionType = UNKNOWN;
        }
    }

    public void setRegionType(final TranscriptRegionType type) { mRegionType = type; }

    public int nextSpliceExonRank()
    {
        if(isExonic())
            return isUpstream() ? ExonUpstream - 1 : ExonDownstream + 1;
        else
            return isUpstream() ? ExonUpstream : ExonDownstream;
    }

    private void setCodingType()
    {
        if(TransData.CodingStart == null || TransData.CodingEnd == null)
        {
            mCodingType = NON_CODING;
            return;
        }

        int position = gene().position();

        if(positionWithin(position, TransData.CodingStart, TransData.CodingEnd))
        {
            mCodingType = CODING;
            return;
        }

        if(TransData.Strand == POS_STRAND)
        {
            mCodingType =  position > TransData.CodingEnd ? UTR_3P : UTR_5P;
        }
        else
        {
            mCodingType =  position < TransData.CodingStart ? UTR_3P : UTR_5P;
        }
    }

    public void setCodingType(final TranscriptCodingType type) { mCodingType = type; }

    public boolean isCanonical() { return TransData.IsCanonical; }

    public Map<Integer,Integer> getAlternativePhasing() { return mAlternativePhasing; }

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

    public boolean isDisruptive() { return mIsDisruptive; }
    public void setIsDisruptive(boolean toggle) { mIsDisruptive = toggle; }

    public boolean reportableDisruption() { return mReportableDisruption; }
    public void setReportableDisruption(boolean toggle) { mReportableDisruption = toggle; }

    public double undisruptedCopyNumber() { return mUndisruptedCopyNumber; }
    public boolean undisruptedCopyNumberSet() { return mUndisruptedCopyNumberSet; }
    public void setUndisruptedCopyNumber(double copyNumber)
    {
        mUndisruptedCopyNumber = copyNumber;
        mUndisruptedCopyNumberSet = true;
    }

    public final String toString()
    {
        return mGene.geneName() + " " + TransData.TransName;
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
