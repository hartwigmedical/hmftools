package com.hartwig.hmftools.common.fusion;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

// class linking an SV breakend to a potentially impacted gene
public class BreakendGeneData
{
    public final String GeneName;
    public final String StableId;
    public final byte Strand;

    private final int mVarId;
    private final boolean mIsStart;
    private boolean mUpstream;

    private GeneData mGeneData;

    private final List<BreakendTransData> mTranscripts;

    @NotNull
    private final String mKaryotypeBand;

    private StructuralVariantType mSvType;
    private String mChromosome;
    private byte mOrientation;
    private int mPosition;
    private double mJunctionCopyNumber;
    private String mInsertSequence;

    public BreakendGeneData(int varId, final boolean isStart, final String geneName, final String stableId,
            final byte strand, @NotNull String karyotypeBand)
    {
        GeneName = geneName;
        StableId = stableId;
        Strand = strand;

        mTranscripts = Lists.newArrayList();

        mVarId = varId;
        mIsStart = isStart;
        mGeneData = null;

        mChromosome = "";
        mOrientation = 0;
        mPosition = -1;
        mJunctionCopyNumber = 0;
        mInsertSequence = "";

        mKaryotypeBand = karyotypeBand;
    }

    public void setGeneData(final GeneData geneData) { mGeneData = geneData; }
    public final GeneData getGeneData() { return mGeneData; }

    public void setPositionalData(final String chromosome, int position, byte orientation)
    {
        mChromosome = chromosome;
        mPosition = position;
        mOrientation = orientation;
        mUpstream = isUpstream(this);
    }

    public void setSvData(final StructuralVariantData var, double jcn)
    {
        mOrientation = mIsStart ? var.startOrientation() : var.endOrientation();
        mJunctionCopyNumber = jcn;
        mPosition = mIsStart ? var.startPosition() : var.endPosition();
        mChromosome = mIsStart ? var.startChromosome() : var.endChromosome();
        mSvType = var.type();
        mInsertSequence = var.insertSequence();
        mUpstream = isUpstream(this);
    }

    public void setType(StructuralVariantType type) { mSvType = type; }

    public int id() { return mVarId; }
    public byte orientation() { return mOrientation; }
    public int position() { return mPosition; }
    public StructuralVariantType type() { return mSvType; }
    public String chromosome() { return mChromosome; }
    public double jcn() { return mJunctionCopyNumber; }
    public String insertSequence() { return mInsertSequence; }
    public boolean isUpstream() { return mUpstream; }

    public boolean isStart() { return mIsStart; }
    public boolean isEnd() { return !mIsStart; }

    public void addTranscript(BreakendTransData transcript) {
        mTranscripts.add(transcript);
    }

    public List<BreakendTransData> transcripts() { return mTranscripts; }

    @Nullable
    public BreakendTransData canonical()
    {
         return mTranscripts.stream().filter(BreakendTransData::isCanonical).findFirst().orElse(null);
    }

    public String karyotypeBand() { return mKaryotypeBand; }

    public static boolean isUpstream(final BreakendGeneData gene)
    {
        return gene.Strand * gene.orientation() > 0;
    }

    public static boolean isDownstream(final BreakendGeneData gene)
    {
        return !isUpstream(gene);
    }

    public boolean breakendWithinGene(int preGeneDistance)
    {
        if(mGeneData == null)
            return false;

        if(mPosition >= mGeneData.GeneStart && mPosition <= mGeneData.GeneEnd)
            return true;

        if(preGeneDistance <= 0)
            return false;

        // return true if the gene has a transcript such that the breakend falls into its pre-gene region
        for(final BreakendTransData trans : mTranscripts)
        {
            // exclude if the position is interrupted by another splice acceptor
            if(!trans.isUpstream() && trans.hasNegativePrevSpliceAcceptorDistance())
                continue;

            int distance = Strand == 1 ? trans.TransData.TransStart - mPosition : mPosition - trans.TransData.TransEnd;

            if(distance > 0 && distance <= preGeneDistance)
                return true;
        }

        return false;
    }

    public String toString()
    {
        return String.format("gene(%s:%s strand=%d) breakend(sv=%d pos=%s:%d:%d) trans(%d)",
                mGeneData.GeneId, mGeneData.GeneName, mGeneData.Strand, mVarId, mChromosome, mPosition, mOrientation, mTranscripts.size());
    }

}
