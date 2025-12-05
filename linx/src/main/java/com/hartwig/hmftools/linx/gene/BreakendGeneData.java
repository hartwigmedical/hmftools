package com.hartwig.hmftools.linx.gene;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.jetbrains.annotations.Nullable;

// class linking an SV breakend to a potentially impacted gene
public class BreakendGeneData
{
    public final GeneData GeneData;

    private final int mSvId;
    private String mSvVcfId;
    private final boolean mIsStart;
    private boolean mUpstream;

    private final List<BreakendTransData> mTranscripts;

    private StructuralVariantType mSvType;
    private String mChromosome;
    private byte mOrientation;
    private int mPosition;
    private double mJunctionCopyNumber;
    private String mInsertSequence;

    public BreakendGeneData(final int svId, final boolean isStart, final GeneData geneData)
    {
        GeneData = geneData;

        mTranscripts = Lists.newArrayList();

        mSvId = svId;
        mSvVcfId = "";
        mIsStart = isStart;

        mChromosome = "";
        mOrientation = 0;
        mPosition = -1;
        mJunctionCopyNumber = 0;
        mInsertSequence = "";
    }

    public String geneName() { return GeneData.GeneName; }
    public String geneId() { return GeneData.GeneId; }
    public byte strand() { return GeneData.Strand; }

    public void setPositionalData(final String chromosome, int position, byte orientation)
    {
        mChromosome = chromosome;
        mPosition = position;
        mOrientation = orientation;
        mUpstream = isUpstream(this);
    }

    public void setSvData(final StructuralVariantData var, double jcn)
    {
        mSvVcfId = mIsStart ? var.vcfIdStart() : var.vcfIdEnd();
        mOrientation = mIsStart ? var.startOrientation() : var.endOrientation();
        mJunctionCopyNumber = jcn;
        mPosition = mIsStart ? var.startPosition() : var.endPosition();
        mChromosome = mIsStart ? var.startChromosome() : var.endChromosome();
        mSvType = var.type();
        mInsertSequence = var.insertSequence();
        mUpstream = isUpstream(this);
    }

    public void setType(StructuralVariantType type) { mSvType = type; }

    public int varId() { return mSvId; }
    public String vcfId() { return mSvVcfId; }
    public byte orientation() { return mOrientation; }
    public int position() { return mPosition; }
    public StructuralVariantType svType() { return mSvType; }
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

    public static boolean isUpstream(final BreakendGeneData gene)
    {
        return gene.strand() * gene.orientation() > 0;
    }

    public static boolean isDownstream(final BreakendGeneData gene)
    {
        return !isUpstream(gene);
    }

    public boolean breakendWithinGene(int preGeneDistance)
    {
        if(GeneData == null)
            return false;

        if(mPosition >= GeneData.GeneStart && mPosition <= GeneData.GeneEnd)
            return true;

        if(preGeneDistance <= 0)
            return false;

        // return true if the gene has a transcript such that the breakend falls into its pre-gene region
        for(final BreakendTransData trans : mTranscripts)
        {
            // exclude if the position is interrupted by another splice acceptor
            if(!trans.isUpstream() && trans.hasNegativePrevSpliceAcceptorDistance())
                continue;

            int distance = GeneData.Strand == POS_STRAND ? trans.TransData.TransStart - mPosition : mPosition - trans.TransData.TransEnd;

            if(distance > 0 && distance <= preGeneDistance)
                return true;
        }

        return false;
    }

    public String coordsStr() { return SvVarData.coordsStr(mChromosome, mPosition, mOrientation); }

    public String toString()
    {
        return String.format("gene(%s:%s strand=%d) breakend(sv=%d coords=%s) trans(%d)",
                GeneData.GeneId, GeneData.GeneName, GeneData.Strand, varId(), coordsStr(), mTranscripts.size());
    }
}
