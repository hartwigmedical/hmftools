package com.hartwig.hmftools.common.variant.structural.annotation;

import java.util.List;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GeneAnnotation {

    public final String GeneName;
    public final String StableId;
    public final int Strand;

    private int mVarId;
    private final boolean mIsStart;
    private boolean mUpstream;

    private EnsemblGeneData mGeneData;

    private List<Transcript> mTranscripts;
    @NotNull
    private final List<String> mSynonyms;
    @NotNull
    private final String mKaryotypeBand;

    private StructuralVariantType mSvType;
    private String mChromosome;
    private byte mOrientation;
    private long mPosition;
    private double mPloidy;
    private String mInsertSequence;

    public GeneAnnotation(int varId, final boolean isStart, final String geneName, final String stableId,
            final int strand, final List<String> synonyms, final String karyotypeBand)
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
        mPloidy = 0;
        mInsertSequence = "";

        mSynonyms = synonyms;
        mKaryotypeBand = karyotypeBand;
    }

    public void setGeneData(final EnsemblGeneData geneData) { mGeneData = geneData; }
    public final EnsemblGeneData getGeneData() { return mGeneData; }

    public void setPositionalData(final String chromosome, long position, byte orientation)
    {
        mChromosome = chromosome;
        mPosition = position;
        mOrientation = orientation;
        mUpstream = isUpstream(this);
    }

    public void setSvData(final StructuralVariantData var)
    {
        mOrientation = mIsStart ? var.startOrientation() : var.endOrientation();
        mPloidy = var.ploidy();
        mPosition = mIsStart ? var.startPosition() : var.endPosition();
        mChromosome = mIsStart ? var.startChromosome() : var.endChromosome();
        mSvType = var.type();
        mInsertSequence = var.insertSequence();
        mUpstream = isUpstream(this);
    }

    public int id() { return mVarId; }
    public byte orientation() { return mOrientation; }
    public long position() { return mPosition; }
    public StructuralVariantType type() { return mSvType; }
    public String chromosome() { return mChromosome; }
    public double ploidy() { return mPloidy; }
    public String insertSequence() { return mInsertSequence; }
    public boolean isUpstream() { return mUpstream; }

    public boolean isStart() { return mIsStart; }
    public boolean isEnd() { return !mIsStart; }

    public void addTranscript(Transcript transcript) {
        mTranscripts.add(transcript);
    }

    public List<Transcript> transcripts() { return mTranscripts; }

    @Nullable
    public Transcript canonical() {
        return mTranscripts.stream().filter(Transcript::isCanonical).findFirst().orElse(null);
    }

    public List<String> synonyms() {
        return ImmutableList.copyOf(mSynonyms);
    }

    public String karyotypeBand() {
        return mKaryotypeBand;
    }

    public static boolean isUpstream(final GeneAnnotation gene)
    {
        return gene.Strand * gene.orientation() > 0;
    }

    public static boolean isDownstream(final GeneAnnotation gene)
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
        for(final Transcript trans : mTranscripts)
        {
            if(trans.prevSpliceAcceptorDistance() < 0) // exclude if the position is interupted by another splice acceptor
                continue;

            long distance = Strand == 1 ? trans.TranscriptStart - mPosition : mPosition - trans.TranscriptEnd;

            if(distance > 0 && distance <= preGeneDistance)
                return true;
        }

        return false;
    }

}
