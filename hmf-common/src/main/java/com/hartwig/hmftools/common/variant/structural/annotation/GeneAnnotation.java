package com.hartwig.hmftools.common.variant.structural.annotation;

import java.util.List;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
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

    private EnsemblGeneData mGeneData;

    private List<Transcript> mTranscripts;
    @NotNull
    private final List<String> mSynonyms;
    @NotNull
    private final List<Integer> mEntrezIds;
    @NotNull
    private final String mKaryotypeBand;

    private StructuralVariantType mSvType;
    private String mChromosome;
    private byte mOrientation;
    private long mPosition;
    private double mPloidy;
    private String mInsertSequence;

    public GeneAnnotation(int varId, final boolean isStart, final String geneName, final String stableId,
            final int strand, final List<String> synonyms, final List<Integer> entrezIds, final String karyotypeBand)
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
        mEntrezIds = entrezIds;
        mKaryotypeBand = karyotypeBand;
    }

    public void setGeneData(final EnsemblGeneData geneData) { mGeneData = geneData; }
    public final EnsemblGeneData getGeneData() { return mGeneData; }

    public void setPositionalData(final String chromosome, long position, byte orientation)
    {
        mChromosome = chromosome;
        mPosition = position;
        mOrientation = orientation;
    }

    public void setSvData(final StructuralVariantData var)
    {
        mOrientation = mIsStart ? var.startOrientation() : var.endOrientation();
        mPloidy = var.ploidy();
        mPosition = mIsStart ? var.startPosition() : var.endPosition();
        mChromosome = mIsStart ? var.startChromosome() : var.endChromosome();
        mSvType = var.type();
        mInsertSequence = var.insertSequence();
    }

    public void setSvData(final EnrichedStructuralVariant var)
    {
        if(var.end() == null && !mIsStart)
            return;

        mOrientation = var.orientation(mIsStart);
        mPloidy = var.ploidy() != null ? var.ploidy() : 0;
        mPosition = var.position(mIsStart);
        mChromosome = var.chromosome(mIsStart);
        mSvType = var.type();
        mInsertSequence = var.insertSequence();
    }

    public int id() { return mVarId; }
    public byte orientation() { return mOrientation; }
    public long position() { return mPosition; }
    public StructuralVariantType type() { return mSvType; }
    public String chromosome() { return mChromosome; }
    public double ploidy() { return mPloidy; }
    public String insertSequence() { return mInsertSequence; }

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

    public List<Integer> entrezIds() {
        return mEntrezIds;
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

}
