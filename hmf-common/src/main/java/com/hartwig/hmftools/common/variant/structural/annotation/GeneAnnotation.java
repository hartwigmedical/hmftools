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

    private EnrichedStructuralVariant mVariant; // CHSH: try to remove this connection

    private final List<Transcript> transcripts = Lists.newArrayList();
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

    private String mPrecedingGeneId;

    public GeneAnnotation(int varId, final boolean isStart, final String geneName, final String stableId,
            final int strand, final List<String> synonyms, final List<Integer> entrezIds,
            final String karyotypeBand)
    {
        GeneName = geneName;
        StableId = stableId;
        Strand = strand;

        mVarId = varId;
        mIsStart = isStart;
        mVariant = null;

        mChromosome = "";
        mOrientation = 0;
        mPosition = -1;
        mPloidy = 0;
        mInsertSequence = "";
        mPrecedingGeneId = "";

        mSynonyms = synonyms;
        mEntrezIds = entrezIds;
        mKaryotypeBand = karyotypeBand;
    }

    public GeneAnnotation(@NotNull final EnrichedStructuralVariant variant, final boolean isStart, final String geneName,
            final String stableId, final int strand, final List<String> synonyms, final List<Integer> entrezIds, final String karyotypeBand)
    {
        GeneName = geneName;
        StableId = stableId;
        Strand = strand;

        mVariant = variant;
        mIsStart = isStart;

        mChromosome = "";
        mOrientation = 0;
        mPosition = -1;
        mPloidy = 0;
        mVarId = -1;
        mInsertSequence = "";
        mPrecedingGeneId = "";

        if(variant != null)
        {
            mVarId = variant.primaryKey() != null ? variant.primaryKey() : 0;
            mOrientation = variant.orientation(mIsStart);
            mPloidy = variant.ploidy() != null ? variant.ploidy() : 0;
            mPosition = variant.position(mIsStart);
            mChromosome = variant.chromosome(mIsStart);
            mSvType = variant.type();
            mInsertSequence = variant.insertSequence();
        }

        mSynonyms = synonyms;
        mEntrezIds = entrezIds;
        mKaryotypeBand = karyotypeBand;
    }

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
        mVariant = var;

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

    public final EnrichedStructuralVariant variant() { return mVariant; }

    public boolean isStart() { return mIsStart; }
    public boolean isEnd() { return !mIsStart; }

    public void addTranscript(@NotNull Transcript transcript) {
        transcripts.add(transcript);
    }

    @NotNull
    public List<Transcript> transcripts() {
        return ImmutableList.copyOf(transcripts);
    }

    @Nullable
    public Transcript canonical() {
        return transcripts.stream().filter(Transcript::isCanonical).findFirst().orElse(null);
    }

    @NotNull
    public List<String> synonyms() {
        return ImmutableList.copyOf(mSynonyms);
    }

    @NotNull
    public List<Integer> entrezIds() {
        return mEntrezIds;
    }

    @NotNull
    public String karyotypeBand() {
        return mKaryotypeBand;
    }

    public final String getmPrecedingGeneId() { return mPrecedingGeneId; }
    public void setPrecedingGeneId(final String geneId) { mPrecedingGeneId = geneId; }

    public static boolean isUpstream(final GeneAnnotation gene)
    {
        return gene.Strand * gene.orientation() > 0;
    }

    public static boolean isDownstream(final GeneAnnotation gene)
    {
        return !isUpstream(gene);
    }

}
