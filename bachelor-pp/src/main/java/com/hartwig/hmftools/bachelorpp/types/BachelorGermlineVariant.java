package com.hartwig.hmftools.bachelorpp.types;

import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import htsjdk.variant.variantcontext.VariantContext;

public class BachelorGermlineVariant {

    private String mPatient;
    private String mSource;
    private String mProgram;
    private String mVariantId;
    private String mGene;
    private String mTranscriptId;
    private String mChromosome;
    private long mPosition;
    private String mRef;
    private String mAlts;
    private String mEffects;

    private int mRefCount;
    private int mAltCount;

    private SomaticVariant mSomaticVariant;
    private VariantContext mVariantContext;
    private EnrichedSomaticVariant mEnrichedVariant;

    public BachelorGermlineVariant(String patient, String source, String program, String varId,
            String gene, String transcriptId, String chromosome, long position,
            String ref, String alts, String effects)
    {
        mPatient = patient;
        mSource = source;
        mProgram = program;
        mVariantId = varId;
        mGene = gene;
        mTranscriptId = transcriptId;
        mChromosome = chromosome;
        mPosition = position;
        mRef = ref;
        mAlts = alts;

        mRef = mRef.replaceAll("\\*", "");
        mAlts = mAlts.replaceAll("\\*", "");

        mEffects = effects;
        mRefCount = 0;
        mAltCount = 0;

        mSomaticVariant = null;
        mVariantContext = null;
        mEnrichedVariant = null;
    }

    public String variantId() { return mVariantId; };
    public String patient() { return mPatient; };
    public String source() { return mSource; };
    public String program() { return mProgram; };
    public String gene() { return mGene; };
    public String transcriptId() { return mTranscriptId; };
    public String chromosome() { return mChromosome; };
    public long position() { return mPosition; };
    public String ref() { return mRef; };
    public String alts() { return mAlts; };
    public String effects() { return mEffects; };
    public int getRefCount() { return mRefCount; }
    public int getAltCount() { return mAltCount; }
    public void setRefCount(int count) { mRefCount = count; }
    public void setAltCount(int count) { mAltCount = count; }

    public final SomaticVariant getSomaticVariant() { return mSomaticVariant; }
    public final VariantContext getVariantContext() { return mVariantContext; }
    public final EnrichedSomaticVariant getEnrichedVariant() { return mEnrichedVariant; }

    public void setSomaticVariant(final SomaticVariant var) { mSomaticVariant = var; }
    public void setVariantContext(final VariantContext var) { mVariantContext = var; }
    public void setEnrichedVariant(final EnrichedSomaticVariant var) { mEnrichedVariant = var; }
}
