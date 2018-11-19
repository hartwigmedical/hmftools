package com.hartwig.hmftools.bachelorpp.types;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantConsequence;

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
    private List<String> mEffectsList;
    private String mAnnotations;
    private int mPhredScore;
    private boolean mIsHomozygous;
    private String mHgvsProtein;
    private String mHgvsCoding;
    private String mMatchType;

    private int mRefCount;
    private int mAltCount;
    private boolean mReadDataSet;

    private double mAdjustedVaf;

    private SomaticVariant mSomaticVariant;
    private VariantContext mVariantContext;
    private EnrichedSomaticVariant mEnrichedVariant;

    public static int PHRED_SCORE_CUTOFF = 150;

    public BachelorGermlineVariant(String patient, String source, String program, String varId,
            String gene, String transcriptId, String chromosome, long position,
            String ref, String alts, String effects, String annotations, String hgvsProtein,
            boolean isHomozygous, int phredScore, String hgvsCoding, String matchType)
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
        mAnnotations = annotations;
        mPhredScore = phredScore;
        mIsHomozygous = isHomozygous;
        mHgvsProtein = hgvsProtein;
        mHgvsCoding = hgvsCoding;
        mMatchType = matchType;

        mRef = mRef.replaceAll("\\*", "");
        mAlts = mAlts.replaceAll("\\*", "");

        mEffects = effects;
        mEffectsList = Arrays.stream(effects.split("&")).collect(Collectors.toList());

        mRefCount = 0;
        mAltCount = 0;
        mReadDataSet = false;
        mAdjustedVaf = 0;

        mSomaticVariant = null;
        mVariantContext = null;
        mEnrichedVariant = null;
    }

    public final String variantId() { return mVariantId; };
    public final String patient() { return mPatient; };
    public final String source() { return mSource; };
    public final String program() { return mProgram; };
    public final String gene() { return mGene; };
    public final String transcriptId() { return mTranscriptId; };
    public final String chromosome() { return mChromosome; };
    public long position() { return mPosition; };
    public final String ref() { return mRef; };
    public final String alts() { return mAlts; };
    public final String effects() { return mEffects; };
    public final List<String> effectsList() { return mEffectsList; }
    public final String annotations() { return mAnnotations; };
    public final String hgvsProtein() { return mHgvsProtein; };
    public final String hgvsCoding() { return mHgvsCoding; };
    public boolean isHomozygous() { return mIsHomozygous; }
    public String matchType() { return mMatchType; }
    public int getRefCount() { return mRefCount; }
    public int getAltCount() { return mAltCount; }
    public int getReadDepth() { return mAltCount + mRefCount; }
    public void setRefCount(int count) { mRefCount = count; }
    public void setAltCount(int count) { mAltCount = count; }

    public void setAltReadData(int altCount, int readDepth)
    {
        mAltCount = altCount;
        mRefCount = readDepth - altCount;
        mReadDataSet = true;
    }

    public boolean hasEffect(final VariantConsequence consequence)
    {
        for(final String effect : mEffectsList)
        {
            if(consequence.isParentTypeOf(effect))
                return true;
        }

        return false;
    }

    public boolean isReadDataSet() { return mReadDataSet; }

    public void setAdjustedVaf(double vaf) { mAdjustedVaf = vaf; }
    public double getAdjustedVaf() { return mAdjustedVaf; }

    public boolean isBiallelic()
    {
        if(mEnrichedVariant == null)
            return false;

        double copyNumber = mEnrichedVariant.adjustedCopyNumber();
        double minorAllelePloidy = copyNumber - (copyNumber * mAdjustedVaf);
        return (minorAllelePloidy < 0.5);
    }

    public boolean isValid()
    {
        return mSomaticVariant != null && mEnrichedVariant != null && mVariantContext != null;
    }

    public int phredScore()
    {
        return mPhredScore;
    }
    public boolean isLowScore()
    {
        return mPhredScore < PHRED_SCORE_CUTOFF && mAdjustedVaf < 0;
    }

    public final SomaticVariant getSomaticVariant() { return mSomaticVariant; }
    public final EnrichedSomaticVariant getEnrichedVariant() { return mEnrichedVariant; }

    public void setSomaticVariant(final SomaticVariant var) { mSomaticVariant = var; }
    public void setVariantContext(final VariantContext var) { mVariantContext = var; }
    public void setEnrichedVariant(final EnrichedSomaticVariant var) { mEnrichedVariant = var; }
}
