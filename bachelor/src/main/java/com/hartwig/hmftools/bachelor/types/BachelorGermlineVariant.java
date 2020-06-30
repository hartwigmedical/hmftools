package com.hartwig.hmftools.bachelor.types;

import static com.hartwig.hmftools.bachelor.types.FilterType.NONE;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class BachelorGermlineVariant implements Comparable<BachelorGermlineVariant>
{
    public final String SampleId;
    public final String Program;
    public final String VariantId;
    public final String Gene;
    public final String TranscriptId;
    public final String Chromosome;
    public final long Position;
    public final String Ref;
    public final String Alts;
    public final CodingEffect CodingEffect;
    public final String Effects;
    public final String Annotations;
    public final int PhredScore;
    public final boolean IsHomozygous;
    public final String HgvsProtein;
    public final String HgvsCoding;
    public final String CodonInfo;

    private FilterType mFilterType;
    private String mMatchType;

    private boolean mClinvarMatch;
    private String mClinvarSig;
    private String mClinvarSigInfo;

    private int mGermlineAltCount;
    private int mGermlineReadDepth;
    private int mTumorAltCount;
    private int mTumorReadDepth;
    private boolean mReadDataSet;

    private double mAdjustedVaf;

    private SomaticVariant mSomaticVariant;
    private VariantContext mVariantContext;
    private EnrichedSomaticVariant mEnrichedVariant;

    public static final String MATCH_TYPE_NONE = "None";
    public static final String MATCH_TYPE_REQUIRED_EFFECT = "RequiredEffect";
    public static final String MATCH_TYPE_WHITELIST = "WhiteList";

    public static final int PHRED_SCORE_CUTOFF = 150;

    public BachelorGermlineVariant(String sampleId, String program, String varId,
            String gene, String transcriptId, String chromosome, long position,
            String ref, String alts, CodingEffect codingEffect, String effects, String annotations, String hgvsProtein,
            boolean isHomozygous, int phredScore, String hgvsCoding, String codonInfo)
    {
        SampleId = sampleId;
        Program = program;
        VariantId = varId;
        Gene = gene;
        TranscriptId = transcriptId;
        Chromosome = chromosome;
        Position = position;
        Annotations = annotations;
        PhredScore = phredScore;
        IsHomozygous = isHomozygous;
        HgvsProtein = hgvsProtein;
        HgvsCoding = hgvsCoding;
        CodonInfo = codonInfo;

        Ref = ref.replaceAll("\\*", "");
        Alts = alts.replaceAll("\\*", "");

        CodingEffect = codingEffect;
        Effects = effects;

        mFilterType = NONE;
        mMatchType = MATCH_TYPE_NONE;

        mGermlineAltCount = 0;
        mTumorAltCount = 0;
        mGermlineReadDepth = 0;
        mTumorReadDepth = 0;
        mReadDataSet = false;
        mAdjustedVaf = 0;

        mClinvarMatch = false;
        mClinvarSig = "";
        mClinvarSigInfo = "";

        mSomaticVariant = null;
        mVariantContext = null;
        mEnrichedVariant = null;
    }

    public FilterType filterType() { return mFilterType; }

    public void overrideFilterType(final FilterType type) { mFilterType = type; }

    public void setFilterType(final FilterType type)
    {
        if(mFilterType != NONE) // leave if already set
            return;

        mFilterType = type;
    }

    public void setMatchType(final String matchType) { mMatchType = matchType; }
    public String getMatchType() { return mMatchType; }

    public void setClinvarData(String clinvarSig, String clinvarSigInfo)
    {
        mClinvarMatch = true;
        mClinvarSig = clinvarSig;
        mClinvarSigInfo = clinvarSigInfo;
    }

    public boolean getClinvarMatch() { return mClinvarMatch; }
    public String getClinvarSig() { return mClinvarSig; }
    public String getClinvarSigInfo() { return mClinvarSigInfo; }

    public String asString()
    {
        return String.format("id(%s) location(%s:%d) gene(%s)", VariantId, Chromosome, Position, Gene);
    }

    public int compareTo(@NotNull BachelorGermlineVariant other)
    {
        // sort based on Chromosome then Position
        if(other.Chromosome.equals(Chromosome))
        {
            return Position < other.Position ? -1 : 1;
        }
        else
        {
            int chr = chromosomeToInt(Chromosome);
            int otherChr = chromosomeToInt(other.Chromosome);

            if(chr > 0 && otherChr > 0)
                return chr < otherChr ? -1 : 1;
            else if(chr > 0)
                return -1;
            else if(otherChr > 0)
                return 1;
            else
                return Chromosome.compareTo(other.Chromosome);
        }
    }

    private static int chromosomeToInt(final String chr)
    {
        try
        {
            return Integer.parseInt(chr);
        }
        catch(Exception e)
        {
            return 0;
        }
    }

    public int getGermlineAltCount() { return mGermlineAltCount; }
    public int getGermlineReadDepth() { return mGermlineReadDepth; }
    public int getTumorAltCount() { return mTumorAltCount; }
    public int getTumorRefCount() { return mTumorReadDepth - mTumorAltCount; }
    public int getTumorReadDepth() { return mTumorReadDepth; }

    public void setGermlineData(int altCount, int readDepth)
    {
        mGermlineAltCount = altCount;
        mGermlineReadDepth = readDepth;
    }

    public void setTumorData(int altCount, int readDepth)
    {
        mTumorAltCount = altCount;
        mTumorReadDepth = readDepth;
    }

    public void setReadData(int tumorCount, int tumorReadDepth)
    {
        mTumorAltCount = tumorCount;
        mTumorReadDepth = tumorReadDepth;
        mReadDataSet = true;
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

    public boolean isValid(boolean requireEnrichment)
    {
        return mSomaticVariant != null && (!requireEnrichment || mEnrichedVariant != null) && mVariantContext != null;
    }

    public boolean isLowScore()
    {
        return PhredScore < PHRED_SCORE_CUTOFF && mAdjustedVaf < 0;
    }

    public final SomaticVariant getSomaticVariant() { return mSomaticVariant; }
    public final EnrichedSomaticVariant getEnrichedVariant() { return mEnrichedVariant; }

    public void setSomaticVariant(final SomaticVariant var) { mSomaticVariant = var; }
    public void setVariantContext(final VariantContext var) { mVariantContext = var; }
    public void setEnrichedVariant(final EnrichedSomaticVariant var) { mEnrichedVariant = var; }
}
