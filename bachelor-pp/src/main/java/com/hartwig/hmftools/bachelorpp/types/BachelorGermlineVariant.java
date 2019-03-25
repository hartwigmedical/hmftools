package com.hartwig.hmftools.bachelorpp.types;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantConsequence;

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
    public final String MatchType;
    public final String CodonInfo;

    public int mGermlineAltCount;
    public int mGermlineReadDepth;
    public int mTumorAltCount;
    public int mTumorReadDepth;
    private boolean mReadDataSet;
    public final List<String> mEffectsList;

    public String mSignificance;
    public String mDiagnosis;

    private double mAdjustedVaf;

    private SomaticVariant mSomaticVariant;
    private VariantContext mVariantContext;
    private EnrichedSomaticVariant mEnrichedVariant;

    public static int PHRED_SCORE_CUTOFF = 150;

    public BachelorGermlineVariant(String sampleId, String program, String varId,
            String gene, String transcriptId, String chromosome, long position,
            String ref, String alts, CodingEffect codingEffect, String effects, String annotations, String hgvsProtein,
            boolean isHomozygous, int phredScore, String hgvsCoding, String matchType, String codonInfo)
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
        MatchType = matchType;
        CodonInfo = codonInfo;

        Ref = ref.replaceAll("\\*", "");
        Alts = alts.replaceAll("\\*", "");

        CodingEffect = codingEffect;
        Effects = effects;

        mEffectsList = Arrays.stream(effects.split("&")).collect(Collectors.toList());

        mGermlineAltCount = 0;
        mTumorAltCount = 0;
        mGermlineReadDepth = 0;
        mTumorReadDepth = 0;
        mReadDataSet = false;
        mAdjustedVaf = 0;

        mSignificance = "";
        mDiagnosis = "";

        mSomaticVariant = null;
        mVariantContext = null;
        mEnrichedVariant = null;
    }

    public int compareTo(final BachelorGermlineVariant other)
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

    public final List<String> effectsList() { return mEffectsList; }
    public int getGermlineAltCount() { return mGermlineAltCount; }
    public int getGermlineRefCount() { return mGermlineReadDepth - mGermlineAltCount; }
    public int getGermlineReadDepth() { return mGermlineReadDepth; }
    public int getTumorAltCount() { return mTumorAltCount; }
    public int getTumorRefCount() { return mTumorReadDepth - mTumorAltCount; }
    public int getTumorReadDepth() { return mTumorReadDepth; }

    public void setTumorData(int altCount, int readDepth)
    {
        mTumorAltCount = altCount;
        mTumorReadDepth = readDepth;
    }

    public void setTumorAltCount(int count) { mTumorAltCount = count; }
    public final String getDiagnosis() { return mDiagnosis; }
    public final String getSignificance() { return mSignificance; }

    public void setDiagnosis(final String text) { mDiagnosis = text; }
    public void setSignificance(final String text) { mSignificance = text; }

    public void setReadData(int glCount, int glReadDepth, int tumorCount, int tumorReadDepth)
    {
        mGermlineAltCount = glCount;
        mGermlineReadDepth = glReadDepth;
        mTumorAltCount = tumorCount;
        mTumorReadDepth = tumorReadDepth;
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
