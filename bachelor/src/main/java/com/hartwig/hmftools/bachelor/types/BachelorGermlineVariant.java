package com.hartwig.hmftools.bachelor.types;

import static com.hartwig.hmftools.bachelor.types.FilterType.PASS;
import static com.hartwig.hmftools.bachelor.types.PathogenicType.CLINVAR_LIKELY_PATHOGENIC;
import static com.hartwig.hmftools.bachelor.types.PathogenicType.CLINVAR_PATHOGENIC;
import static com.hartwig.hmftools.bachelor.types.PathogenicType.UNANNOTATED;
import static com.hartwig.hmftools.bachelor.types.PathogenicType.WHITE_LIST;

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
    private PathogenicType mPathogenicType;
    private boolean mMatchesRequiredEffect;

    private String mClinvarSig;
    private String mClinvarSigInfo;

    private int mGermlineReadDepth;
    private int mTumorAltCount;
    private int mTumorReadDepth;
    private boolean mReadDataSet;

    private boolean mHasEnrichmentData;
    private double mAdjustedCopyNumber;
    private double mAdjustedVaf;
    private double mMinorAlleleCopyNumber;

    private SomaticVariant mSomaticVariant;
    private VariantContext mVariantContext;

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

        mFilterType = FilterType.NONE;
        mPathogenicType = PathogenicType.NONE;
        mMatchesRequiredEffect = false;

        mTumorAltCount = 0;
        mGermlineReadDepth = 0;
        mTumorReadDepth = 0;
        mReadDataSet = false;

        mHasEnrichmentData = false;
        mAdjustedVaf = 0;
        mAdjustedCopyNumber = 0;
        mMinorAlleleCopyNumber = 0;

        mClinvarSig = "";
        mClinvarSigInfo = "";

        mSomaticVariant = null;
        mVariantContext = null;
    }

    public FilterType filterType() { return mFilterType; }
    public void setFilterType(final FilterType type) { mFilterType = type; }

    public PathogenicType pathogenicType() { return mPathogenicType; }
    public void setPathogenicType(final PathogenicType type) { mPathogenicType = type; }

    public void setMatchRequiredEffect() { mMatchesRequiredEffect = true; }

    public boolean isReportable()
    {
        // PASS and
        // * Pathogenicity in ('WHITE_LIST','CLINVAR_PATHOGENIC','CLINVAR_LIKELY_PATHOGENIC')
        //* Pathogenicity = 'UNANNOTATED' and effect is configured as a known snpeffect
        if(mFilterType != PASS)
            return false;

        if(mPathogenicType == CLINVAR_PATHOGENIC || mPathogenicType == CLINVAR_LIKELY_PATHOGENIC)
            return true;

        if(mMatchesRequiredEffect)
        {
            return (mPathogenicType == UNANNOTATED);
        }
        else
        {
            return (mPathogenicType == WHITE_LIST);
        }
    }

    public void setClinvarData(String clinvarSig, String clinvarSigInfo)
    {
        mClinvarSig = clinvarSig;
        mClinvarSigInfo = clinvarSigInfo;
    }

    public String getClinvarConsolidatedInfo()
    {
        if(mClinvarSigInfo.isEmpty())
            return mClinvarSig;
        else
            return mClinvarSig + ";" + mClinvarSigInfo;
    }

    public String asString()
    {
        return String.format("gene(%s) location(%s:%d)", Gene, Chromosome, Position);
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

    public String refStatus() { return IsHomozygous ? "HOM" : "HET"; }

    public int getGermlineReadDepth() { return mGermlineReadDepth; }
    public int getTumorAltCount() { return mTumorAltCount; }
    public int getTumorRefCount() { return mTumorReadDepth - mTumorAltCount; }
    public int getTumorReadDepth() { return mTumorReadDepth; }

    public void setGermlineData(int readDepth)
    {
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

    public void setEnrichmentData(double adjustedVaf, double maCopyNumber, double adjustedCopyNumber)
    {
        mHasEnrichmentData = true;
        mMinorAlleleCopyNumber = maCopyNumber;
        mAdjustedVaf = adjustedVaf;
        mAdjustedCopyNumber = adjustedCopyNumber;
    }

    public double getAdjustedVaf() { return mAdjustedVaf; }
    public double getMinorAlleleCopyNumber() { return mMinorAlleleCopyNumber; }
    public double getAdjustedCopyNumber() { return mAdjustedCopyNumber; }

    public boolean isBiallelic()
    {
        // Purple annotation hasn't been applied in this case, assume no variant is biallelic.
        if (!mHasEnrichmentData)
        {
            return false;
        }

        double minorAllelePloidy = mAdjustedCopyNumber - (mAdjustedCopyNumber * mAdjustedVaf);
        return (minorAllelePloidy < 0.5);
    }

    public boolean isValid(boolean requireEnrichment)
    {
        if(mSomaticVariant == null)
            return false;

        return !requireEnrichment || mHasEnrichmentData;
    }

    public boolean isLowScore()
    {
        return PhredScore < PHRED_SCORE_CUTOFF && mAdjustedVaf < 0;
    }

    public final SomaticVariant getSomaticVariant() { return mSomaticVariant; }
    public void setSomaticVariant(final SomaticVariant var) { mSomaticVariant = var; }

    public void setVariantContext(final VariantContext var) { mVariantContext = var; }
    public final VariantContext getVariantContext() { return mVariantContext; }
}
