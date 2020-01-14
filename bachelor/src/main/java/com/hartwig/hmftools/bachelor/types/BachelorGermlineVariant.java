package com.hartwig.hmftools.bachelor.types;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

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
    public final String MatchType;
    public final String CodonInfo;

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
    public static final String MATCH_TYPE_HOTSPOT = "HotSpot";
    public static final String MATCH_TYPE_WHITELIST = "WhiteList";

    public static final int PHRED_SCORE_CUTOFF = 150;

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

    public String asCsv(boolean shortForm)
    {
        String output = "";

        if(shortForm)
        {
            output = String.format("%s,%s,%s,%s,%s",
                    SampleId, Program, VariantId, Gene, TranscriptId);

            output += String.format(",%s,%d,%s,%s,%s,%s,%s",
                    Chromosome, Position, Ref, Alts, CodingEffect, Effects, Annotations);

            output += String.format(",%s,%s,%d,%s,%s,%s,%d,%d,%d,%d,%s",
                    HgvsProtein, IsHomozygous, PhredScore, HgvsCoding, MatchType,
                    mReadDataSet, mGermlineAltCount, mGermlineReadDepth, mTumorAltCount, mTumorReadDepth, CodonInfo);

            output += String.format(",%s,%s,%s",
                    mClinvarMatch, mClinvarSig, mClinvarSigInfo);

        }

        return output;
    }

    // SampleId,Program,Id,Gene,TranscriptId,Chromosome,Position,Ref,Alt
    // CodingEffect,Effect,Annotations,HgvsProtein,IsHomozygous,PhredScore,HgvsCoding
    // MatchType,HasDepthInfo,GermlineAltCount,GermlineReadDepth,TumorAltCount,TumorReadDepth,CodonInfo
    public static final int COL_INDEX_SAMPLE = 0;
    private static final int COL_INDEX_PROGRAM = 1;
    private static final int COL_INDEX_SV_ID = 2;
    private static final int COL_INDEX_GENE = 3;
    private static final int COL_INDEX_TRAN_ID = 4;
    private static final int COL_INDEX_CHR = 5;
    private static final int COL_INDEX_POS = 6;
    private static final int COL_INDEX_REF = 7;
    private static final int COL_INDEX_ALT = 8;

    private static final int COL_INDEX_CODING_EFFECT = 9;
    private static final int COL_INDEX_EFFECTS = 10;
    public static final int COL_INDEX_ANNOTS = 11;
    private static final int COL_INDEX_PROTEIN = 12;
    private static final int COL_INDEX_HZ = 13;
    private static final int COL_INDEX_PHRED = 14;
    private static final int COL_INDEX_CODING = 15;
    private static final int COL_INDEX_MATCH_TYPE = 16;
    private static final int COL_INDEX_HAS_DEPTH = 17;
    private static final int COL_INDEX_GL_ALT_COUNT = 18;
    private static final int COL_INDEX_GL_READ_DEPTH = 19;
    private static final int COL_INDEX_TUMOR_ALT_COUNT = 20;
    private static final int COL_INDEX_TUMOR_READ_DEPTH = 21;
    private static final int COL_INDEX_CODON_INFO = 22;
    private static final int COL_INDEX_CLINVAR_MATCH = 23;
    private static final int COL_INDEX_CLINVAR_SIG = 24;
    public static final int COL_INDEX_CLINVAR_SIG_INFO = 25;

    public static final int BACHELOR_CSV_FIELD_COUNT = COL_INDEX_CLINVAR_SIG_INFO + 1;

    public static BachelorGermlineVariant fromCsv(final String[] items)
    {
        com.hartwig.hmftools.common.variant.CodingEffect codingEffect
                = com.hartwig.hmftools.common.variant.CodingEffect.valueOf(items[COL_INDEX_CODING_EFFECT]);

        BachelorGermlineVariant bachRecord = new BachelorGermlineVariant(
                items[COL_INDEX_SAMPLE],
                items[COL_INDEX_PROGRAM],
                items[COL_INDEX_SV_ID],
                items[COL_INDEX_GENE],
                items[COL_INDEX_TRAN_ID],
                items[COL_INDEX_CHR],
                Long.parseLong(items[COL_INDEX_POS]),
                items[COL_INDEX_REF],
                items[COL_INDEX_ALT],
                codingEffect,
                items[COL_INDEX_EFFECTS],
                items[COL_INDEX_ANNOTS],
                items[COL_INDEX_PROTEIN],
                Boolean.parseBoolean(items[COL_INDEX_HZ]),
                Integer.parseInt(items[COL_INDEX_PHRED]),
                items[COL_INDEX_CODING],
                items[COL_INDEX_MATCH_TYPE],
                items[COL_INDEX_CODON_INFO]);

        boolean hasClinvarInfo = Boolean.parseBoolean(items[COL_INDEX_CLINVAR_MATCH]);

        if(hasClinvarInfo)
        {
            bachRecord.setClinvarData(items[COL_INDEX_CLINVAR_SIG], items[COL_INDEX_CLINVAR_SIG_INFO]);
        }

        boolean hasDepthInfo = Boolean.parseBoolean(items[COL_INDEX_HAS_DEPTH]);

        int glAltCount = Integer.parseInt(items[COL_INDEX_GL_ALT_COUNT]);
        int glReadDepth = Integer.parseInt(items[COL_INDEX_GL_READ_DEPTH]);
        bachRecord.setGermlineData(glAltCount, glReadDepth);

        if (hasDepthInfo)
        {
            int tumorAltCount = Integer.parseInt(items[COL_INDEX_TUMOR_ALT_COUNT]);
            int tumorReadDepth = Integer.parseInt(items[COL_INDEX_TUMOR_READ_DEPTH]);

            bachRecord.setReadData(tumorAltCount, tumorReadDepth);
        }

        return bachRecord;
    }
}
