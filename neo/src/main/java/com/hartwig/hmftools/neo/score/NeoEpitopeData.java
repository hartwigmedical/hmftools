package com.hartwig.hmftools.neo.score;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.rna.RnaExpressionMatrix.INVALID_EXP;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.score.NeoRnaData.NO_TPM_VALUE;
import static com.hartwig.hmftools.neo.score.TpmCalculator.calcEffectiveTpm;
import static com.hartwig.hmftools.neo.score.TpmMediansCache.CANCER_VALUE;
import static com.hartwig.hmftools.neo.score.TpmMediansCache.PAN_CANCER_VALUE;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.neo.NeoEpitopeType;
import com.hartwig.hmftools.common.neo.RnaNeoEpitope;
import com.hartwig.hmftools.common.rna.RnaExpressionMatrix;

public class NeoEpitopeData
{
    public final int Id;
    public final NeoEpitopeType VariantType;
    public final String VariantInfo;
    // public final String GeneId;
    public final String GeneName;
    public final String UpAminoAcids;
    public final String NovelAminoAcids;
    public final String DownAminoAcids;
    public final List<String>[] Transcripts;
    public final int[] NmdBases;
    public final int[] CodingBasesLength;
    public final int FusedIntronLength;
    public final int[] SkippedAcceptorsDonors;
    public final double VariantCopyNumber;
    public final double CopyNumber;
    public final double SubclonalLikelihood;

    public final NeoRnaData RnaData;

    public final List<PeptideScoreData> mPeptides;

    public double mRawEffectiveTpm;
    public double mEffectiveTpm;
    public double mExpectedTpm;
    public double mCopyNumberFactor;

    public double effectiveTpm() { return mEffectiveTpm; }
    public double rawEffectiveTpm() { return mRawEffectiveTpm; }
    public double expectedTpm() { return mExpectedTpm; }
    public double copyNumberFactor() { return mCopyNumberFactor; }

    public NeoEpitopeData(
            final int id, final NeoEpitopeType variantType, final String variantInfo, final String geneId, final String geneName,
            final String upAminoAcids, final String novelAminoAcids, final String downAminoAcids,
            final List<String> transNamesUp, final List<String> transNamesDown, int nmdBasesMin, int nmdBasesMax,
            final double variantCopyNumber, final double copyNumber, final double subclonalLikelihood,
            int codingBasesLengthMin, int codingBasesLengthMax, int fusedIntronLength, int skippedDonors, int skippedAcceptors)
    {
        Id = id;
        VariantType = variantType;
        VariantInfo = variantInfo;
        // GeneId = geneId;
        GeneName = geneName;
        UpAminoAcids = upAminoAcids;
        NovelAminoAcids = novelAminoAcids;
        DownAminoAcids = downAminoAcids;
        Transcripts = new List[FS_PAIR];
        Transcripts[FS_UP] = transNamesUp;
        Transcripts[FS_DOWN] = transNamesDown;
        CopyNumber = copyNumber;
        VariantCopyNumber = variantCopyNumber;
        SubclonalLikelihood = subclonalLikelihood;

        FusedIntronLength = fusedIntronLength;
        SkippedAcceptorsDonors = new int[] { skippedDonors, skippedAcceptors };
        NmdBases = new int[] { nmdBasesMin, nmdBasesMax };
        CodingBasesLength = new int[] { codingBasesLengthMin, codingBasesLengthMax };

        RnaData = new NeoRnaData();

        mPeptides = Lists.newArrayList();
        mRawEffectiveTpm = NO_TPM_VALUE;
        mExpectedTpm = NO_TPM_VALUE;
    }

    public List<PeptideScoreData> peptides() { return mPeptides; }

    public void addPeptides(final List<PeptideScoreData> peptides)
    {
        mPeptides.addAll(peptides);
    }

    public void setExpressionData(
            final SampleData sampleData, final Map<String,Double> sampleTPMs, final RnaExpressionMatrix transExpressionCache,
            final TpmMediansCache cohortTpmMediansCache)
    {
        if(transExpressionCache == null && sampleTPMs.isEmpty())
            return;

        double[] sampleExpression = {0, 0};
        double[] panCancerTpm = {0, 0};
        double[] cancerTpm = {0, 0};

        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            for(String transName : Transcripts[fs])
            {
                double transExpression = INVALID_EXP;

                if(!sampleTPMs.isEmpty())
                {
                    if(sampleTPMs.containsKey(transName))
                        transExpression = sampleTPMs.get(transName);
                    else
                        NE_LOGGER.warn("sample({}) missing transcript({}) TPM", sampleData.TumorId, transName);
                }
                else
                {
                    transExpression = transExpressionCache.getExpression(transName, sampleData.TumorId);
                }

                // distinguish non-existent expression vs zero TPM
                if(transExpression != INVALID_EXP)
                    sampleExpression[fs] += transExpression;

                double[] tpmValues = cohortTpmMediansCache.getTranscriptTpm(transName, sampleData.CancerType);

                if(tpmValues[PAN_CANCER_VALUE] != NO_TPM_VALUE)
                    panCancerTpm[fs] += tpmValues[PAN_CANCER_VALUE];

                if(tpmValues[CANCER_VALUE] != NO_TPM_VALUE)
                    cancerTpm[fs] += tpmValues[CANCER_VALUE];
            }
        }

        RnaData.setExpression(sampleExpression);
        RnaData.setCohortValues(cancerTpm, panCancerTpm);
    }

    public void setFusionRnaSupport(final List<RnaNeoEpitope> rnaNeoDataList)
    {
        if(!VariantType.isFusion())
            return;

        RnaNeoEpitope matchedRnaData = rnaNeoDataList.stream()
                .filter(x -> x.Id == Id && x.VariantInfo.equals(VariantInfo)).findFirst().orElse(null);

        if(matchedRnaData != null)
        {
            RnaData.setCoverage(matchedRnaData.FragmentCount, matchedRnaData.BaseDepth[FS_UP], matchedRnaData.BaseDepth[FS_DOWN]);
        }
    }

    public void setMutationRnaSupport(final List<SomaticVariant> somaticVariants)
    {
        if(!VariantType.isPointMutation())
            return;

        SomaticVariant matchedRnaData = somaticVariants.stream()
                .filter(x -> x.variantInfo().equals(VariantInfo)).findFirst().orElse(null);

        if(matchedRnaData != null)
        {
            RnaData.setCoverage(matchedRnaData.RnaFragments, matchedRnaData.RnaDepth, matchedRnaData.RnaDepth);
        }
    }

    public boolean matchingTranscripts(final NeoEpitopeData other, int stream)
    {
        if(Transcripts[stream].size() != other.Transcripts[stream].size())
            return false;

        return Transcripts[stream].stream().allMatch(x -> other.Transcripts[stream].contains(x));
    }

    public void setCalculatedTpmValues(
            final double tpmNormalisationFactor, final double samplePloidy, final double probLowValue, final double probHighValue)
    {
        mRawEffectiveTpm = tpmNormalisationFactor != NO_TPM_VALUE ? RnaData.fragmentSupport() * tpmNormalisationFactor : NO_TPM_VALUE;

        double scLikelihood = VariantType.isPointMutation() ? SubclonalLikelihood : 1; // note use for fusions

        double variantCn = max(VariantCopyNumber, 0);
        double segmentCn = RnaData.hasExpression() ? CopyNumber : samplePloidy;
        double adjustedCopyNumber = segmentCn > 0 ? variantCn / segmentCn : 1;
        mCopyNumberFactor = 1 - scLikelihood * (1 - min(1.0, adjustedCopyNumber));

        mExpectedTpm = RnaData.getTPM(FS_UP) * mCopyNumberFactor;

        if(tpmNormalisationFactor == NO_TPM_VALUE || !RnaData.hasExpression() || mRawEffectiveTpm == NO_TPM_VALUE)
        {
            mEffectiveTpm = mExpectedTpm;
        }
        else
        {
            mEffectiveTpm = calcEffectiveTpm(tpmNormalisationFactor, mRawEffectiveTpm, mExpectedTpm, probLowValue, probHighValue);
        }
    }

    public String toString()
    {
        return format("%d: %s:%s gene(%s)", Id, VariantType, VariantInfo, GeneName);
    }
}
