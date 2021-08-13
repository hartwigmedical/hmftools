package com.hartwig.hmftools.neo.cohort;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.rna.RnaExpressionMatrix.INVALID_EXP;
import static com.hartwig.hmftools.neo.cohort.DataLoader.loadAlleleCoverage;
import static com.hartwig.hmftools.neo.cohort.DataLoader.loadNeoEpitopes;
import static com.hartwig.hmftools.neo.cohort.DataLoader.loadPredictionData;
import static com.hartwig.hmftools.neo.cohort.DataLoader.loadRnaNeoData;

import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.neo.RnaNeoEpitope;
import com.hartwig.hmftools.common.rna.RnaExpressionMatrix;
import com.hartwig.hmftools.neo.bind.BindScorer;

public class NeoSampleBindTask implements Callable
{
    private final String mSampleId;
    private final NeoCohortConfig mConfig;
    private final BindScorer mScorer;
    private final RnaExpressionMatrix mTransExpression;
    private final CohortWriters mWriters;

    public NeoSampleBindTask(
            final String sampleId, final NeoCohortConfig config, final BindScorer scorer,
            final RnaExpressionMatrix transExpression, final CohortWriters writers)
    {
        mSampleId = sampleId;
        mConfig = config;
        mScorer = scorer;
        mTransExpression = transExpression;
        mWriters = writers;
    }

    @Override
    public Long call()
    {
        processSample();
        return (long)1;
    }

    public void processSample()
    {
        Map<Integer,NeoEpitopeData> neoEpitopeDataMap = loadNeoEpitopes(mSampleId, mConfig.NeoDataDir);

        List<AlleleCoverage> alleleCoverages = loadAlleleCoverage(mSampleId, mConfig.LilacDataDir);

        // List<Boolean> geneLostStatus = getGeneStatus(alleleCoverages);

        List<BindingPredictionData> allPredictions = loadPredictionData(mSampleId, mConfig.McfPredictionsDir);

        List<RnaNeoEpitope> rnaNeoDataList = loadRnaNeoData(mSampleId, mConfig.IsofoxDataDir);

        // organise by neoepitope
        Map<Integer,List<BindingPredictionData>> neoPredictions = Maps.newHashMap();

        // Map<String,List<BindingPredictionData>> peptidePredictions = Maps.newHashMap();

        // organise into map by peptide, avoiding repeated peptides
        for(BindingPredictionData predData : allPredictions)
        {
            double score = mScorer.calcScore(predData.Allele, predData.Peptide);
            double rankPerc = mScorer.calcScoreRank(predData.Allele, predData.Peptide, score);
            double likelihood = mScorer.calcLikelihood(predData.Allele, predData.Peptide, rankPerc);
            predData.setScoreData(score, rankPerc, likelihood);

            List<BindingPredictionData> predictions = neoPredictions.get(predData.NeId);

            if(predictions == null)
            {
                predictions = Lists.newArrayList();
                neoPredictions.put(predData.NeId, predictions);
            }

            predictions.add(predData);

            /*
            List<BindingPredictionData> predictions = peptidePredictions.get(predData.Peptide);

            if(predictions == null)
            {
                predictions = Lists.newArrayList();
                peptidePredictions.put(predData.Peptide, predictions);
            }
            else
            {
                if(predictions.stream().anyMatch(x -> x.Allele.equals(predData.Allele)))
                    continue;
            }

            predictions.add(predData);
             */
        }

        for(Map.Entry<Integer,List<BindingPredictionData>> entry : neoPredictions.entrySet())
        {
            NeoEpitopeData neoData = neoEpitopeDataMap.get(entry.getKey());

            if(mTransExpression.hasSampleId(mSampleId))
            {
                for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
                {
                    neoData.TransExpression[fs] = 0;

                    for(String transName : neoData.Transcripts[fs])
                    {
                        double expression = mTransExpression.getExpression(transName, mSampleId);

                        // distinguish non-existent expression vs zero TPM
                        if(expression != INVALID_EXP)
                            neoData.TransExpression[fs] += expression;
                    }
                }
            }

            RnaNeoEpitope rnaNeoData = rnaNeoDataList.stream()
                    .filter(x -> x.Id == neoData.Id && x.VariantInfo.equals(neoData.VariantInfo)).findFirst().orElse(null);

            if(rnaNeoData != null)
            {
                neoData.RnaNovelFragments = rnaNeoData.FragmentCount;
                neoData.RnaBaseDepth[FS_UP] = rnaNeoData.BaseDepth[FS_UP];
                neoData.RnaBaseDepth[FS_DOWN] = rnaNeoData.BaseDepth[FS_DOWN];
            }

            NeoPredictionData neoPredData = new NeoPredictionData(neoData.Id);

            Map<String,NeoPredictionData> allelePredictions = Maps.newHashMap();
            alleleCoverages.forEach(x -> allelePredictions.put(x.Allele, new NeoPredictionData(neoData.Id)));

            for(BindingPredictionData predData : entry.getValue())
            {
                neoPredData.processPredictionData(predData, mConfig.McfSumFactor);

                NeoPredictionData allelePredData = allelePredictions.get(predData.Allele);
                allelePredData.processPredictionData(predData, mConfig.McfSumFactor);

                if(mConfig.WriteTypes.contains(CohortWriteType.ALLELE_PEPTIDE))
                {
                    AlleleCoverage alleleCoverage = alleleCoverages.stream().filter(x -> x.Allele.equals(predData.Allele)).findFirst().orElse(null);

                    if(alleleCoverage != null)
                        mWriters.writePeptideData(mSampleId, neoData, predData, alleleCoverage);
                }
            }

            for(AlleleCoverage alleleCoverage : alleleCoverages)
            {
                NeoPredictionData allelePredData = allelePredictions.get(alleleCoverage.Allele);
                mWriters.writeNeoData(mSampleId, neoData, allelePredData, alleleCoverage);
            }
        }

        /*
        boolean allValid = true;
        for(Map.Entry<String,List<BindingPredictionData>> entry : peptidePredictions.entrySet())
        {
            String peptide = entry.getKey();
            List<BindingPredictionData> predictions = entry.getValue();

            expandHomozygous(predictions);

            if(predictions.size() != EXPECTED_ALLELE_COUNT)
            {
                NE_LOGGER.error("peptide({}) has incorrect allele prediction count({})", peptide, predictions.size());
                continue;
            }

            // process the 6 alleles using the coverage
            double maxAffinity = predictions.stream().mapToDouble(x -> x.affinity()).max().orElse(0);
            double minPresentation = predictions.stream().mapToDouble(x -> x.presentation()).min().orElse(0);

            for(int alleleIndex = 0; alleleIndex < alleleCoverages.size(); ++alleleIndex)
            {
                AlleleCoverage alleleCoverage = alleleCoverages.get(alleleIndex);

                BindingPredictionData predData = predictions.stream()
                        .filter(x -> x.Allele.equals(alleleCoverage.Allele)).findFirst().orElse(null);

                if(predData == null)
                {
                    NE_LOGGER.error("peptide({}) missing allele({}) prediction", peptide, alleleCoverage.Allele);
                    allValid = false;
                    break;
                }

                NeoEpitopeData neoData = neoEpitopeDataMap.get(predData.NeId);

                if(neoData == null)
                {
                    NE_LOGGER.error("neId({}) neoepitope data not found", predData.NeId);
                    continue;
                }

                mWriters.writePeptideData(mSampleId, neoData, predData, alleleCoverage);
            }
        }

        if(!allValid)
        {
            NE_LOGGER.warn("sampleId({}) has inconsistent allele coverage vs prediction alleles", mSampleId);
            return;
        }
        */

    }

}
