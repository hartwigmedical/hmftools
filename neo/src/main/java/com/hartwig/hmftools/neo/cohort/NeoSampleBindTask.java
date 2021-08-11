package com.hartwig.hmftools.neo.cohort;

import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.cohort.AlleleCoverage.EXPECTED_ALLELE_COUNT;
import static com.hartwig.hmftools.neo.cohort.AlleleCoverage.getGeneStatus;
import static com.hartwig.hmftools.neo.cohort.DataLoader.loadAlleleCoverage;
import static com.hartwig.hmftools.neo.cohort.DataLoader.loadNeoEpitopes;
import static com.hartwig.hmftools.neo.cohort.DataLoader.loadPredictionData;
import static com.hartwig.hmftools.neo.cohort.BindingPredictionData.expandHomozygous;

import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.neo.bind.BindScorer;

public class NeoSampleBindTask implements Callable
{
    private final String mSampleId;
    private final NeoCohortConfig mConfig;
    private final BindScorer mScorer;
    private final CohortWriters mWriters;

    public NeoSampleBindTask(
            final String sampleId, final NeoCohortConfig config, final BindScorer scorer, final CohortWriters writers)
    {
        mSampleId = sampleId;
        mConfig = config;
        mScorer = scorer;
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

        List<Boolean> geneLostStatus = getGeneStatus(alleleCoverages);

        List<BindingPredictionData> allPredictions = loadPredictionData(mSampleId, mConfig.McfPredictionsDir);

        Map<String,PeptideScores> peptideScores = Maps.newHashMap();

        Map<String,List<BindingPredictionData>> peptidePredictions = Maps.newHashMap();

        // organise into map by peptide, avoiding repeated peptides
        for(BindingPredictionData predData : allPredictions)
        {
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

            double score = mScorer.calcScore(predData.Allele, predData.Peptide);
            double rankPerc = mScorer.calcScoreRank(predData.Allele, predData.Peptide, score);
            double likelihood = mScorer.calcLikelihood(predData.Allele, predData.Peptide, rankPerc);
            predData.setScoreData(score, rankPerc, likelihood);

            predictions.add(predData);
        }

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

            PeptideScores scores = new PeptideScores(peptide, maxAffinity, minPresentation);

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

                /*
                scores.Affinity[NORMAL] = min(scores.Affinity[NORMAL], allelePrediction.Affinity);
                scores.Presentation[NORMAL] = max(scores.Presentation[NORMAL], allelePrediction.Presentation);

                if(!alleleCoverage.isLost())
                {
                    scores.Affinity[TUMOR] = min(scores.Affinity[TUMOR], allelePrediction.Affinity);
                    scores.Presentation[TUMOR] = max(scores.Presentation[TUMOR], allelePrediction.Presentation);
                }

                if(alleleCoverage.isLost() || !geneLostStatus.get(alleleIndex))
                {
                    scores.Affinity[SIM_TUMOR] = min(scores.Affinity[SIM_TUMOR], allelePrediction.Affinity);
                    scores.Presentation[SIM_TUMOR] = max(scores.Presentation[SIM_TUMOR], allelePrediction.Presentation);
                }
                */
            }

            /*
            if(!allValid)
                break;

            peptideScores.put(peptide, scores);
            */
        }

        if(!allValid)
        {
            NE_LOGGER.warn("sampleId({}) has inconsistent allele coverage vs prediction alleles", mSampleId);
            return;
        }

        /*
        SampleSummary sampleSummary = new SampleSummary(peptideScores.size());

        for(Map.Entry<String,PeptideScores> entry : peptideScores.entrySet())
        {
            String peptide = entry.getKey();
            PeptideScores scores = entry.getValue();

            for(int i = 0; i < STATUS_MAX; ++i)
            {
                sampleSummary.Results[i].AffinityTotal += pow(1 / scores.Affinity[i], mConfig.McfSumFactor);

                sampleSummary.Results[i].PresentationTotal += pow(scores.Presentation[i], mConfig.McfSumFactor);

                if(scores.Presentation[i] >= 0.95)
                    ++sampleSummary.Results[i].PresentationCount;
            }

            mWriters.writePeptideData(mSampleId, peptide, scores);
        }

        mWriters.writeSampleSummary(mSampleId, sampleSummary);
        */
    }

}
