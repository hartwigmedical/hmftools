package com.hartwig.hmftools.amber.purity;

import static java.lang.String.format;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.hartwig.hmftools.amber.PositionEvidence;

import org.apache.commons.lang3.tuple.Pair;

public class PeakSearch
{
    private final List<CandidatePeakEvaluationResult> mPeaks;

    public PeakSearch(final List<PositionEvidence> evidence, final int threads)
    {
        List<Pair<Double, Double>> searchValues = new SearchGrid().searchValuesAndSteps();
        List<CandidatePeakEvaluation> evaluations = new ArrayList<>();
        for(Pair<Double, Double> pair : searchValues)
        {
            CandidatePeak level = new CandidatePeak(pair.getLeft(), pair.getRight());
            evaluations.add(new CandidatePeakEvaluation(level, evidence));
        }
        try
        {
            ExecutorService executor = Executors.newFixedThreadPool(threads);
            executor.invokeAll(evaluations);
            executor.shutdown();
        }
        catch(InterruptedException e)
        {
            AMB_LOGGER.error("peak search interrupted", e);
        }
        List<CandidatePeakEvaluationResult> results = evaluations.stream()
                .map(CandidatePeakEvaluation::result)
                .toList();
        mPeaks = new LocalMaximaFinder<>(results).maxima();
        for(CandidatePeakEvaluationResult peak : mPeaks)
        {
            AMB_LOGGER.debug(format("actual peak at %.3f with score: %.3f ", peak.candidatePeak().vaf(), peak.score()));
        }
    }

    public List<CandidatePeakEvaluationResult> peaks()
    {
        return mPeaks;
    }
}
