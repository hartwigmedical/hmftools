package com.hartwig.hmftools.amber.purity;

import static java.lang.String.format;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.hartwig.hmftools.amber.PositionEvidence;
import com.hartwig.hmftools.amber.contamination.SearchGrid;

import org.apache.commons.lang3.tuple.Pair;

public class PeakSearch
{
    private final List<VafLevelEvaluationResult> Peaks;
    private final List<VafLevel> ContaminationPeaks = new ArrayList<>();
    private final List<VafLevel> CopyNumberPeaks = new ArrayList<>();

    public PeakSearch(List<PositionEvidence> evidence)
    {
        List<Pair<Double, Double>> searchValues = new SearchGrid().searchValuesAndSteps();
        List<VafLevelEvaluation> evaluations = new ArrayList<>();
        for(Pair<Double, Double> pair : searchValues)
        {
            VafLevel level = new VafLevel(pair.getLeft(), pair.getRight());
            evaluations.add(new VafLevelEvaluation(level, evidence));
        }
        ExecutorService executor = Executors.newFixedThreadPool(4); // TODO thread count
        try
        {
            executor.invokeAll(evaluations);
        }
        catch(InterruptedException e)
        {
            AMB_LOGGER.error("Peak search interrupted", e);
        }
        executor.shutdown();
        List<VafLevelEvaluationResult> results = evaluations.stream()
                .filter(VafLevelEvaluation::hasScore)
                .map(VafLevelEvaluation::result)
                .toList();
        for(VafLevelEvaluationResult result : results)
        {
            AMB_LOGGER.debug(format("Potential peak at %.3f with score: %.3f ", result.Vaf().vaf(), result.Score()));
        }
        Peaks = new LocalMaximaFinder<>(results).maxima();
        for(VafLevelEvaluationResult peak : Peaks)
        {
            AMB_LOGGER.debug(format("Actual peak at %.3f with score: %.3f ", peak.Vaf().vaf(), peak.Score()));
        }
    }

    public List<VafLevelEvaluationResult> peaks()
    {
        return Peaks;
    }

}
