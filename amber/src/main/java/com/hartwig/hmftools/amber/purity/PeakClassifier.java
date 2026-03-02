package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.hartwig.hmftools.amber.PositionEvidence;
import com.hartwig.hmftools.amber.contamination.SearchGrid;

public class PeakClassifier
{
    private final List<VafLevelEvaluationResult> ContaminationPeaks = new ArrayList<>();
    private final List<VafLevelEvaluationResult> CopyNumberPeaks = new ArrayList<>();

    public PeakClassifier(List<PositionEvidence> evidence)
    {
        List<Double> searchValues = new SearchGrid().searchValues();
        List<VafLevelEvaluation> evaluations =
                searchValues.stream().map(value -> new VafLevelEvaluation(new VafLevel(value), evidence)).toList();
        ExecutorService executor = Executors.newFixedThreadPool(4); // TODO thread count
        try
        {
            executor.invokeAll(evaluations);
        }
        catch(InterruptedException e)
        {
            AMB_LOGGER.error("Peak search interrupted", e);
        }

        List<VafLevelEvaluationResult> results = evaluations.stream()
                .filter(VafLevelEvaluation::hasScore)
                .map(VafLevelEvaluation::result)
                .toList();
        List<VafLevelEvaluationResult> peaks = new LocalMaximaFinder<>(results).maxima();
        for(VafLevelEvaluationResult peak : peaks)
        {

        }
    }
}
