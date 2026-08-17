package com.hartwig.hmftools.purple.reseg;

import static java.lang.String.format;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.purple.region.ObservedRegion;

public final class Resegmentation
{
    public static List<ObservedRegion> run(final List<ObservedRegion> inputRegions, int threads)
    {
        double segmentationPenalty = PenaltyCalculator.calculatePenalty(inputRegions);
        PPL_LOGGER.trace("resegmentation penalty: {}", segmentationPenalty);

        RatioPeakResult peakResult = RatioPeakAnalyser.findRatioPeak(inputRegions);

        List<ObservedRegion> normalisedSegments = GcNormaliser.normaliseRegions(inputRegions, peakResult);

        SupersegmentResult supersegmentResult = SupersegmentIdentifier.identify(normalisedSegments);

        PPL_LOGGER.trace("identified {} supersegments, {} excluded segments",
                supersegmentResult.Supersegments.size(), supersegmentResult.ExcludedSegments.size());

        List<ResegmenterTask> resegmenterTasks = Lists.newArrayList();

        for(int i = 0; i < supersegmentResult.Supersegments.size(); i++)
        {
            resegmenterTasks.add(new ResegmenterTask(supersegmentResult.Supersegments.get(i), segmentationPenalty));
        }

        List<Callable<Void>> threadTasks = resegmenterTasks.stream().collect(Collectors.toList());

        if(!TaskExecutor.executeTasks(threadTasks, threads))
            throw new RuntimeException("resegmentation task execution failed");

        List<ObservedRegion> finalRegions = Lists.newArrayList(supersegmentResult.ExcludedSegments);

        resegmenterTasks.forEach(x -> finalRegions.addAll(x.SupersegmentResults));

        Collections.sort(finalRegions);

        PPL_LOGGER.debug("resegmentation segments({} -> {}) penalty({}) superSegments({} excluded={})",
                inputRegions.size(), finalRegions.size(), format("%.4f", segmentationPenalty),
                        supersegmentResult.Supersegments.size(), supersegmentResult.ExcludedSegments.size());

        return finalRegions;
    }

    private static class ResegmenterTask implements Callable<Void>
    {
        public final List<ObservedRegion> SupersegmentResults;
        private final Supersegment mSupersegment;
        private final double mSegmentationPenalty;

        public ResegmenterTask(final Supersegment supersegment, double segmentationPenalty)
        {
            mSupersegment = supersegment;
            mSegmentationPenalty = segmentationPenalty;
            SupersegmentResults = Lists.newArrayList();
        }

        @Override
        public Void call()
        {
            List<ObservedRegion> resegmentedRegions = SupersegmentResegmenter.resegment(mSupersegment.BothNoneMembers, mSegmentationPenalty);
            SupersegmentResults.addAll(SkippedSegmentReinjector.reinjectSkippedRegions(resegmentedRegions, mSupersegment));
            return null;
        }
    }
}
