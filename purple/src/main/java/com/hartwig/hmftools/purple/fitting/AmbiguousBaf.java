package com.hartwig.hmftools.purple.fitting;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.purple.PurpleConstants.AMBIGUOUS_BAF_CONSOLIDATE_RATIO;
import static com.hartwig.hmftools.purple.PurpleConstants.AMBIGUOUS_BAF_EDGE_RATIO;
import static com.hartwig.hmftools.purple.PurpleConstants.AMBIGUOUS_BAF_MIN;
import static com.hartwig.hmftools.purple.PurpleConstants.AMBIGUOUS_BAF_MIN_VALUE;
import static com.hartwig.hmftools.purple.PurpleConstants.AMBIGUOUS_BAF_OBSERVED_RATIO_MAX;
import static com.hartwig.hmftools.purple.PurpleConstants.AMBIGUOUS_BAF_THRESHOLD;
import static com.hartwig.hmftools.purple.PurpleConstants.BAF_PNT_5;
import static com.hartwig.hmftools.purple.PurpleConstants.AMBIGUOUS_BAF_MAX;
import static com.hartwig.hmftools.purple.PurpleConstants.STANDARD_BAF_TARGET_PERC;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.reseg.PeakTroughFinder;
import com.hartwig.hmftools.purple.reseg.RatioBucketSeries;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class AmbiguousBaf
{
    public static double depthBasedBAF(int averageDepth)
    {
        int minDepth = averageDepth / 2;
        int maxDepth = averageDepth * 3 / 2;
        double totalProbability = 0;

        for(int i = minDepth; i < maxDepth; i++)
        {
            final BinomialDistribution distribution = new BinomialDistribution(i, BAF_PNT_5);
            double probability = 1d * distribution.inverseCumulativeProbability(STANDARD_BAF_TARGET_PERC) / i;
            totalProbability += probability;
        }

        double expectedBaf = totalProbability / (maxDepth - minDepth);
        return max(expectedBaf, AMBIGUOUS_BAF_THRESHOLD);
    }

    private static final double BAF_ROUND_FACTOR = 0.005;

    public static double findAmbiguousBaf(final List<ObservedRegion> segments, int averageTumorDepth)
    {
        // calculated for reference purposes
        double depthBasedBaf = depthBasedBAF(averageTumorDepth);

        RatioBucketSeries rawTumorRatios = new RatioBucketSeries(0.5, 0.6, BAF_ROUND_FACTOR);

        for(ObservedRegion segment : segments)
        {
            if(segment.germlineStatus() != GermlineStatus.DIPLOID || segment.bafCount() <= 0)
                continue;

            if(segment.observedTumorRatio() <= 0 || segment.observedTumorRatio() >= AMBIGUOUS_BAF_OBSERVED_RATIO_MAX)
                continue;

            double ratio = round(segment.observedBAF() / BAF_ROUND_FACTOR) * BAF_ROUND_FACTOR;
            rawTumorRatios.addCount(ratio, segment.bafCount());
        }

        List<RatioBucketSeries.Bucket> observedBAFs = rawTumorRatios.buildSmoothedSeries();

        if(observedBAFs.size() < 2)
            return -1;

        List<Integer> peakIndices = PeakTroughFinder.findLocalPeakOrTroughIndices(observedBAFs, true);

        PeakCandidate strongSupportPeak = findClosestStronglySupportedPeak(observedBAFs, peakIndices);

        if(strongSupportPeak != null)
        {
            // take the average of this peak and its upper edge
            for(int index = strongSupportPeak.index + 1; index < observedBAFs.size(); ++index)
            {
                RatioBucketSeries.Bucket observedBAF = observedBAFs.get(index);

                if(observedBAF.value() < AMBIGUOUS_BAF_EDGE_RATIO * strongSupportPeak.observedBAF.value())
                {
                    double ambiguousBaf = (observedBAF.level() + strongSupportPeak.observedBAF.level()) * 0.5;

                    ambiguousBaf = max(min(ambiguousBaf, AMBIGUOUS_BAF_MAX), AMBIGUOUS_BAF_MIN);

                    PPL_LOGGER.debug(format("average tumor depth(%d) ambiguous BAF(%.3f from peak=%.3f upper=%.3f std=%.3f)",
                            averageTumorDepth, ambiguousBaf, strongSupportPeak.observedBAF.level(), observedBAF.level(),
                            depthBasedBaf));

                    // take the higher of the depth-based method and the peak method
                    return max(ambiguousBaf, depthBasedBaf);
                }
            }
        }

        PPL_LOGGER.debug(format("average tumor depth(%d) ambiguous BAF(std=%.3f)", averageTumorDepth, depthBasedBaf));
        return depthBasedBaf;
    }

    private static PeakCandidate findClosestStronglySupportedPeak(
            final List<RatioBucketSeries.Bucket> observedBAFs, final List<Integer> peakIndices)
    {
        // find all raw peaks (i.e. BAF levels where the next and previous BAF has a lower average) with a moving average >= 100
        List<PeakCandidate> peakCandidates = Lists.newArrayList();

        for(int peakIndex : peakIndices)
        {
            RatioBucketSeries.Bucket observedBAF = observedBAFs.get(peakIndex);

            if(observedBAF.value() < AMBIGUOUS_BAF_MIN_VALUE)
                continue;

            // for each peak, find the BAF level before and after it where the average drops below 2/3 of the peak value. If either BAF level
            // doesn’t exist, the peak is discarded. Call these BAF values the ‘left edge’ and ‘right edge' of the peak
            Double lowerEdge = null;
            Double upperEdge = null;

            for(int i = 0; i <= 1; ++i)
            {
                boolean searchUp = (i == 0);

                int index = peakIndex + (searchUp ? 1 : -1);

                while(index >= 0 && index < observedBAFs.size())
                {
                    RatioBucketSeries.Bucket nextBAF = observedBAFs.get(index);

                    if(nextBAF.value() < observedBAF.value() * AMBIGUOUS_BAF_EDGE_RATIO)
                    {
                        if(searchUp)
                            upperEdge = nextBAF.level();
                        else
                            lowerEdge = nextBAF.level();

                        break;
                    }
                    else if(nextBAF.value() >= observedBAF.value())
                    {
                        break;
                    }

                    index += searchUp ? 1 : -1;
                }
            }

            if(lowerEdge != null && upperEdge != null)
                peakCandidates.add(new PeakCandidate(peakIndex, observedBAF, lowerEdge, upperEdge));
        }

        if(peakCandidates.isEmpty())
            return null;

        consolidatePeaks(observedBAFs, peakCandidates);

        PeakCandidate closestBAF = null;
        double closestDiff = 0;

        for(PeakCandidate candidate : peakCandidates)
        {
            double diff = abs(candidate.observedBAF.level() - BAF_PNT_5);
            if(closestBAF == null || diff < closestDiff)
            {
                closestDiff = diff;
                closestBAF = candidate;
            }
        }

        return closestBAF != null ? closestBAF : null;
    }

    private static void consolidatePeaks(final List<RatioBucketSeries.Bucket> observedBAFs, final List<PeakCandidate> peakCandidates)
    {
        // Consolidate adjacent peaks if one peak BAF overlaps the left-edge to right-edge BAF range of the other,
        // find the local trough (i.e. minimum moving average of all BAF bins between the two peak BAF levels). If both peak averages
        // exceed 1.3x the trough average, both peaks are retained. Otherwise, the smaller peak is dropped. Continue this process until
        // all peaks are processed with no consolidation done

        while(peakCandidates.size() > 1)
        {
            boolean foundOverlap = false;

            int index = 0;
            while(index < peakCandidates.size() - 1)
            {
                PeakCandidate candidate = peakCandidates.get(index);
                PeakCandidate nextCandidate = peakCandidates.get(index + 1);

                if(candidate.upperEdge < nextCandidate.lowerEdge)
                {
                    ++index;
                    continue;
                }

                double lowerPeakValue = candidate.observedBAF.value();
                double upperPeakValue = nextCandidate.observedBAF.value();
                double lowestBaf = max(lowerPeakValue, upperPeakValue);

                for(int i = candidate.index + 1; i < nextCandidate.index; ++i)
                {
                    RatioBucketSeries.Bucket observedBAF = observedBAFs.get(i);
                    lowestBaf = min(lowestBaf, observedBAF.value());
                }

                boolean keepBoth = lowerPeakValue > lowestBaf * AMBIGUOUS_BAF_CONSOLIDATE_RATIO
                        && upperPeakValue > lowestBaf * AMBIGUOUS_BAF_CONSOLIDATE_RATIO;

                if(keepBoth)
                {
                    ++index;
                }
                else
                {
                    foundOverlap = true;

                    if(candidate.observedBAF.value() < nextCandidate.observedBAF.value())
                    {
                        peakCandidates.remove(index);
                    }
                    else
                    {
                        peakCandidates.remove(index + 1);
                    }
                }
            }

            if(!foundOverlap)
                return;
        }
    }

    private record PeakCandidate(int index, RatioBucketSeries.Bucket observedBAF, double lowerEdge, double upperEdge) {}
}
