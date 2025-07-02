package com.hartwig.hmftools.common.mappability;

import static java.lang.Math.exp;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.secondsSinceNow;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Reads the data written by GeneUtils ProbeQualityProfiler and provides utilities for calculating probe quality scores based on the data.
// See ProbeQualityProfiler class for more context.
public class ProbeQualityProfile
{
    private final Map<String, List<ProbeQualityWindow>> mWindows;
    // Precomputed region coverage so we can tell if the profile covers a probe region. In typical usage this should cover the whole genome.
    private final Map<String, List<BaseRegion>> mCoverage;

    public static final String PROBE_QUALITY_FILE_CONFIG = "probe_quality_profile";
    public static final String PROBE_QUALITY_FILE_DESC = "Genome regions to probe quality";

    // Must match the config used for generating the file
    private static final int BASE_WINDOW_LENGTH = 40;
    private static final int BASE_WINDOW_SPACING = 20;

    private static final String QUALITY_SCORE_FIELD = "QualityScore";

    private static final Logger LOGGER = LogManager.getLogger(ProbeQualityProfile.class);

    public ProbeQualityProfile(String filePath)
    {
        mWindows = loadProbeQualityWindows(filePath);
        mCoverage = computeWindowCoverage(mWindows);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(PROBE_QUALITY_FILE_CONFIG, false, PROBE_QUALITY_FILE_DESC);
    }

    private static Map<String, List<ProbeQualityWindow>> loadProbeQualityWindows(final String filePath)
    {
        LOGGER.debug("Loading probe quality profile file: {}", filePath);
        long startTimeMs = System.currentTimeMillis();
        Map<String, List<ProbeQualityWindow>> result = new HashMap<>();
        try(DelimFileReader reader = new DelimFileReader(filePath))
        {
            int chromosomeField = reader.getColumnIndex(FLD_CHROMOSOME);
            int startField = reader.getColumnIndex(FLD_POSITION_START);
            int qualityScoreField = reader.getColumnIndex(QUALITY_SCORE_FIELD);
            for(DelimFileReader.Row row : reader)
            {
                String chromosome = row.get(chromosomeField);
                int start = row.getInt(startField);
                double qualityScore = row.getDouble(qualityScoreField);
                int end = start + BASE_WINDOW_LENGTH;
                ProbeQualityWindow window = new ProbeQualityWindow(start, end, (float) qualityScore);
                List<ProbeQualityWindow> windows = result.get(chromosome);
                // Not using computeIfAbsent() because that is much slower.
                if(windows == null)
                {
                    windows = new ArrayList<>();
                    result.put(chromosome, windows);
                }
                windows.add(window);
            }
        }

        // Store windows sorted, helps computations later.
        result.forEach((chromosome, windows) -> Collections.sort(windows));

        result.forEach((chromosome, windows) -> LOGGER.debug("Loaded chromosome {} with {} windows", chromosome, windows.size()));
        LOGGER.debug("Loading complete, secs({})", secondsSinceNow(startTimeMs));
        return result;
    }

    private static Map<String, List<BaseRegion>> computeWindowCoverage(Map<String, List<ProbeQualityWindow>> windows)
    {
        LOGGER.debug("Computing profile window coverage");
        long startTimeMs = System.currentTimeMillis();
        Map<String, List<BaseRegion>> result = windows.entrySet().stream().collect(
                Collectors.toMap(Map.Entry::getKey, entry -> computeWindowCoverage(entry.getValue())));
        LOGGER.debug("Coverage complete, secs({})", secondsSinceNow(startTimeMs));
        return result;
    }

    private static List<BaseRegion> computeWindowCoverage(List<ProbeQualityWindow> windows)
    {
        List<BaseRegion> regions = windows.stream().map(BaseRegion::clone).collect(Collectors.toList());
        regions = BaseRegion.checkMergeOverlapsFast(regions, false);
        return regions;
    }

    // Compute the final quality score from the probe quality profile.
    // Returns empty optional if the profile doesn't cover the probe region.
    public Optional<Double> computeQualityScore(ChrBaseRegion probe)
    {
        // If the profile doesn't completely cover the probe then we say we can't assess its quality
        // (since the uncovered region could affect the quality significantly).
        if(!coversRegion(probe))
        {
            return Optional.empty();
        }

        List<ProbeQualityWindow> windows = mWindows.get(probe.chromosome());
        assert windows != null;
        int windowsStart = findFirstWindowCovering(windows, probe.start());
        Stream<ProbeQualityWindow> overlappingWindows =
                windows.stream().skip(windowsStart).takeWhile(window -> window.start() <= probe.end());
        double qualityScore = aggregateQualityScore(overlappingWindows, probe.baseRegion());
        return Optional.of(qualityScore);
    }

    // Checks if this probe quality profile fully covers the specified region.
    private boolean coversRegion(ChrBaseRegion region)
    {
        List<BaseRegion> coverage = mCoverage.get(region.chromosome());
        return coverage != null && coverage.stream().anyMatch(r -> r.containsRegion(region));
    }

    // Efficiently finds the index of the first window that overlaps a position.
    private static int findFirstWindowCovering(List<ProbeQualityWindow> windows, int position)
    {
        // We are able to find the first window that overlaps using Collections.binarySearch() because, since the windows are of equal size,
        // sorting by start (done previously) implies sorting by end.
        int index = Collections.binarySearch(
                windows,
                new BaseRegion(0, position),   // Dummy interval, only end is used
                Comparator.comparingInt(BaseRegion::end)
        );
        if(index < 0)
        {
            index = -index - 1;
        }
        return index;
    }

    // Compute the final quality score from windows which overlap the probe.
    private static double aggregateQualityScore(Stream<ProbeQualityWindow> windows, BaseRegion probe)
    {
        // Using a soft minimum function to aggregate the scores proved to be a good estimator in experiment.
        final double AGGREGATE_SHARPNESS = 10;
        final double EXP_NORMALISATION = -AGGREGATE_SHARPNESS / BASE_WINDOW_LENGTH;
        double[] sums = new double[2];
        windows.forEach(window ->
        {
            float value = window.getQualityScore();
            int overlap = min(window.end(), probe.end()) - max(window.start(), probe.start());
            double weight = exp(EXP_NORMALISATION * overlap * value);
            sums[0] += value * weight;
            sums[1] += weight;
        });
        double qualityScore = sums[0] / sums[1];
        return qualityScore;
    }

    // TODO: best interface?
    // Covers a target region with nonoverlapping probes such that the probe quality scores are optimised.
    public Optional<List<BaseRegion>> createBestProbesForRegion(ChrBaseRegion target, int probeLength, float minQualityScore)
    {
        if(probeLength < 1)
        {
            throw new IllegalArgumentException("probeLength must >= 1");
        }
        if(target.baseLength() < 1)
        {
            return Optional.of(List.of());
        }
        // If the profile doesn't completely cover the target region then we say we can't assess the quality of possible probes
        // (since the uncovered region could affect the quality significantly).
        if(!coversRegion(target))
        {
            return Optional.empty();
        }

        // Algorithm is to try probe tilings starting from the leftmost possible position. There will be probeLength-1 possible tilings.
        // Pick the one with the highest sum of quality scores where no quality score is lower than minQualityScore.

        int initialFirstProbeStart = target.start() - probeLength + 1;
        List<BaseRegion> bestProbes = null;
        double bestScore = 0;
        Integer windowsStart = null;
        Integer windowsEnd = null;
        for(int firstProbeStart = initialFirstProbeStart; firstProbeStart <= target.start(); ++firstProbeStart)
        {
            List<BaseRegion> probes = tileProbes(firstProbeStart, target.end(), probeLength);
            List<Optional<Double>> qualityScores = new ArrayList<>(probes.size());
            for (int i = 0; i < probes.size(); ++i) {
                // TODO: efficiently calculate quality scores
            }
            boolean qualityOk = qualityScores.stream().allMatch(quality -> quality.isPresent() && quality.get() == minQualityScore);
            if (qualityOk) {
                double score = qualityScores.stream().mapToDouble(Optional::get).sum();
                if (score > bestScore) {
                    bestProbes = probes;
                    bestScore = score;
                }
            }
        }

        if(bestProbes == null)
        {
            return Optional.empty();
        }
        else
        {
            return Optional.of(bestProbes);
        }
    }

    // Create consecutive probes starting exactly at `start` and covering `end`.
    private List<BaseRegion> tileProbes(int start, int end, int probeLength)
    {
        return IntStream.iterate(start, i -> i + probeLength)
                .mapToObj(s -> new BaseRegion(s, s + probeLength))
                .takeWhile(probe -> probe.start() <= end).toList();
    }
}
