package com.hartwig.hmftools.common.mappability;

import static java.lang.Float.parseFloat;
import static java.lang.Integer.parseInt;
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
import java.util.OptionalDouble;
import java.util.OptionalInt;
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
    // Keyed by chromosome.
    public final Map<String, List<ProbeQualityWindow>> mWindows;

    public static final String CFG_PROBE_QUALITY_FILE = "probe_quality_profile";
    private static final String DESC_PROBE_QUALITY_FILE = "Genome regions to probe quality";

    // Must match the config used for generating the file
    private static final int BASE_WINDOW_LENGTH = 40;
    private static final int BASE_WINDOW_SPACING = 20;
    private static final String FLD_QUALITY_SCORE = "QualityScore";

    private static final double AGGREGATE_SHARPNESS = 5;

    private static final Logger LOGGER = LogManager.getLogger(ProbeQualityProfile.class);

    public ProbeQualityProfile(final String filePath)
    {
        this(loadProbeQualityWindows(filePath));
    }

    private ProbeQualityProfile(final Map<String, List<ProbeQualityWindow>> windows)
    {
        mWindows = windows;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(CFG_PROBE_QUALITY_FILE, true, DESC_PROBE_QUALITY_FILE);
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
            int qualityScoreField = reader.getColumnIndex(FLD_QUALITY_SCORE);

            String curChromosome = null;
            List<ProbeQualityWindow> curWindows = null;

            for(DelimFileReader.Row row : reader)
            {
                String chromosome = row.getRawValue(chromosomeField);
                int start = parseInt(row.getRawValue(startField));
                double qualityScore = parseFloat(row.getRawValue(qualityScoreField));
                int end = start + BASE_WINDOW_LENGTH;
                ProbeQualityWindow window = new ProbeQualityWindow(start, end, (float) qualityScore);

                // File is typically sorted by chromosome so we can improve performance by only doing the hash table lookup when the
                // chromosome changes.
                if(!chromosome.equals(curChromosome))
                {
                    curChromosome = chromosome;
                    curWindows = result.get(curChromosome);
                    // Not using computeIfAbsent() because that is much slower.
                    if(curWindows == null)
                    {
                        curWindows = new ArrayList<>();
                        result.put(curChromosome, curWindows);
                    }
                }

                curWindows.add(window);
            }
        }

        // Store windows sorted, helps computations later.
        result.forEach((chromosome, windows) -> Collections.sort(windows));

        result.forEach((chromosome, windows) ->
                LOGGER.trace("Loaded chromosome {} with {} windows", chromosome, windows.size()));
        LOGGER.debug("Loading complete, secs({})", secondsSinceNow(startTimeMs));
        return result;
    }

    // Compute the final quality score from the probe quality profile.
    // Returns empty optional if the profile doesn't cover the probe region.
    public OptionalDouble computeQualityScore(final ChrBaseRegion probe)
    {
        if(probe.baseLength() < 1)
        {
            throw new IllegalArgumentException("probe length must be >= 1");
        }

        // If the profile doesn't completely cover the probe then we say we can't assess its quality
        // (since the uncovered region could affect the quality significantly).

        List<ProbeQualityWindow> windows = mWindows.get(probe.chromosome());
        if(windows == null)
        {
            // Probe chromosome not covered at all.
            return OptionalDouble.empty();
        }
        OptionalInt windowsStart = findFirstWindowContaining(windows, probe.start());
        if(windowsStart.isEmpty())
        {
            // Probe start position not covered.
            return OptionalDouble.empty();
        }
        // For some reason using Stream skip() and takeWhile() here is extremely slow compared to List.subList().
        int windowsEnd = scanWhileOverlap(windows, windowsStart.getAsInt(), probe.baseRegion());   // Inclusive
        List<ProbeQualityWindow> overlappingWindows = windows.subList(windowsStart.getAsInt(), windowsEnd + 1);
        if(!overlappingWindows.get(overlappingWindows.size() - 1).containsPosition(probe.end()))
        {
            // Probe middle or end not covered.
            return OptionalDouble.empty();
        }
        double qualityScore = aggregateQualityScore(overlappingWindows.stream(), probe.baseRegion());
        return OptionalDouble.of(qualityScore);
    }

    // Efficiently finds the index of the first window that contains a position.
    private static OptionalInt findFirstWindowContaining(final List<ProbeQualityWindow> windows, int position)
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
        if(index < windows.size() && windows.get(index).containsPosition(position))
        {
            return OptionalInt.of(index);
        }
        else
        {
            return OptionalInt.empty();
        }
    }

    // Compute the final quality score from windows which overlap the probe.
    // `windows` should not be empty.
    private static double aggregateQualityScore(final Stream<ProbeQualityWindow> windows, final BaseRegion probe)
    {
        // Using a soft minimum function to aggregate the scores proved to be a good estimator in experiment.
        double[] sums = new double[2];
        windows.forEach(window ->
        {
            // Idea is to scale the risk approximately linearly with the overlap.
            // This biases the quality score upwards progressively as the overlap decreases.
            int overlap = min(window.end(), probe.end()) - max(window.start(), probe.start());
            double overlapFrac = (double) overlap / BASE_WINDOW_LENGTH;
            double adjustedQuality = 1 / ((1 / window.getQualityScore() - 1) * overlapFrac + 1);
            double weight = exp(-AGGREGATE_SHARPNESS * adjustedQuality);
            sums[0] += adjustedQuality * weight;
            sums[1] += weight;
        });
        if(sums[1] > 0)
        {
            return sums[0] / sums[1];
        }
        else
        {
            // Shouldn't occur if windows is nonempty but don't want to return a NaN or inf.
            return 0;
        }
    }

    // Gets the last index >= `index` in `windows` which overlaps `region`.
    private static int scanWhileOverlap(final List<ProbeQualityWindow> windows, int index, final BaseRegion region)
    {
        while(index + 1 < windows.size() && windows.get(index + 1).overlaps(region))
        {
            ++index;
        }
        return index;
    }
}
