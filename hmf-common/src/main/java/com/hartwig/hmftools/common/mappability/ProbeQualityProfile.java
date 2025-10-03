package com.hartwig.hmftools.common.mappability;

import static java.lang.Double.isFinite;
import static java.lang.Float.parseFloat;
import static java.lang.Integer.parseInt;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;
import static java.util.Objects.requireNonNull;

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

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Reads the data written by ProbeQualityProfiler and provides utilities for calculating probe quality scores based on the data.
// See ProbeQualityProfiler class for more context.
public class ProbeQualityProfile
{
    private final int mBaseWindowLength;
    private final int mBaseWindowSpacing;
    private final int mMatchScoreThreshold;
    private final int mMatchScoreOffset;
    // Keyed by chromosome.
    protected final Map<String, List<ProbeQualityWindow>> mWindows;

    public static final String CFG_PROBE_QUALITY_FILE = "probe_quality_profile";
    private static final String DESC_PROBE_QUALITY_FILE = "Genome regions to probe quality";

    // Must match the config used for generating the file.
    private static final int RESOURCE_BASE_WINDOW_LENGTH = 40;
    private static final int RESOURCE_BASE_WINDOW_SPACING = 20;
    private static final int RESOURCE_MATCH_SCORE_THRESHOLD = 24;
    private static final int RESOURCE_MATCH_SCORE_OFFSET = 26;

    private static final String FLD_QUALITY_SCORE = "QualityScore";

    private static final Logger LOGGER = LogManager.getLogger(ProbeQualityProfile.class);

    public ProbeQualityProfile(final Map<String, List<ProbeQualityWindow>> windows, int baseWindowLength, int baseWindowSpacing,
            int matchScoreThreshold, int matchScoreOffset)
    {
        if(baseWindowLength < 1)
        {
            throw new IllegalArgumentException("baseWindowLength must be >= 1");
        }
        if(!(1 <= baseWindowSpacing && baseWindowSpacing <= baseWindowLength))
        {
            throw new IllegalArgumentException("baseWindowSpacing must be in range [1, baseWindowLength]");
        }
        windows.forEach((chromosome, chrWindows) -> chrWindows.forEach(window ->
        {
            // The code will probably still work but better to be safe than sorry.
            // No reason why the windows should not be aligned to the grid we expect.
            if((window.region().start() - 1) % baseWindowSpacing != 0)
            {
                String error = format("Invalid window start: %s:%d", chromosome, window.region().start());
                LOGGER.error(error);
                throw new RuntimeException(error);
            }
            if(window.region().baseLength() != baseWindowLength)
            {
                String error = format("Invalid window length: %s:%s", chromosome, window.region());
                LOGGER.error(error);
                throw new RuntimeException(error);
            }
        }));
        mBaseWindowLength = baseWindowLength;
        mBaseWindowSpacing = baseWindowSpacing;
        mMatchScoreThreshold = matchScoreThreshold;
        mMatchScoreOffset = matchScoreOffset;
        mWindows = windows;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(CFG_PROBE_QUALITY_FILE, true, DESC_PROBE_QUALITY_FILE);
    }

    public static ProbeQualityProfile loadFromResourceFile(final String filePath)
    {
        return new ProbeQualityProfile(
                loadProbeQualityWindows(filePath, RESOURCE_BASE_WINDOW_LENGTH),
                RESOURCE_BASE_WINDOW_LENGTH, RESOURCE_BASE_WINDOW_SPACING,
                RESOURCE_MATCH_SCORE_THRESHOLD, RESOURCE_MATCH_SCORE_OFFSET);
    }

    private static Map<String, List<ProbeQualityWindow>> loadProbeQualityWindows(final String filePath, int baseWindowLength)
    {
        LOGGER.debug("Loading probe quality profile file: {}", filePath);

        long startTimeMs = System.currentTimeMillis();
        Map<String, List<ProbeQualityWindow>> result = new HashMap<>();

        try(DelimFileReader reader = new DelimFileReader(filePath))
        {
            int chromosomeField = requireNonNull(reader.getColumnIndex(FLD_CHROMOSOME));
            int startField = requireNonNull(reader.getColumnIndex(FLD_POSITION_START));
            int qualityScoreField = requireNonNull(reader.getColumnIndex(FLD_QUALITY_SCORE));

            String curChromosome = null;
            List<ProbeQualityWindow> curWindows = null;

            for(DelimFileReader.Row row : reader)
            {
                String chromosome = row.getRawValue(chromosomeField);
                int start = parseInt(row.getRawValue(startField));
                double qualityScore = parseFloat(row.getRawValue(qualityScoreField));
                int end = start + baseWindowLength - 1;
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
        LOGGER.trace("Loading complete, secs({})", secondsSinceNow(startTimeMs));
        return result;
    }

    public int matchScoreThreshold()
    {
        return mMatchScoreThreshold;
    }

    public int matchScoreOffset()
    {
        return mMatchScoreOffset;
    }

    // Compute the final quality score from the probe quality profile.
    // Returns empty optional if the profile doesn't cover the probe region.
    public OptionalDouble computeQualityScore(final ChrBaseRegion probe)
    {
        if(probe.baseLength() < mBaseWindowLength)
        {
            // The aggregation of the quality score is design to work for probes that cover multiple windows.
            throw new IllegalArgumentException(format("probe length must be >= %d", mBaseWindowLength));
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
        if(!overlappingWindows.get(overlappingWindows.size() - 1).region().containsPosition(probe.end()))
        {
            // Probe middle or end not covered.
            return OptionalDouble.empty();
        }
        double qualityScore = aggregateQualityScore(overlappingWindows, probe.baseRegion());
        return OptionalDouble.of(qualityScore);
    }

    // Efficiently finds the index of the first window that contains a position.
    private static OptionalInt findFirstWindowContaining(final List<ProbeQualityWindow> windows, int position)
    {
        // We are able to find the first window that overlaps using Collections.binarySearch() because, since the windows are of equal size,
        // sorting by start (done previously) implies sorting by end.
        int index = Collections.binarySearch(
                windows,
                new ProbeQualityWindow(0, position, 0),   // Dummy interval, only end is used
                Comparator.comparingInt(window -> window.region().end())
        );
        if(index < 0)
        {
            index = -index - 1;
        }
        if(index < windows.size() && windows.get(index).region().containsPosition(position))
        {
            return OptionalInt.of(index);
        }
        else
        {
            return OptionalInt.empty();
        }
    }

    // Compute the final quality score from windows which overlap the probe.
    // `windows` must not be empty.
    private double aggregateQualityScore(final List<ProbeQualityWindow> windows, final BaseRegion probe)
    {
        // Pick the minimum quality score, however need to take into account partial overlap of the windows on the edge of the probe.
        // Windows fully overlapping the probe are used as-is.
        // Windows partially overlapping the probe have their quality score interpolated towards windows which overlap more.

        double centreQuality = windows.get(windows.size() / 2).qualityScore();

        double minQuality = centreQuality;

        // Interpolate towards left side.
        double higherOverlapQuality = centreQuality;
        for(int i = windows.size() / 2 - 1; i >= 0; --i)
        {
            ProbeQualityWindow window = windows.get(i);
            double quality = adjustedWindowQualityScore(window, higherOverlapQuality, probe);
            minQuality = min(minQuality, quality);
            higherOverlapQuality = quality;
        }

        // Interpolate towards right side.
        higherOverlapQuality = centreQuality;
        for(int i = windows.size() / 2 + 1; i < windows.size(); ++i)
        {
            ProbeQualityWindow window = windows.get(i);
            double quality = adjustedWindowQualityScore(window, higherOverlapQuality, probe);
            minQuality = min(minQuality, quality);
            higherOverlapQuality = quality;
        }

        return minQuality;
    }

    private double adjustedWindowQualityScore(final ProbeQualityWindow window, double adjacentQualityScore, final BaseRegion probe)
    {
        int overlap = min(window.region().end(), probe.end()) - max(window.region().start(), probe.start());
        double qualityScore = window.qualityScore();
        if(overlap < mBaseWindowLength)
        {
            double overlapFract = (double) overlap / mBaseWindowLength;
            qualityScore = interpolateQualityScores(window.qualityScore(), adjacentQualityScore, overlapFract);
        }
        return qualityScore;
    }

    private static double interpolateQualityScores(double q1, double q2, double t)
    {
        // Interpolating in the "risk space" because it's the risk that's proportional to the overlap. (The quality is a reciprocal value.)
        double value = q1 * q2 / (q1 + t * q2 - t * q1);
        if(isFinite(value))
        {
            return value;
        }
        else
        {
            // If the input qualities are near 0 the result could be NaN. Desired result is 0.
            return 0;
        }
    }

    // Gets the last index >= `index` in `windows` which overlaps `region`.
    private static int scanWhileOverlap(final List<ProbeQualityWindow> windows, int index, final BaseRegion region)
    {
        while(index + 1 < windows.size() && windows.get(index + 1).region().overlaps(region))
        {
            ++index;
        }
        return index;
    }
}
