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

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.OptionalInt;
import java.util.stream.IntStream;

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
    protected final Map<String, ProbeQualityWindows> mWindows;

    public static final String CFG_PROBE_QUALITY_FILE = "probe_quality_profile";
    private static final String DESC_PROBE_QUALITY_FILE = "Genome regions to probe quality";

    // Must match the config used for generating the file.
    private static final int RESOURCE_BASE_WINDOW_LENGTH = 40;
    private static final int RESOURCE_BASE_WINDOW_SPACING = 20;
    private static final int RESOURCE_MATCH_SCORE_THRESHOLD = 24;
    private static final int RESOURCE_MATCH_SCORE_OFFSET = 26;

    private static final String FLD_QUALITY_SCORE = "QualityScore";

    private static final Logger LOGGER = LogManager.getLogger(ProbeQualityProfile.class);

    private ProbeQualityProfile(final Map<String, ProbeQualityWindows> windows, int baseWindowLength, int baseWindowSpacing,
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
                loadProbeQualityWindows(filePath, RESOURCE_BASE_WINDOW_LENGTH, RESOURCE_BASE_WINDOW_SPACING),
                RESOURCE_BASE_WINDOW_LENGTH, RESOURCE_BASE_WINDOW_SPACING,
                RESOURCE_MATCH_SCORE_THRESHOLD, RESOURCE_MATCH_SCORE_OFFSET);
    }

    protected static class ProbeQualityWindows
    {
        // To save memory, only store the end position and the length, since every window is the same size.
        // Store the end rather than the start because later we need to search on the end.
        private final int mBaseWindowLength;
        // Also store the data in raw arrays to avoid the memory overhead of boxing with ArrayList.
        private int[] mEndPositions;
        private float[] mQualityScores;
        private int mSize;
        private int mCapacity;

        // Overall, the memory optimisations here approximately halve the memory usage.
        // For the full genome profile, the reduction is from 8GB to 4GB.

        public ProbeQualityWindows(int baseWindowLength)
        {
            if(baseWindowLength < 1)
            {
                throw new IllegalArgumentException("baseWindowLength must be >= 1");
            }
            mBaseWindowLength = baseWindowLength;
            mCapacity = 0;
            mSize = 0;
            mEndPositions = new int[mCapacity];
            mQualityScores = new float[mCapacity];
        }

        public int size()
        {
            return mSize;
        }

        public BaseRegion getRegion(int index)
        {
            int endPosition = mEndPositions[index];
            int startPosition = endPosition - mBaseWindowLength + 1;
            return new BaseRegion(startPosition, endPosition);
        }

        public float getQualityScore(int index)
        {
            return mQualityScores[index];
        }

        public void add(int startPosition, float qualityScore)
        {
            int endPosition = startPosition + mBaseWindowLength - 1;
            if(mCapacity <= mSize)
            {
                // Cap growth to prevent huge overallocation.
                int newCapacity = min((mCapacity + 1) * 2, mCapacity + 10_000_000);

                int[] newEndPositions = new int[newCapacity];
                System.arraycopy(mEndPositions, 0, newEndPositions, 0, mCapacity);
                mEndPositions = newEndPositions;

                float[] newQualityScores = new float[newCapacity];
                System.arraycopy(mQualityScores, 0, newQualityScores, 0, mCapacity);
                mQualityScores = newQualityScores;

                mCapacity = newCapacity;
            }
            mEndPositions[mSize] = endPosition;
            mQualityScores[mSize] = qualityScore;
            ++mSize;
        }

        public void sort()
        {
            List<Integer> indices = IntStream.range(0, mSize)
                    .boxed()
                    .sorted((i, j) -> Integer.compare(mEndPositions[i], mEndPositions[j]))
                    .toList();

            int[] newEndPositions = new int[mSize];
            for(int i = 0; i < mSize; ++i)
            {
                newEndPositions[i] = mEndPositions[indices.get(i)];
            }
            mEndPositions = newEndPositions;

            float[] newQualityScores = new float[mSize];
            for(int i = 0; i < mSize; ++i)
            {
                newQualityScores[i] = mQualityScores[indices.get(i)];
            }
            mQualityScores = newQualityScores;
        }
    }

    private static Map<String, ProbeQualityWindows> loadProbeQualityWindows(final String filePath, int baseWindowLength,
            int baseWindowSpacing)
    {
        LOGGER.debug("Loading probe quality profile file: {}", filePath);

        long startTimeMs = System.currentTimeMillis();
        Map<String, ProbeQualityWindows> result = new HashMap<>();

        try(DelimFileReader reader = new DelimFileReader(filePath))
        {
            int chromosomeField = requireNonNull(reader.getColumnIndex(FLD_CHROMOSOME));
            int startField = requireNonNull(reader.getColumnIndex(FLD_POSITION_START));
            int qualityScoreField = requireNonNull(reader.getColumnIndex(FLD_QUALITY_SCORE));

            String curChromosome = null;
            ProbeQualityWindows curWindows = null;

            for(DelimFileReader.Row row : reader)
            {
                String chromosome = row.getRawValue(chromosomeField);
                int start = parseInt(row.getRawValue(startField));
                double qualityScore = parseFloat(row.getRawValue(qualityScoreField));

                if((start - 1) % baseWindowSpacing != 0)
                {
                    String error = format("Invalid window start: %s:%d", chromosome, start);
                    LOGGER.error(error);
                    throw new RuntimeException(error);
                }

                // File is typically sorted by chromosome so we can improve performance by only doing the hash table lookup when the
                // chromosome changes.
                if(!chromosome.equals(curChromosome))
                {
                    curChromosome = chromosome;
                    curWindows = result.get(curChromosome);
                    // Not using computeIfAbsent() because that is much slower.
                    if(curWindows == null)
                    {
                        curWindows = new ProbeQualityWindows(baseWindowLength);
                        result.put(curChromosome, curWindows);
                    }
                }

                curWindows.add(start, (float) qualityScore);
            }
        }

        // Store windows sorted, helps computations later.
        result.forEach((chromosome, windows) -> windows.sort());

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

        ProbeQualityWindows windows = mWindows.get(probe.chromosome());
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
        int windowsEnd = scanWhileOverlap(windows, windowsStart.getAsInt(), probe.baseRegion());   // Inclusive
        if(!windows.getRegion(windowsEnd).containsPosition(probe.end()))
        {
            // Probe middle or end not covered.
            return OptionalDouble.empty();
        }
        double qualityScore = aggregateQualityScore(windows, windowsStart.getAsInt(), windowsEnd + 1, probe.baseRegion());
        return OptionalDouble.of(qualityScore);
    }

    // Efficiently finds the index of the first window that contains a position.
    private static OptionalInt findFirstWindowContaining(final ProbeQualityWindows windows, int position)
    {
        // We are able to find the first window that overlaps using Collections.binarySearch() because, since the windows are of equal size,
        // sorting by start (done previously) implies sorting by end.
        int index = Arrays.binarySearch(windows.mEndPositions, position);
        if(index < 0)
        {
            index = -index - 1;
        }
        if(index < windows.size() && windows.getRegion(index).containsPosition(position))
        {
            return OptionalInt.of(index);
        }
        else
        {
            return OptionalInt.empty();
        }
    }

    // Compute the final quality score from windows which overlap the probe.
    private double aggregateQualityScore(final ProbeQualityWindows windows, int overlapStart, int overlapEnd, final BaseRegion probe)
    {
        // Pick the minimum quality score, however need to take into account partial overlap of the windows on the edge of the probe.
        // Windows fully overlapping the probe are used as-is.
        // Windows partially overlapping the probe have their quality score interpolated towards windows which overlap more.

        int overlapCount = overlapEnd - overlapStart;

        double centreQuality = windows.getQualityScore(overlapCount / 2 + overlapStart);

        double minQuality = centreQuality;

        // Interpolate towards left side.
        double higherOverlapQuality = centreQuality;
        for(int i = overlapCount / 2 - 1; i >= 0; --i)
        {
            double quality = adjustedWindowQualityScore(
                    windows.getRegion(i + overlapStart), windows.getQualityScore(i + overlapStart),
                    higherOverlapQuality, probe);
            minQuality = min(minQuality, quality);
            higherOverlapQuality = quality;
        }

        // Interpolate towards right side.
        higherOverlapQuality = centreQuality;
        for(int i = overlapCount / 2 + 1; i < overlapCount; ++i)
        {
            double quality = adjustedWindowQualityScore(
                    windows.getRegion(i + overlapStart), windows.getQualityScore(i + overlapStart),
                    higherOverlapQuality, probe);
            minQuality = min(minQuality, quality);
            higherOverlapQuality = quality;
        }

        return minQuality;
    }

    private double adjustedWindowQualityScore(final BaseRegion region, double qualityScore, double adjacentQualityScore,
            final BaseRegion probe)
    {
        int overlap = min(region.end(), probe.end()) - max(region.start(), probe.start());
        if(overlap < mBaseWindowLength)
        {
            double overlapFract = (double) overlap / mBaseWindowLength;
            qualityScore = interpolateQualityScores(qualityScore, adjacentQualityScore, overlapFract);
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
    private static int scanWhileOverlap(final ProbeQualityWindows windows, int index, final BaseRegion region)
    {
        while(index + 1 < windows.size() && windows.getRegion(index + 1).overlaps(region))
        {
            ++index;
        }
        return index;
    }
}
