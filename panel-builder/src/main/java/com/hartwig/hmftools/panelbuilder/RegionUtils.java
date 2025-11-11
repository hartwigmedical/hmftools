package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Miscellaneous utility functionality for base regions.
public class RegionUtils
{
    // Compute regions within `region1` which do not overlap any of `regions`.
    public static List<BaseRegion> regionNegatedIntersection(final BaseRegion region1, Stream<BaseRegion> regions)
    {
        if(!region1.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
        regions = regions.peek(region ->
        {
            if(!region.hasValidPositions())
            {
                throw new IllegalArgumentException("Invalid region");
            }
        });

        // Ignore positions which don't overlap region1, since they can never be part of the result.
        regions = regions.filter(region1::overlaps);

        // Sort by start position ascending, then end position ascending.
        regions = regions.sorted(Comparator.comparing(BaseRegion::start).thenComparing(BaseRegion::end));

        List<BaseRegion> result = new ArrayList<>();

        Iterator<BaseRegion> iterator = regions.iterator();
        // Setting this to just before region1 simplifies handling of the first and last unoverlapped regions.
        int prevOverlapPos = region1.start() - 1;

        while(iterator.hasNext())
        {
            BaseRegion region = iterator.next();
            // Possibilities:
            //   - Current region starts at same position as previous:
            //     - And ends >= previous end: nothing to do.
            //   - Current region starts after previous start:
            //     - And ends <= previous end + 1: nothing to do.
            //     - And ends > previous end + 1: unoverlapped region in between.
            if(region.start() > prevOverlapPos + 1)
            {
                result.add(new BaseRegion(prevOverlapPos + 1, min(region.start() - 1, region1.end())));
            }
            prevOverlapPos = max(prevOverlapPos, region.end());
        }

        if(prevOverlapPos < region1.end())
        {
            // Overlapped regions end before the end of region1, so there is an unoverlapped region afterward.
            result.add(new BaseRegion(prevOverlapPos + 1, region1.end()));
        }

        return result;
    }

    // Checks if all the bases within `queryRegion` are also within the union of `regions`.
    public static boolean isFullyOverlappedBy(final ChrBaseRegion queryRegion, Stream<ChrBaseRegion> regions)
    {
        if(!queryRegion.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
        regions = regions.peek(region ->
        {
            if(!region.hasValidPositions())
            {
                throw new IllegalArgumentException("Invalid region");
            }
        });

        // Similar to regionNegatedIntersection() but without tracking the resulting regions.
        Iterator<ChrBaseRegion> iterator = regions
                .filter(queryRegion::overlaps)
                .sorted(Comparator.comparing(ChrBaseRegion::start).thenComparing(ChrBaseRegion::end))
                .iterator();
        int prevOverlapped = queryRegion.start() - 1;
        if(!iterator.hasNext())
        {
            // No overlapping regions.
            return false;
        }
        while(iterator.hasNext())
        {
            ChrBaseRegion region = iterator.next();
            if(region.start() > prevOverlapped + 1)
            {
                // Since the regions are sorted by start, a gap must mean an unoverlapped region.
                return false;
            }
            prevOverlapped = max(prevOverlapped, region.end());
        }
        // Could be a gap at the end, if not then the whole region must be overlapped.
        return prevOverlapped >= queryRegion.end();
    }

    public static double regionCentreFloat(final BaseRegion region)
    {
        if(!region.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
        return (region.start() + region.end()) / 2.0;
    }

    public static int regionCentre(final BaseRegion region)
    {
        return (int) regionCentreFloat(region);
    }

    public static BasePosition regionCentre(final ChrBaseRegion region)
    {
        // Must be inverse of regionCenteredAt()

        return new BasePosition(region.chromosome(), regionCentre(region.baseRegion()));
    }

    public static BaseRegion regionStartingAt(int startPosition, int length)
    {
        if(length < 1)
        {
            throw new IllegalArgumentException("Invalid length");
        }
        BaseRegion region = new BaseRegion(startPosition, startPosition + length - 1);
        return region;
    }

    public static ChrBaseRegion regionStartingAt(final String chromosome, int startPosition, int length)
    {
        return ChrBaseRegion.from(chromosome, regionStartingAt(startPosition, length));
    }

    public static int regionCentreStartOffset(int length)
    {
        return (length / 2) - (1 - length % 2);
    }

    public static BaseRegion regionCenteredAt(int centrePosition, int length)
    {
        // Must be inverse of regionCentre()

        if(length < 1)
        {
            throw new IllegalArgumentException("Invalid length");
        }
        int start = centrePosition - regionCentreStartOffset(length);
        int end = start + length - 1;
        BaseRegion region = new BaseRegion(start, end);
        return region;
    }

    public static ChrBaseRegion regionCenteredAt(final BasePosition position, int length)
    {
        return ChrBaseRegion.from(position.Chromosome, regionCenteredAt(position.Position, length));
    }

    public static BaseRegion regionEndingAt(int endPosition, int length)
    {
        if(length < 1)
        {
            throw new IllegalArgumentException("Invalid length");
        }
        BaseRegion region = new BaseRegion(endPosition - length + 1, endPosition);
        return region;
    }

    public static ChrBaseRegion regionEndingAt(final String chromosome, int endPosition, int length)
    {
        return ChrBaseRegion.from(chromosome, regionEndingAt(endPosition, length));
    }

    public static Optional<BaseRegion> regionIntersection(final BaseRegion region1, final BaseRegion region2)
    {
        if(!region1.hasValidPositions() || !region2.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
        BaseRegion result = new BaseRegion(max(region1.start(), region2.start()), min(region1.end(), region2.end()));
        if(result.hasValidPositions())
        {
            return Optional.of(result);
        }
        else
        {
            return Optional.empty();
        }
    }

    public static Optional<ChrBaseRegion> regionIntersection(final ChrBaseRegion region1, final ChrBaseRegion region2)
    {
        if(region1.chromosome().equals(region2.chromosome()))
        {
            return regionIntersection(region1.baseRegion(), region2.baseRegion())
                    .map(intersection -> ChrBaseRegion.from(region1.chromosome(), intersection));
        }
        else
        {
            return Optional.empty();
        }
    }

    public static boolean regionOverlapsOrAdjacent(final BaseRegion region1, final BaseRegion region2)
    {
        return region1.overlaps(region2) || region1.end() + 1 == region2.start() || region2.end() + 1 == region1.start();
    }

    public static List<ChrBaseRegion> mergeOverlapAndAdjacentRegions(Stream<ChrBaseRegion> regions)
    {
        regions = regions.sorted();
        List<ChrBaseRegion> result = new ArrayList<>();
        regions.forEach(region ->
        {
            ChrBaseRegion prev = result.isEmpty() ? null : result.get(result.size() - 1);
            if(prev != null && region.chromosome().equals(prev.chromosome()) && region.start() <= prev.end() + 1)
            {
                if(region.end() > prev.end())
                {
                    // Don't mutate in place because we borrowed the object from the input stream.
                    result.set(result.size() - 1, new ChrBaseRegion(prev.chromosome(), prev.start(), region.end()));
                }
            }
            else
            {
                result.add(region);
            }
        });
        return result;
    }

    public static boolean isRegionValid(final ChrBaseRegion region, final Map<String, Integer> chromosomeLengths)
    {
        Integer chromosomeLength = chromosomeLengths.get(region.chromosome());
        return region.hasValidPositions() && chromosomeLength != null && region.end() <= chromosomeLength;
    }

    public static boolean isPositionValid(final BasePosition position, final Map<String, Integer> chromosomeLengths)
    {
        Integer chromosomeLength = chromosomeLengths.get(position.Chromosome);
        return position.Position >= 1 && chromosomeLength != null && position.Position <= chromosomeLength;
    }

    // Gets subregion indexed by [start, end) taking into account orientation.
    public static ChrBaseRegion getSubregion(final ChrBaseRegion region, final Orientation orientation, int start, int end)
    {
        if(!region.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
        if(!(start >= 0 && end <= region.baseLength()))
        {
            throw new IllegalArgumentException("Invalid start and end");
        }

        int length = end - start;
        ChrBaseRegion result = orientation.isForward()
                ? regionStartingAt(region.chromosome(), region.start() + start, length)
                : regionEndingAt(region.chromosome(), region.end() - start, length);
        return result;
    }
}
