package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Miscellaneous utility functionality for base regions.
public class RegionUtils
{
    // Compute regions within `targetRegion` which do not overlap `coveredRegions`.
    public static List<BaseRegion> computeUncoveredRegions(final BaseRegion targetRegion, Stream<BaseRegion> coveredRegions)
    {
        if(!targetRegion.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
        coveredRegions = coveredRegions.peek(coveredRegion ->
        {
            if(!coveredRegion.hasValidPositions())
            {
                throw new IllegalArgumentException("Invalid region");
            }
        });

        // Ignore covered positions which don't overlap the target region, since they can never produce an uncovered region.
        coveredRegions = coveredRegions.filter(targetRegion::overlaps);

        // Sort by start position ascending, then end position ascending.
        coveredRegions = coveredRegions.sorted(Comparator.comparing(BaseRegion::start).thenComparing(BaseRegion::end));

        List<BaseRegion> uncoveredRegions = new ArrayList<>();

        Iterator<BaseRegion> iterator = coveredRegions.iterator();
        // Setting this to just before the target region simplifies handling of first and last uncovered regions.
        int prevCoveredPos = targetRegion.start() - 1;

        // Remaining covered regions.
        while(iterator.hasNext())
        {
            BaseRegion coveredRegion = iterator.next();
            // Possibilities:
            //   - Current region starts at same position as previous:
            //     - And ends >= previous end: nothing to do.
            //   - Current region starts after previous start:
            //     - And ends <= previous end + 1: nothing to do.
            //     - And ends > previous end + 1: uncovered region in between.
            if(coveredRegion.start() > prevCoveredPos + 1)
            {
                int uncoveredStart = prevCoveredPos + 1;
                int uncoveredEnd = min(coveredRegion.start() - 1, targetRegion.end());
                uncoveredRegions.add(new BaseRegion(uncoveredStart, uncoveredEnd));
            }
            prevCoveredPos = max(prevCoveredPos, coveredRegion.end());
        }

        if(prevCoveredPos < targetRegion.end())
        {
            // Covered regions end before the end of the target region, so there is an uncovered region afterward.
            uncoveredRegions.add(new BaseRegion(prevCoveredPos + 1, targetRegion.end()));
        }

        return uncoveredRegions;
    }

    // Checks if `targetRegion` is fully covered by `coveredRegions`.
    public static boolean isCoveredBy(final ChrBaseRegion targetRegion, Stream<ChrBaseRegion> coveredRegions)
    {
        // Similar to computeUncoveredRegions() but without tracking the uncovered regions.
        coveredRegions = coveredRegions.filter(targetRegion::overlaps);
        coveredRegions = coveredRegions.sorted(Comparator.comparing(ChrBaseRegion::start).thenComparing(ChrBaseRegion::end));
        Iterator<ChrBaseRegion> iterator = coveredRegions.iterator();
        int prevCovered = targetRegion.start() - 1;
        if(!iterator.hasNext())
        {
            // No overlapping regions, therefore can't be covered.
            return false;
        }
        while(iterator.hasNext())
        {
            ChrBaseRegion coveredRegion = iterator.next();
            if(coveredRegion.start() > prevCovered + 1)
            {
                // Since the regions are sorted by start, a gap must mean an uncovered region.
                return false;
            }
            prevCovered = max(prevCovered, coveredRegion.end());
        }
        // Could be a gap at the end, if not then the whole region must be covered.
        return prevCovered >= targetRegion.end();
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
        return -(length / 2) + (1 - length % 2);
    }

    public static BaseRegion regionCenteredAt(int centrePosition, int length)
    {
        // Must be inverse of regionCentre()

        if(length < 1)
        {
            throw new IllegalArgumentException("Invalid length");
        }
        int start = centrePosition + regionCentreStartOffset(length);
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
            if(prev != null && region.start() <= prev.end() + 1)
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
}
