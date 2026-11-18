package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionStartingAt;

import java.util.ArrayList;
import java.util.List;
import java.util.OptionalInt;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

// An ordered list of non-overlapping reference genome regions treated as contiguous in probe-sequence space, 0-indexed [0, length) in list
// order. Lets probe generation continue past a region boundary into the adjacent region (the whole chromosome for DNA, the next exon for
// RNA). Genome-forward and orientation-agnostic; strand handling is the caller's concern.
public class RegionMapping
{
    private final List<ChrBaseRegion> mRegions;
    private final int mLength;

    public RegionMapping(final List<ChrBaseRegion> regions)
    {
        if(regions.isEmpty())
        {
            throw new IllegalArgumentException("Region mapping has no regions");
        }
        for(ChrBaseRegion region : regions)
        {
            if(!region.hasValidPositions())
            {
                throw new IllegalArgumentException("Invalid region: " + region);
            }
        }
        // Regions must be in ascending genome order and not overlap.
        for(int i = 1; i < regions.size(); ++i)
        {
            ChrBaseRegion prev = regions.get(i - 1);
            ChrBaseRegion next = regions.get(i);
            if(prev.overlaps(next))
            {
                throw new IllegalArgumentException("Region mapping regions must not overlap: " + prev + " " + next);
            }
            if(prev.compareTo(next) >= 0)
            {
                throw new IllegalArgumentException("Region mapping regions must be in ascending order: " + prev + " " + next);
            }
        }

        mRegions = List.copyOf(regions);
        mLength = mRegions.stream().mapToInt(ChrBaseRegion::baseLength).sum();
    }

    // Identity mapping - a whole chromosome as one region, so probe-space is genome position - 1 (the DNA case).
    public static RegionMapping wholeChromosome(final String chromosome, int length)
    {
        return new RegionMapping(List.of(new ChrBaseRegion(chromosome, 1, length)));
    }

    public int length()
    {
        return mLength;
    }

    // Whether a probe-space position coincides with a region boundary (start, end, or a junction between adjacent regions). Used to decide
    // when to pin a tiling edge flush to an exon boundary.
    public boolean isRegionBoundary(int probeSpacePosition)
    {
        if(probeSpacePosition == 0)
        {
            return true;
        }
        int cumulative = 0;
        for(ChrBaseRegion region : mRegions)
        {
            cumulative += region.baseLength();
            if(cumulative == probeSpacePosition)
            {
                return true;
            }
            if(cumulative > probeSpacePosition)
            {
                return false;
            }
        }
        return false;
    }

    // Converts a genome position to its probe-space position, or empty if the position is not within any mapped region.
    public OptionalInt toProbeSpacePosition(final String chromosome, int position)
    {
        int regionStart = 0;
        for(ChrBaseRegion region : mRegions)
        {
            if(region.chromosome().equals(chromosome) && position >= region.start() && position <= region.end())
            {
                return OptionalInt.of(regionStart + (position - region.start()));
            }
            regionStart += region.baseLength();
        }
        return OptionalInt.empty();
    }

    // Converts a probe-space position to its genome position. The position must be within [0, length).
    public BasePosition toGenomePosition(int spacePosition)
    {
        if(!(spacePosition >= 0 && spacePosition < mLength))
        {
            throw new IllegalArgumentException("Invalid probe-space position: " + spacePosition);
        }
        int regionStart = 0;
        for(ChrBaseRegion region : mRegions)
        {
            int regionEnd = regionStart + region.baseLength();
            if(spacePosition < regionEnd)
            {
                return new BasePosition(region.chromosome(), region.start() + (spacePosition - regionStart));
            }
            regionStart = regionEnd;
        }
        throw new IllegalStateException("Probe-space position not found");
    }

    // Converts a probe-space range [spaceStart, spaceEnd) to the ordered genome regions it covers, split at region boundaries: one region if
    // within a single mapped region, several if it crosses junctions.
    public List<ChrBaseRegion> toGenomeRegions(int spaceStart, int spaceEnd)
    {
        if(!(spaceStart >= 0 && spaceStart < spaceEnd && spaceEnd <= mLength))
        {
            throw new IllegalArgumentException("Invalid probe-space range");
        }

        List<ChrBaseRegion> result = new ArrayList<>();
        int regionStart = 0;
        for(ChrBaseRegion region : mRegions)
        {
            int regionEnd = regionStart + region.baseLength();
            int intersectionStart = max(regionStart, spaceStart);
            int intersectionEnd = min(regionEnd, spaceEnd);
            if(intersectionStart < intersectionEnd)
            {
                result.add(regionStartingAt(
                        region.chromosome(), region.start() + (intersectionStart - regionStart), intersectionEnd - intersectionStart));
            }
            regionStart = regionEnd;
        }
        return result;
    }
}
