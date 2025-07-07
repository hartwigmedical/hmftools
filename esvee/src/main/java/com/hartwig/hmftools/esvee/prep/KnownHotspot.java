package com.hartwig.hmftools.esvee.prep;

import static java.lang.String.format;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.prep.types.JunctionData;
import com.hartwig.hmftools.esvee.prep.types.PrepRead;
import com.hartwig.hmftools.esvee.prep.types.ReadGroup;

public class KnownHotspot
{
    public final ChrBaseRegion RegionStart;
    public final Orientation OrientStart;
    public final ChrBaseRegion RegionEnd;
    public final Orientation OrientEnd;
    public final String GeneInfo;

    public KnownHotspot(
            final ChrBaseRegion regionStart, final Orientation orientStart, final ChrBaseRegion regionEnd, final Orientation orientEnd,
            final String geneInfo)
    {
        RegionStart = regionStart;
        OrientStart = orientStart;
        RegionEnd = regionEnd;
        OrientEnd = orientEnd;
        GeneInfo = geneInfo;
    }

    public boolean matches(
            final String chrStart, final String chrEnd, int posStart, int posEnd, Orientation orientStart, Orientation orientEnd)
    {
        if(RegionStart.containsPosition(chrStart, posStart) && OrientStart == orientStart
        && RegionEnd.containsPosition(chrEnd, posEnd) && OrientEnd == orientEnd)
        {
            return true;
        }

        if(RegionStart.containsPosition(chrEnd, posEnd) && OrientStart == orientEnd
        && RegionEnd.containsPosition(chrStart, posStart) && OrientEnd == orientStart)
        {
            return true;
        }

        return false;
    }

    public boolean sameOrientation() { return OrientStart == OrientEnd;}

    public static boolean readGroupMatchesHotspot(final List<KnownHotspot> knownHotspots, final ReadGroup readGroup)
    {
        for(KnownHotspot knownHotspot : knownHotspots)
        {
            boolean matchesStart = false;
            boolean matchesEnd = false;

            for(PrepRead read : readGroup.reads())
            {
                matchesStart |= knownHotspot.RegionStart.overlaps(read.Chromosome, read.start(), read.end());
                matchesStart |= knownHotspot.RegionStart.containsPosition(read.MateChromosome, read.record().getMateAlignmentStart());
                matchesEnd |= knownHotspot.RegionEnd.overlaps(read.Chromosome, read.start(), read.end());
                matchesEnd |= knownHotspot.RegionEnd.containsPosition(read.MateChromosome, read.record().getMateAlignmentStart());

                if(!matchesStart || !matchesEnd)
                    continue;

                // must match the orientation pairing of the hotspot
                boolean readSameOrientation = read.orientation() == read.mateOrientation();

                if(readSameOrientation == knownHotspot.sameOrientation())
                    return true;
            }
        }

        return false;
    }

    private static boolean readGroupMatchesRegion(final ReadGroup readGroup, final ChrBaseRegion region)
    {
        return readGroup.reads().stream().anyMatch(x -> region.containsPosition(x.MateChromosome, x.record().getMateAlignmentStart()));
    }

    public static boolean junctionMatchesHotspot(final List<KnownHotspot> knownHotspots, final JunctionData junctionData)
    {
        for(KnownHotspot knownHotspot : knownHotspots)
        {
            if(knownHotspot.RegionStart.containsPosition(junctionData.Position))
            {
                // check for reads spanning to the other end
                for(ReadGroup readGroup : junctionData.junctionGroups())
                {
                    if(readGroupMatchesRegion(readGroup, knownHotspot.RegionEnd))
                        return true;
                }
            }
            else if(knownHotspot.RegionEnd.containsPosition(junctionData.Position))
            {
                // check for reads spanning to the other end
                for(ReadGroup readGroup : junctionData.junctionGroups())
                {
                    if(readGroupMatchesRegion(readGroup, knownHotspot.RegionStart))
                        return true;
                }
            }
        }

        return false;
    }

    public String toString()
    {
        return format("genes(%s) start(%s:%d) end(%s:%d)",  GeneInfo, RegionStart, OrientStart.asByte(), RegionEnd, OrientEnd.asByte());
    }
}
