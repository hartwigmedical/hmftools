package com.hartwig.hmftools.esvee.pon_gen;

import static com.hartwig.hmftools.common.genome.region.Orientation.fromByteStr;
import static com.hartwig.hmftools.esvee.pon_gen.LocationCounter.getOrAddLocation;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.Orientation;

public class PonLocations
{
    public String LocationId;
    public final List<LocationCounter> Locations;

    public PonLocations(final String locationId)
    {
        LocationId = locationId;
        Locations = Lists.newArrayList();
    }

    public int locationCount()
    {
        return Locations.stream().mapToInt(x -> x.getNextLocations() != null ? x.getNextLocations().size() : 1).sum();
    }

    public static String[] locationValues(final String key) { return key.split("_"); }

    // methods for SVs
    public synchronized void addPosition(final int posStart, final int posEnd)
    {
        LocationCounter locationCounter = getOrAddLocation(Locations, posStart, true);
        LocationCounter endLocationCounter = getOrAddLocation(locationCounter.getNextLocations(), posEnd, false);
        endLocationCounter.incrementCount();
    }

    public static String formSvLocationId(final String chrStart, final String chrEnd, final byte orientStart, final byte orientEnd)
    {
        return String.format("%s_%s_%d_%d", chrStart, chrEnd, orientStart, orientEnd);
    }

    public static String chrStart(final String[] locValues) { return locValues[0]; }
    public static String chrEnd(final String[] locValues) { return locValues[1]; }

    public static Orientation orientStart(final String[] locValues) { return fromByteStr(locValues[2]); }
    public static Orientation orientEnd(final String[] locValues) { return fromByteStr(locValues[3]); }

    // methods for SGLs
    public synchronized void addPosition(final int posStart)
    {
        LocationCounter locationCounter = getOrAddLocation(Locations, posStart, false);
        locationCounter.incrementCount();
    }

    public static String formSglLocationId(final String chromosome, final byte orient)
    {
        return String.format("%s_%d", chromosome, orient);
    }

    public static String chromosome(final String[] locValues) { return locValues[0]; }
    public static Orientation orientation(final String[] locValues) { return fromByteStr(locValues[1]); }
}
