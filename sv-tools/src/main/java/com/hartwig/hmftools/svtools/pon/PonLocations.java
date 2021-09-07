package com.hartwig.hmftools.svtools.pon;

import static com.hartwig.hmftools.svtools.pon.LocationCounter.getOrAddLocation;

import java.util.List;

import com.google.common.collect.Lists;

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

    public synchronized void addPosition(final int posStart, final int posEnd)
    {
        LocationCounter locationCounter = getOrAddLocation(Locations, posStart, true);
        LocationCounter endLocationCounter = getOrAddLocation(locationCounter.getNextLocations(), posEnd, false);
        endLocationCounter.incrementCount();
    }

    public synchronized void addPosition(final int posStart)
    {
        LocationCounter locationCounter = getOrAddLocation(Locations, posStart, false);
        locationCounter.incrementCount();
    }

    public static String[] locationValues(final String key) { return key.split("_"); }

    public static String formLocationId(final String chrStart, final String chrEnd, final byte orientStart, final byte orientEnd)
    {
        return String.format("%s_%s_%d_%d", chrStart, chrEnd, orientStart, orientEnd);
    }

    public static String chrStart(final String[] locValues) { return locValues[0]; }
    public static String chrEnd(final String[] locValues) { return locValues[1]; }

    public static byte orientStart(final String[] locValues) { return Byte.parseByte(locValues[2]); }
    public static byte orientEnd(final String[] locValues) { return Byte.parseByte(locValues[3]); }

    public static String formLocationId(final String chromosome, final byte orient)
    {
        return String.format("%s_%d", chromosome, orient);
    }

    public static String chromosome(final String[] locValues) { return locValues[0]; }
    public static byte orientation(final String[] locValues) { return Byte.parseByte(locValues[1]); }
}
