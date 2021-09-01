package com.hartwig.hmftools.svtools.pon;

import static com.hartwig.hmftools.svtools.pon.PonBuilder.PON_LOGGER;
import static com.hartwig.hmftools.svtools.pon.PonLocations.formLocationId;

import java.util.Map;

import com.google.common.collect.Maps;

public class PonStore
{
    private final Map<String,PonLocations> mSvLocations;
    private final Map<String,PonLocations> mSglLocations;

    public PonStore()
    {
        mSvLocations = Maps.newHashMap();
        mSglLocations = Maps.newHashMap();
    }

    public Map<String,PonLocations> getSvLocations() { return mSvLocations; }
    public Map<String,PonLocations> getSglLocations() { return mSglLocations; }

    public void addLocation(
            final String chrStart, final String chrEnd, final byte orientStart, final byte orientEnd, final int posStart, final int posEnd)
    {
        PonLocations locations = getOrAddSvLocation(formLocationId(chrStart, chrEnd, orientStart, orientEnd));
        locations.addPosition(posStart, posEnd);
    }

    public void addLocation(
            final String chromosome, final byte orientation, final int position)
    {
        PonLocations locations = getOrAddSglLocation(formLocationId(chromosome, orientation));
        locations.addPosition(position);
    }

    private synchronized PonLocations getOrAddSvLocation(final String locationId)
    {
        PonLocations locations = mSvLocations.get(locationId);

        if(locations != null)
            return locations;

        locations = new PonLocations(locationId);
        mSvLocations.put(locationId, locations);
        return locations;
    }

    private synchronized PonLocations getOrAddSglLocation(final String locationId)
    {
        PonLocations locations = mSglLocations.get(locationId);

        if(locations != null)
            return locations;

        locations = new PonLocations(locationId);
        mSglLocations.put(locationId, locations);
        return locations;
    }

    public synchronized String statsString()
    {
        return String.format("PON stats: SVs(loc=%d variants=%d) SGLs(loc=%d variants=%d)",
                svLocationCount(), svPonCount(), sglLocationCount(), sglPonCount());
    }

    public int svLocationCount() { return mSvLocations.size(); }

    public int svPonCount()
    {
        return mSvLocations.values().stream().mapToInt(x -> x.locationCount()).sum() ;
    }

    public int sglPonCount()
    {
        return mSglLocations.values().stream().mapToInt(x -> x.locationCount()).sum() ;
    }

    public int sglLocationCount() { return mSglLocations.size(); }


    // testing
    public PonLocations getLocation(final String chrStart, final String chrEnd, final byte orientStart, final byte orientEnd)
    {
        return mSvLocations.get(formLocationId(chrStart, chrEnd, orientStart, orientEnd));
    }

    public PonLocations getLocation(final String chromosome, final byte orientation)
    {
        return mSglLocations.get(formLocationId(chromosome, orientation));
    }

    public LocationCounter getLocationCounter(
            final String chrStart, final String chrEnd, final byte orientStart, final byte orientEnd, final int posStart, final int posEnd)
    {
        PonLocations locations = mSvLocations.get(formLocationId(chrStart, chrEnd, orientStart, orientEnd));
        if(locations == null)
            return null;

        LocationCounter startLoc = locations.Locations.stream().filter(x -> x.Position == posStart).findFirst().orElse(null);
        if(startLoc == null)
            return null;

        return startLoc.getNextLocations().stream().filter(x -> x.Position == posEnd).findFirst().orElse(null);
    }

    public LocationCounter getLocationCounter(
            final String chromosome, final byte orientation, final int position)
    {
        PonLocations locations = mSglLocations.get(formLocationId(chromosome, orientation));
        if(locations == null)
            return null;

        return locations.Locations.stream().filter(x -> x.Position == position).findFirst().orElse(null);
    }

}
