package com.hartwig.hmftools.ctdna;

import static java.lang.Math.abs;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public class ProximateLocations
{
    private final Map<String, List<Location>> mRegisteredLocations;

    private final byte NO_ORIENTATION = 0;

    public ProximateLocations()
    {
        mRegisteredLocations = Maps.newHashMap();
    }

    public static final int NEAR_DISTANCE = 50;

    public boolean isNearRegisteredLocation(final String chromosome, final int position)
    {
        return isNearRegisteredLocation(chromosome, position, NO_ORIENTATION);
    }

    public boolean isNearRegisteredLocation(final String chromosome, final int position, byte orientation)
    {
        List<Location> positions = mRegisteredLocations.get(chromosome);

        if(positions == null)
            return false;

        return positions.stream()
                .anyMatch(x -> abs(x.Position - position) <= NEAR_DISTANCE && x.Orientations.contains(orientation));
    }

    public void addRegisteredLocation(final String chromosome, final int position)
    {
        addRegisteredLocation(chromosome, position, NO_ORIENTATION);
    }

    public void addRegisteredLocation(final String chromosome, final int position, byte orientation)
    {
        List<Location> positions = mRegisteredLocations.get(chromosome);

        if(positions == null)
        {
            positions = Lists.newArrayList(new Location(position, orientation));
            mRegisteredLocations.put(chromosome, positions);
            return;
        }

        int index = 0;
        while(index < positions.size())
        {
            Location location = positions.get(index);
            if(location.Position == position)
            {
                location.Orientations.add(orientation);
                return;
            }

            if(location.Position > position)
                ++index;

            break;
        }

        positions.add(index, new Location(position, orientation));
    }

    private class Location
    {
        public final int Position;
        public final Set<Byte> Orientations;

        public Location(int position, byte orientation)
        {
            Position = position;
            Orientations = Sets.newHashSet();
            Orientations.add(orientation);
        }
    }
}
