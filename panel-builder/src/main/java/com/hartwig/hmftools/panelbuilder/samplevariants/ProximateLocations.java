package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.Math.abs;
import static java.lang.String.format;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class ProximateLocations
{
    private final Map<String, ArrayList<Location>> mRegisteredLocations;

    private static final byte NO_ORIENTATION = 0;
    private static final int NEAR_DISTANCE = 50;

    public ProximateLocations()
    {
        mRegisteredLocations = new HashMap<>();
    }

    public boolean isNearRegisteredLocation(final String chromosome, final int position)
    {
        return isNearRegisteredLocation(chromosome, position, NO_ORIENTATION);
    }

    public boolean isNearRegisteredLocation(final String chromosome, final int position, byte orientation)
    {
        List<Location> positions = mRegisteredLocations.get(chromosome);

        if(positions == null)
        {
            return false;
        }

        return positions.stream()
                .anyMatch(x -> abs(x.Position - position) <= NEAR_DISTANCE && x.Orientations.contains(orientation));
    }

    public void addRegisteredLocation(final String chromosome, final int position)
    {
        addRegisteredLocation(chromosome, position, NO_ORIENTATION);
    }

    public void addRegisteredLocation(final String chromosome, final int position, byte orientation)
    {
        ArrayList<Location> positions = mRegisteredLocations.get(chromosome);

        if(positions == null)
        {
            positions = new ArrayList<>();
            positions.add(new Location(position, orientation));
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
            {
                ++index;
            }

            break;
        }

        positions.add(index, new Location(position, orientation));
    }

    private static class Location
    {
        public final int Position;
        public final Set<Byte> Orientations;

        public Location(int position, byte orientation)
        {
            Position = position;
            Orientations = new HashSet<>();
            Orientations.add(orientation);
        }

        public String toString()
        {
            return format("%d orients(%d)", Position, Orientations.size());
        }
    }
}
