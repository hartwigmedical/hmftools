package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.Math.abs;
import static java.util.Collections.emptyMap;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class ProximateLocations
{
    // chromosome -> position -> orientations
    private final Map<String, Map<Integer, Set<Byte>>> mLocations;

    private static final byte NO_ORIENTATION = 0;
    private static final int NEAR_DISTANCE = 50;

    public ProximateLocations()
    {
        mLocations = new HashMap<>();
    }

    public record Location(
            String chromosome,
            int position,
            byte orientation
    )
    {
        public Location(String chromosome, int position)
        {
            this(chromosome, position, NO_ORIENTATION);
        }
    }

    // TODO: desired behaviour for NO_ORIENTATION?

    public boolean isNearLocation(final Location location)
    {
        return mLocations.getOrDefault(location.chromosome(), emptyMap())
                .entrySet().stream()
                .anyMatch(entry ->
                        abs(entry.getKey() - location.position()) <= NEAR_DISTANCE && entry.getValue().contains(location.orientation()));
    }

    public void addLocation(final Location location)
    {
        Map<Integer, Set<Byte>> positions = mLocations.computeIfAbsent(location.chromosome(), k -> new HashMap<>());
        Set<Byte> orientations = positions.computeIfAbsent(location.position(), k -> new HashSet<>());
        orientations.add(location.orientation());
    }

    public void addLocations(final List<Location> locations)
    {
        locations.forEach(this::addLocation);
    }
}
