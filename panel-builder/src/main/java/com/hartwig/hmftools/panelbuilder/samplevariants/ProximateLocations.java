package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.Math.abs;
import static java.util.Collections.emptyList;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ProximateLocations
{
    // chromosome -> list of positions
    private final Map<String, List<Integer>> mLocations;
    private final int mNearDistance;

    public ProximateLocations(int nearDistance)
    {
        mLocations = new HashMap<>();
        mNearDistance = nearDistance;
    }

    public record Location(
            String chromosome,
            int position
    )
    {
    }

    public boolean isNearLocation(final Location location)
    {
        return mLocations.getOrDefault(location.chromosome(), emptyList())
                .stream()
                .anyMatch(position -> abs(position - location.position()) <= mNearDistance);
    }

    public void addLocation(final Location location)
    {
        List<Integer> positions = mLocations.computeIfAbsent(location.chromosome(), k -> new ArrayList<>());
        positions.add(location.position());
    }

    public void addLocations(final List<Location> locations)
    {
        locations.forEach(this::addLocation);
    }
}
