package com.hartwig.hmftools.svtools.pon;

import java.util.List;

import com.google.common.collect.Lists;

public class LocationCounter
{
    public final int Position;

    private short mCount;

    private final List<LocationCounter> mNextLocations;

    public LocationCounter(int position, boolean initLinks)
    {
        Position = position;
        mCount = 0;

        mNextLocations = initLinks ? Lists.newArrayList() : null;
    }

    public int getCount() { return mCount; }
    public void incrementCount()
    {
        if(mCount < Short.MAX_VALUE - 1)
            ++mCount;
    }

    public final List<LocationCounter> getNextLocations() { return mNextLocations; }

    public String toString()
    {
        if(mNextLocations != null)
            return String.format("pos(%d) links(%d)", Position, mNextLocations.size());
        else
            return String.format("pos(%d) count(%d)", Position, mCount);
    }

    public synchronized static LocationCounter getOrAddLocation(final List<LocationCounter> locations, int position, boolean requireLinks)
    {
        if(locations.isEmpty() || position < locations.get(0).Position)
        {
            // add to start
            LocationCounter location = new LocationCounter(position, requireLinks);
            locations.add(0, location);
            return location;
        }

        int locationCount = locations.size();
        if(position > locations.get(locationCount - 1).Position)
        {
            // add to end
            LocationCounter location = new LocationCounter(position, requireLinks);
            locations.add(location);
            return location;
        }

        // binary search to find position to insert
        int currentIndex = locationCount / 2;
        int lowerIndex = 0;
        int upperIndex = locationCount - 1;
        while(true)
        {
            LocationCounter location = locations.get(currentIndex);

            if(location.Position == position)
                return location;

            if(location.Position < position)
            {
                // new position is higher than the current index
                lowerIndex = currentIndex;
            }
            else
            {
                upperIndex = currentIndex;
            }

            // lower = 12, upper = 14, next = 13
            if(lowerIndex >= upperIndex - 1)
            {
                // check whether to insert below or above the current position
                LocationCounter lowerLocation = locations.get(lowerIndex);
                LocationCounter upperLocation = locations.get(upperIndex);

                if(position == lowerLocation.Position)
                    return lowerLocation;
                else if(position == upperLocation.Position)
                    return upperLocation;

                LocationCounter newLocation = new LocationCounter(position, requireLinks);

                if(position > lowerLocation.Position)
                    locations.add(lowerIndex + 1, newLocation);
                else
                    locations.add(lowerIndex, newLocation);

                return newLocation;
            }

            int nextIndex = (lowerIndex + upperIndex) / 2;

            currentIndex = nextIndex;
        }
    }

    public static boolean isValid(final List<LocationCounter> locations)
    {
        // test for incrementing positions and no duplicates
        for(int i = 0; i < locations.size() - 1; ++i)
        {
            LocationCounter current = locations.get(i);
            LocationCounter next = locations.get(i + 1);

            if(current.Position >= next.Position)
                return false;
        }

        return true;
    }
}
