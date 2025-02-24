package com.hartwig.hmftools.common.sequencing;

import static java.lang.Math.sqrt;
import static java.lang.String.format;

import org.jetbrains.annotations.Nullable;

public class IlluminaBamUtils
{
    public record IlluminaReadNameAttributes(
            String instrumentName, int runId, String flowcellId, int flowcellLane, int tileNumber, int xCoord, int yCoord)
    {
        public String laneKey()
        {
            return format("%s:%d:%s:%d", instrumentName, runId, flowcellId, flowcellLane);
        }

        public String tileKey()
        {
            return format("%s:%d", laneKey(), tileNumber);
        }

        public TileCoord tileCoord()
        {
            return new TileCoord(xCoord, yCoord);
        }
    }

    public record TileCoord(int xCoord, int yCoord)
    {
        public double distance(final TileCoord tileCoord)
        {
            int xDiff = xCoord - tileCoord.xCoord;
            int yDiff = yCoord - tileCoord.yCoord;
            return sqrt(xDiff * xDiff + yDiff * yDiff);
        }
    }

    @Nullable
    public static IlluminaReadNameAttributes getReadNameAttributes(final String readName)
    {
        String[] components = readName.split(":");
        if(components.length < 7)
            return null;

        try
        {
            return new IlluminaReadNameAttributes(components[0], Integer.parseInt(components[1]), components[2], Integer.parseInt(components[3]), Integer.parseInt(components[4]), Integer.parseInt(components[5]), Integer.parseInt(components[6]));
        }
        catch(NumberFormatException e)
        {
            return null;
        }
    }
}
