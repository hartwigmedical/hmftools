package com.hartwig.hmftools.common.sequencing;

import org.jetbrains.annotations.Nullable;

public class IlluminaBamUtils
{
    public record IlluminaReadNameAttributes(
            String instrumentName, int runId, String flowcellId, int flowcellLane, int tileNumber, int xCoord, int yCoord)
    {
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
