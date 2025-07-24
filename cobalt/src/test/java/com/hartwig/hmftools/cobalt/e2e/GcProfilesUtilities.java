package com.hartwig.hmftools.cobalt.e2e;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Preconditions;

public class GcProfilesUtilities
{
    public void writeDefaultGcProfile(int start, int stop, File gcProfile) throws IOException
    {
        int interval = stop - start;
        Preconditions.checkArgument(interval > 0);
        Preconditions.checkArgument(interval % 1000 == 0);

        int steps = interval / 1000;
        List<String> lines = new ArrayList<>();

        int position = start;
        for (int i=0; i<steps; i++)
        {
            lines.add(String.format("chr1\t%d\t%.2f\t1\t1", position, 0.5));
            position += 1000;
        }
        Files.write(gcProfile.toPath(), lines);
    }
}
