package com.hartwig.hmftools.cobalt.e2e;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Preconditions;

public class PanelFileUtilities
{
    public void writePanelNormalisationFile(int start, int stop, File destination) throws IOException
    {
        int interval = stop - start;
        Preconditions.checkArgument(interval > 0);
        Preconditions.checkArgument(interval % 1000 == 0);

        int steps = (interval / 1000) + 1;
        List<String> lines = new ArrayList<>();
        lines.add("chromosome\tposition\trelativeEnrichment");
        int position = start;
        for (int i=0; i<steps; i++)
        {
            int position1Based = position + 1;
            lines.add(String.format("chr1\t%d\t%.4f", position1Based, 1.0001));
            position += 1000;
        }
        Files.write(destination.toPath(), lines);
    }
}
