package com.hartwig.hmftools.cobalt.exclusions;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class ExcludedRegionsFile
{
    private final List<ChrBaseRegion> Regions = new ArrayList<>();

    public ExcludedRegionsFile(BufferedReader reader) throws IOException
    {
        reader.readLine(); // Skip header line

        String line;
        while((line = reader.readLine()) != null)
        {
            String[] fields = line.split("\t");
            {
                String chromosome = fields[0];
                int start = Integer.parseInt(fields[1]);
                int end = Integer.parseInt(fields[2]);
                Regions.add(new ChrBaseRegion(chromosome, start, end));
            }
        }
    }

    public List<ChrBaseRegion> regions()
    {
        return Regions;
    }
}