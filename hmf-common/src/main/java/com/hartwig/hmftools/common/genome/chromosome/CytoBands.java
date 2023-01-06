package com.hartwig.hmftools.common.genome.chromosome;

import static java.lang.String.format;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.logging.log4j.core.util.IOUtils;

public class CytoBands
{
    public static String resourceAsString(final RefGenomeVersion refGenomeVersion)
    {
        try
        {
            String resourceFilename = resourceFilename(refGenomeVersion);
            return readResource(resourceFilename);
        }
        catch(IOException e)
        {
            return "";
        }
    }

    public static String resourceFilename(final RefGenomeVersion refGenomeVersion)
    {
        return format("/refgenome/cytoBands.%s.tsv", refGenomeVersion.identifier());
    }

    private static String readResource(final String resource) throws IOException
    {
        InputStream in = CytoBands.class.getResourceAsStream(resource);
        BufferedReader reader = new BufferedReader(new InputStreamReader(in));
        return IOUtils.toString(reader);
    }


}
