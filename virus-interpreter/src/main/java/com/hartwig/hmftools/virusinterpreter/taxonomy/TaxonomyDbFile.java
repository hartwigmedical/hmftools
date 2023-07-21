package com.hartwig.hmftools.virusinterpreter.taxonomy;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.virusinterpreter.VirusInterpreterApplication.VI_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

public final class TaxonomyDbFile
{
    public static TaxonomyDb loadFromTsv(final String taxonomyDbTsv) throws IOException
    {
        List<String> linesTaxonomyDb = Files.readAllLines(new File(taxonomyDbTsv).toPath());

        Map<Integer, String> taxidToNameMap = Maps.newHashMap();

        for(String line : linesTaxonomyDb)
        {
            String[] parts = line.split(TSV_DELIM);
            if(parts.length == 2)
            {
                int taxid = Integer.parseInt(parts[0].trim());
                String name = parts[1].trim();
                taxidToNameMap.put(taxid, name);
            }
            else
            {
                VI_LOGGER.warn("Suspicious line detected in taxonomy db tsv: {}", line);
            }
        }

        return new TaxonomyDb(taxidToNameMap);
    }
}
