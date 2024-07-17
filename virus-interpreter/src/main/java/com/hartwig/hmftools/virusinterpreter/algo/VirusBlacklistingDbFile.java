package com.hartwig.hmftools.virusinterpreter.algo;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.virusinterpreter.VirusInterpreterApplication.VI_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.jetbrains.annotations.NotNull;

public final class VirusBlacklistingDbFile
{
    public static List<Integer> loadFromTsv(@NotNull String virusBlacklistingTsv) throws IOException
    {
        List<String> linesVirusBlacklistingDb = Files.readAllLines(new File(virusBlacklistingTsv).toPath());

        List<Integer> blacklistedTaxids = new ArrayList<>();

        for(String line : linesVirusBlacklistingDb.subList(1, linesVirusBlacklistingDb.size()))
        {
            String[] parts = line.split(TSV_DELIM);
            if(parts.length == 4)
            {
                blacklistedTaxids.add(Integer.parseInt(parts[0].trim()));
            }
            else
            {
                VI_LOGGER.warn("Suspicious line detected in virus blacklisting db tsv: {}", line);
            }
        }

        return blacklistedTaxids;
    }
}