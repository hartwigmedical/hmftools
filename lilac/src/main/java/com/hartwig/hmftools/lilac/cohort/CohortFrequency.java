package com.hartwig.hmftools.lilac.cohort;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class CohortFrequency
{
    private final Map<HlaAllele, Double> mAlleleFrequencies;

    private static final String ALLELE_FREQUENCY_FILE = "allele_freq_file";

    public CohortFrequency(final CommandLine cmd)
    {
        mAlleleFrequencies = Maps.newHashMap();

        if(cmd != null)
            loadData(cmd.getOptionValue(ALLELE_FREQUENCY_FILE));
    }

    public Map<HlaAllele, Double> getAlleleFrequencies()
    {
        return mAlleleFrequencies;
    }

    public double getAlleleFrequency(final HlaAllele allele)
    {
        Double frequency = mAlleleFrequencies.get(allele);
        return frequency != null ? frequency : 0;
    }

    public static void addCmdLineOptions(final Options options)
    {
        options.addOption(ALLELE_FREQUENCY_FILE, true, "File with cohort allele frequencies");
    }

    private void loadData(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            for(final String line : fileData)
            {
                if(line.startsWith("Allele"))
                    continue;

                String[] items = line.split(",");

                if(items.length != 2)
                {
                    LL_LOGGER.error("invalid allele frequency record: {}", line);
                    return;
                }

                String alleleStr = items[0];
                double frequency = Double.parseDouble(items[1]) / 100.0;
                mAlleleFrequencies.put(HlaAllele.fromString(alleleStr), frequency);
            }

            LL_LOGGER.info("loaded {} allele frequencies from file({})", mAlleleFrequencies.size(), filename);
        }
        catch (IOException e)
        {
            LL_LOGGER.error("failed to read cohort allele frequency file({}): {}", filename, e.toString());
        }
    }
}
