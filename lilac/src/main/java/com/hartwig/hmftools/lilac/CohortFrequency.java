package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

public class CohortFrequency
{
    private final Map<HlaAllele,Double> mAlleleFrequencies;

    public CohortFrequency(final String freqFile)
    {
        mAlleleFrequencies = Maps.newHashMap();

        if(!freqFile.isEmpty())
            loadData(freqFile);
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

    private void loadData(final String filename)
    {
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
                double frequency = Double.parseDouble(items[1]);
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
