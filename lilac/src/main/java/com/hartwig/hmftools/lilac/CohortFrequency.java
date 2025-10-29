package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.RARE_ALLELES_FREQ_CUTOFF;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import org.jetbrains.annotations.Nullable;

public class CohortFrequency
{
    private final GeneSelector mGenes;
    private final Map<HlaAllele, Double> mAlleleFrequencies;

    public CohortFrequency(@Nullable final GeneSelector genes, final String freqFile)
    {
        mGenes = genes;
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
        if(!allele.Gene.hasFrequencies())
            return RARE_ALLELES_FREQ_CUTOFF;

        Double frequency = mAlleleFrequencies.get(allele);
        return frequency != null ? frequency : 0.0;
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
                HlaAllele allele = HlaAllele.fromString(alleleStr);
                if(mGenes != null && !mGenes.contains(allele.Gene))
                    continue;

                mAlleleFrequencies.put(allele, frequency);
            }

            LL_LOGGER.info("loaded {} allele frequencies from file({})", mAlleleFrequencies.size(), filename);
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to read cohort allele frequency file({}): {}", filename, e.toString());
        }
    }
}
