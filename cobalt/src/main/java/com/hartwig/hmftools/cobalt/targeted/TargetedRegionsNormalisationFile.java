package com.hartwig.hmftools.cobalt.targeted;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

public class TargetedRegionsNormalisationFile
{
    private final String mPath;

    public TargetedRegionsNormalisationFile(final String mPath)
    {
        this.mPath = mPath;
    }

    public ListMultimap<Chromosome, TargetRegionEnrichment> load()
    {
        ListMultimap<Chromosome, TargetRegionEnrichment> result = ArrayListMultimap.create();
        DelimFileReader reader = new DelimFileReader(mPath);
        reader.stream().forEach(row ->
        {
            Chromosome chr = HumanChromosome.fromString(row.get(0));
            int position = Integer.parseInt(row.get(1));
            double enrichment = Double.parseDouble(row.get(2));
            result.put(chr, new TargetRegionEnrichment(chr, position, enrichment));
        });
        return result;
    }
}
