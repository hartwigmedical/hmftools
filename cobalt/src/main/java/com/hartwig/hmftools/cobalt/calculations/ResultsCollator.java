package com.hartwig.hmftools.cobalt.calculations;

import java.util.List;

import com.google.common.base.Preconditions;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class ResultsCollator
{
    private final RefGenomeVersion GenomeVersion;

    public ResultsCollator(RefGenomeVersion genomeVersion)
    {
        this.GenomeVersion = genomeVersion;
    }

    public ListMultimap<Chromosome, CobaltRatio> collateResults(ListMultimap<Chromosome, BamRatio> tumorResults,
            ListMultimap<Chromosome, BamRatio> referenceResults)
    {
        Preconditions.checkArgument(tumorResults.keySet().equals(referenceResults.keySet()));

        final ListMultimap<Chromosome, CobaltRatio> finalResults = ArrayListMultimap.create();
        tumorResults.keySet().forEach((chromosome) -> {
            List<BamRatio> tumorRatiosForChromosome = tumorResults.get(chromosome);
            List<BamRatio> referenceRatiosForChromosome = referenceResults.get(chromosome);
            Preconditions.checkState(tumorRatiosForChromosome.size() == referenceRatiosForChromosome.size());
            for (int i = 0; i < tumorRatiosForChromosome.size(); i++)
            {
                BamRatio tumorRatio = tumorRatiosForChromosome.get(i);
                BamRatio referenceRatio = referenceRatiosForChromosome.get(i);
                finalResults.put(chromosome, merge(referenceRatio, tumorRatio));
            }
        });
        return finalResults;
    }

    private CobaltRatio merge(BamRatio referenceRatio, BamRatio tumorRatio)
    {
        Preconditions.checkArgument(referenceRatio.mChromosome.equals(tumorRatio.mChromosome));
        Preconditions.checkArgument(referenceRatio.Position == tumorRatio.Position);
        return new CobaltRatio(GenomeVersion.versionedChromosome(referenceRatio.mChromosome), referenceRatio.Position,
                referenceRatio.readDepth(), referenceRatio.ratio(), referenceRatio.gcContent(), referenceRatio.getDiploidAdjustedRatio(),
                tumorRatio.readDepth(), tumorRatio.ratio(), tumorRatio.gcContent());
    }
}
