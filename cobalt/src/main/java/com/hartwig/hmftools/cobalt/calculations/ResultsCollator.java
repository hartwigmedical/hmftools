package com.hartwig.hmftools.cobalt.calculations;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkState;

import java.util.List;
import java.util.Set;

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
        checkArgument(!tumorResults.isEmpty() || !referenceResults.isEmpty());
        checkArgument(tumorResults.isEmpty() || referenceResults.isEmpty() || tumorResults.keySet().equals(referenceResults.keySet()));

        final ListMultimap<Chromosome, CobaltRatio> finalResults = ArrayListMultimap.create();
        Set<Chromosome> chromosomes = tumorResults.isEmpty() ? referenceResults.keySet() : tumorResults.keySet();
        chromosomes.forEach((chromosome) ->
        {
            List<BamRatio> tumorRatiosForChromosome = tumorResults.get(chromosome);
            List<BamRatio> referenceRatiosForChromosome = referenceResults.get(chromosome);
            checkState(tumorRatiosForChromosome.isEmpty() ||
                    referenceRatiosForChromosome.isEmpty() ||
                    tumorRatiosForChromosome.size() == referenceRatiosForChromosome.size());
            int size = tumorRatiosForChromosome.isEmpty() ? referenceRatiosForChromosome.size() : tumorRatiosForChromosome.size();
            for(int i = 0; i < size; i++)
            {
                BamRatio tumorRatio = tumorRatiosForChromosome.isEmpty() ? null : tumorRatiosForChromosome.get(i);
                BamRatio referenceRatio = referenceRatiosForChromosome.isEmpty() ? null : referenceRatiosForChromosome.get(i);
                finalResults.put(chromosome, merge(referenceRatio, tumorRatio));
            }
        });
        return finalResults;
    }

    private CobaltRatio merge(BamRatio referenceRatio, BamRatio tumorRatio)
    {
        checkArgument(referenceRatio != null || tumorRatio != null);
        checkArgument(referenceRatio == null || tumorRatio == null || referenceRatio.mChromosome.equals(tumorRatio.mChromosome));
        checkArgument(referenceRatio == null || tumorRatio == null || referenceRatio.Position == tumorRatio.Position);
        if (tumorRatio == null)
        {
            return referenceOnlyCobaltRatio(referenceRatio);
        }
        if (referenceRatio == null)
        {
            return tumorOnlyCobaltRatio(tumorRatio);
        }
        return new CobaltRatio(GenomeVersion.versionedChromosome(referenceRatio.mChromosome), referenceRatio.Position,
                referenceRatio.readDepth(), referenceRatio.ratio(), referenceRatio.gcContent(), referenceRatio.getDiploidAdjustedRatio(),
                tumorRatio.readDepth(), tumorRatio.ratio(), tumorRatio.gcContent());
    }

    private CobaltRatio tumorOnlyCobaltRatio(BamRatio tumorRatio)
    {
        checkArgument(tumorRatio != null);
        return new CobaltRatio(GenomeVersion.versionedChromosome(tumorRatio.mChromosome), tumorRatio.Position,
                -1.0, -1.0, -1.0, -1.0,
                tumorRatio.readDepth(), tumorRatio.ratio(), tumorRatio.gcContent());
    }
    private CobaltRatio referenceOnlyCobaltRatio(BamRatio referenceRatio)
    {
        checkArgument(referenceRatio != null);
        return new CobaltRatio(GenomeVersion.versionedChromosome(referenceRatio.mChromosome), referenceRatio.Position,
                referenceRatio.readDepth(), referenceRatio.ratio(), referenceRatio.gcContent(), referenceRatio.getDiploidAdjustedRatio(),
                -1.0, -1.0, -1.0);
    }
}
