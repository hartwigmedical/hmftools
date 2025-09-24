package com.hartwig.hmftools.cobalt.calculations;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.ContigComparator;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class DiploidNormaliser implements ResultsNormaliser
{
    private final int mRollingMedianMaxDistance;
    private final int mRollingMedianMinCoverage;
    private final Map<Chromosome, DiploidRatioNormaliser> chromosomeToNormaliser = new HashMap<>();

    public DiploidNormaliser(final int rollingMedianMaxDistance, final int rollingMedianMinCoverage)
    {
        this.mRollingMedianMaxDistance = rollingMedianMaxDistance;
        this.mRollingMedianMinCoverage = rollingMedianMinCoverage;
    }

    @Override
    public void recordValue(final BamRatio bamRatio)
    {
        DiploidRatioNormaliser normaliserForChromosome = chromosomeToNormaliser.computeIfAbsent(bamRatio.mChromosome,
                k -> new DiploidRatioNormaliser(mRollingMedianMaxDistance, mRollingMedianMinCoverage));
        normaliserForChromosome.recordRatio(bamRatio);
    }

    @Override
    public void recordsAllAdded()
    {
        chromosomeToNormaliser.values().forEach(DiploidRatioNormaliser::dataCollectionFinished);
    }

    @Override
    public void applyNormalisation(BamRatio bamRatio)
    {
        if (bamRatio.mChromosome.equals(HumanChromosome._Y))
        {
            bamRatio.setDiploidAdjustedRatio(bamRatio.ratio());
        }
        else
        {
            double normalised = chromosomeToNormaliser.get(bamRatio.mChromosome).normalise(bamRatio);
            bamRatio.setDiploidAdjustedRatio(normalised);
        }
    }

    public List<MedianRatio> medianRatios(RefGenomeVersion refGenome)
    {
        return chromosomeToNormaliser.entrySet().stream().map(entry ->
                {
                    String chromosome = refGenome.versionedChromosome(entry.getKey());
                    DiploidRatioNormaliser value = entry.getValue();
                    return new MedianRatio(chromosome, value.median(), (int) value.count());
                }
        ).sorted((m, n) -> ContigComparator.INSTANCE.compare(m.Chromosome, n.Chromosome)).toList();
    }
}
