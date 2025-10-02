package com.hartwig.hmftools.cobalt.calculations;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.chromosome.ContigComparator;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class DiploidNormaliser implements ResultsNormaliser
{
    private final int mRollingMedianMaxDistance;
    private final int mRollingMedianMinCoverage;
    private final RefGenomeVersion mRefGenomeVersion;
    private final Map<Chromosome, DiploidRatioNormaliser> chromosomeToNormaliser = new HashMap<>();
    private List<MedianRatio> mMedianRatios;
    private final Map<Chromosome, Double> mChromosomeToExpectedValue = new HashMap<>();

    public DiploidNormaliser(final int rollingMedianMaxDistance, final int rollingMedianMinCoverage,
            final RefGenomeVersion mRefGenomeVersion)
    {
        this.mRollingMedianMaxDistance = rollingMedianMaxDistance;
        this.mRollingMedianMinCoverage = rollingMedianMinCoverage;
        this.mRefGenomeVersion = mRefGenomeVersion;
    }

    @Override
    public void recordValue(final BamRatio bamRatio)
    {
        DiploidRatioNormaliser normaliserForChromosome = chromosomeToNormaliser.computeIfAbsent(bamRatio.mChromosome,
                k -> new DiploidRatioNormaliser(mRollingMedianMaxDistance, mRollingMedianMinCoverage));
        normaliserForChromosome.recordRatio(bamRatio);
    }

    @Override
    public void dataCollectionFinished()
    {
        chromosomeToNormaliser.values().forEach(DiploidRatioNormaliser::dataCollectionFinished);
        mMedianRatios = chromosomeToNormaliser.entrySet().stream().map(entry ->
                {
                    String chromosome = mRefGenomeVersion.versionedChromosome(entry.getKey());
                    DiploidRatioNormaliser value = entry.getValue();
                    return new MedianRatio(chromosome, value.median(), (int) value.count());
                }
        ).sorted((m, n) -> ContigComparator.INSTANCE.compare(m.Chromosome, n.Chromosome)).toList();
        for(CobaltChromosome cobaltChromosome : new CobaltChromosomes(mMedianRatios).chromosomes())
        {
            double expectedRatio;
            if(cobaltChromosome.humanChromosome().equals(HumanChromosome._Y))
            {
                expectedRatio = 1.0;
            }
            else
            {
                expectedRatio = cobaltChromosome.actualRatio();
            }
            chromosomeToNormaliser.get(cobaltChromosome.humanChromosome()).setmExpectedRatio(expectedRatio);
        }
    }

    @Override
    public void normalise(BamRatio bamRatio)
    {
        if(bamRatio.mChromosome.equals(HumanChromosome._Y))
        {
            bamRatio.setDiploidAdjustedRatio(bamRatio.ratio());
        }
        else
        {
            double normalised = chromosomeToNormaliser.get(bamRatio.mChromosome).normalise(bamRatio);
            bamRatio.setDiploidAdjustedRatio(normalised);
        }
    }

    public List<MedianRatio> medianRatios()
    {
        return mMedianRatios;
    }
}
