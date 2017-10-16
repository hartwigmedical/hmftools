package com.hartwig.hmftools.common.cobalt;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CobaltPositionFactory {

    private static final Logger LOGGER = LogManager.getLogger(CobaltPositionFactory.class);
    private final Multimap<Chromosome, ModifiableCobaltPosition> result = ArrayListMultimap.create();


    public void addReadCounts(@NotNull final Multimap<String, ReadCount> referenceCount, @NotNull final Multimap<String, ReadCount> tumorCount) {
        final GenomePositionSelector<ReadCount> tumorCountSelector = GenomePositionSelectorFactory.create(tumorCount);

        for (String chromosomeName : referenceCount.keySet()) {
            if (HumanChromosome.contains(chromosomeName)) {
                Chromosome chromosome = HumanChromosome.fromString(chromosomeName);
                for (ReadCount referenceReadCount : referenceCount.get(chromosomeName)) {
                    ModifiableCobaltPosition position = create(referenceReadCount);
                    tumorCountSelector.select(position).ifPresent(x -> position.setTumorReadCount(x.readCount()));
                    result.put(chromosome, position);
                }
            } else {
                LOGGER.info("Excluding unsupported {} chromosome from read ratios", chromosomeName);

            }
        }
    }

    public void addRatios(@NotNull final Multimap<String, ReadRatio> referenceGCRatio,
            @NotNull final Multimap<String, ReadRatio> tumorGCRatio, @NotNull final Multimap<String, ReadRatio> referenceGCDiploidRatio) {

        final GenomePositionSelector<ReadRatio> referenceGCRatioSelector = GenomePositionSelectorFactory.create(referenceGCRatio);
        final GenomePositionSelector<ReadRatio> referenceGCDiploidRatioSelector =
                GenomePositionSelectorFactory.create(referenceGCDiploidRatio);
        final GenomePositionSelector<ReadRatio> tumorGCRatioSelector = GenomePositionSelectorFactory.create(tumorGCRatio);
        for (Chromosome chromosome : result.keySet()) {
            for (ModifiableCobaltPosition cobaltPosition : result.get(chromosome)) {
                referenceGCRatioSelector.select(cobaltPosition).ifPresent(x -> cobaltPosition.setReferenceGCRatio(x.ratio()));
                referenceGCDiploidRatioSelector.select(cobaltPosition).ifPresent(x -> cobaltPosition.setReferenceGCDiploidRatio(x.ratio()));
                tumorGCRatioSelector.select(cobaltPosition).ifPresent(x -> cobaltPosition.setTumorGCRatio(x.ratio()));
            }
        }
    }

    public Multimap<Chromosome, ? extends CobaltPosition> build() {
        return result;
    }

    private ModifiableCobaltPosition create(@NotNull final ReadCount reference) {
        return ModifiableCobaltPosition.create()
                .from(reference)
                .setReferenceReadCount(reference.readCount())
                .setTumorReadCount(0)
                .setReferenceGCRatio(-1)
                .setReferenceGCDiploidRatio(-1)
                .setTumorGCRatio(-1);
    }
}
