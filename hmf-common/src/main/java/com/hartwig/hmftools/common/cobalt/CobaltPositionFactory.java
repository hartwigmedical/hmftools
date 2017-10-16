package com.hartwig.hmftools.common.cobalt;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CobaltPositionFactory {

    private static final Logger LOGGER = LogManager.getLogger(CobaltPositionFactory.class);
    private final Multimap<String, ModifiableCobaltPosition> result = ArrayListMultimap.create();

    public CobaltPositionFactory(@NotNull final Multimap<String, ReadCount> referenceCount) {
        for (String chromosomeName : referenceCount.keySet()) {
            if (HumanChromosome.contains(chromosomeName)) {
                for (ReadCount referenceReadCount : referenceCount.get(chromosomeName)) {
                    result.put(chromosomeName, create(referenceReadCount));
                }
            } else {
                LOGGER.info("Excluding unsupported {} chromosome from read ratios", chromosomeName);

            }
        }
    }

    public void addTumorCount(@NotNull final Multimap<String, ReadCount> tumorCount) {
        final GenomePositionSelector<ReadCount> tumorCountSelector = GenomePositionSelectorFactory.create(tumorCount);
        for (String chromosomeName : result.keySet()) {
            for (ModifiableCobaltPosition cobaltPosition : result.get(chromosomeName)) {
                tumorCountSelector.select(cobaltPosition).ifPresent(x -> cobaltPosition.setTumorReadCount(x.readCount()));
            }
        }
    }

    public void addRatios(@NotNull final Multimap<String, ReadRatio> referenceGCRatio,
            @NotNull final Multimap<String, ReadRatio> tumorGCRatio, @NotNull final Multimap<String, ReadRatio> referenceGCDiploidRatio) {

        final GenomePositionSelector<ReadRatio> referenceGCRatioSelector = GenomePositionSelectorFactory.create(referenceGCRatio);
        final GenomePositionSelector<ReadRatio> referenceGCDiploidRatioSelector =
                GenomePositionSelectorFactory.create(referenceGCDiploidRatio);
        final GenomePositionSelector<ReadRatio> tumorGCRatioSelector = GenomePositionSelectorFactory.create(tumorGCRatio);
        for (String chromosomeName : result.keySet()) {
            for (ModifiableCobaltPosition cobaltPosition : result.get(chromosomeName)) {
                referenceGCRatioSelector.select(cobaltPosition).ifPresent(x -> cobaltPosition.setReferenceGCRatio(x.ratio()));
                referenceGCDiploidRatioSelector.select(cobaltPosition).ifPresent(x -> cobaltPosition.setReferenceGCDiploidRatio(x.ratio()));
                tumorGCRatioSelector.select(cobaltPosition).ifPresent(x -> cobaltPosition.setTumorGCRatio(x.ratio()));
            }
        }
    }

    public Multimap<String, ? extends CobaltPosition> build() {
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
