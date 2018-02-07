package com.hartwig.hmftools.common.gene;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimaps;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.zipper.RegionZipper;

import org.jetbrains.annotations.NotNull;

public class GeneCopyNumberFactory {

    @NotNull
    public static List<GeneCopyNumber> geneCopyNumbers(@NotNull final List<HmfGenomeRegion> genes,
            @NotNull final List<PurpleCopyNumber> somaticCopyNumbers, @NotNull final List<PurpleCopyNumber> germlineDeletions,
            @NotNull final List<PurityAdjustedSomaticVariant> enrichedSomatics) {

        final ListMultimap<String, PurityAdjustedSomaticVariant> variantMap = Multimaps.index(enrichedSomatics, SomaticVariant::gene);

        final List<GeneCopyNumber> result = Lists.newArrayList();
        for (HmfGenomeRegion gene : genes) {
            final List<PurityAdjustedSomaticVariant> variants =
                    variantMap.containsKey(gene.gene()) ? variantMap.get(gene.gene()) : Collections.EMPTY_LIST;

            final GeneCopyNumberBuilder builder = new GeneCopyNumberBuilder(gene);
            RegionZipper.zip(somaticCopyNumbers, gene.exome(), builder);
            RegionZipper.zip(germlineDeletions, gene.exome(), builder);
            variants.forEach(builder::somatic);

            GeneCopyNumber geneCopyNumber = builder.build();
            if (geneCopyNumber.totalRegions() > 0) {
                result.add(geneCopyNumber);
            }

        }
        return result;
    }
}
