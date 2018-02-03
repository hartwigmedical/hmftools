package com.hartwig.hmftools.common.gene;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimaps;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.region.hmfslicer.HmfExonRegion;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.zipper.RegionZipper;
import com.hartwig.hmftools.common.zipper.RegionZipperHandler;

import org.jetbrains.annotations.NotNull;

public class GeneCopyNumberFactory implements RegionZipperHandler<PurpleCopyNumber, HmfExonRegion> {

    @NotNull
    private final GeneCopyNumberBuilder builder;

    @NotNull
    public static List<GeneCopyNumber> geneCopyNumbers(@NotNull final List<HmfGenomeRegion> genes,
            @NotNull final List<PurpleCopyNumber> copyNumbers, final List<PurityAdjustedSomaticVariant> enrichedSomatics) {

        final ListMultimap<String, PurityAdjustedSomaticVariant> geneVariants = Multimaps.index(enrichedSomatics, SomaticVariant::gene);

        final List<GeneCopyNumber> result = Lists.newArrayList();
        for (HmfGenomeRegion gene : genes) {
            final List<PurityAdjustedSomaticVariant> somatics =
                    geneVariants.containsKey(gene.gene()) ? geneVariants.get(gene.gene()) : Collections.EMPTY_LIST;

            GeneCopyNumber geneCopyNumber =  new GeneCopyNumberFactory(gene, copyNumbers, somatics).geneCopyNumber();
            if (geneCopyNumber.totalRegions() > 0) {
                result.add(geneCopyNumber);
            }

        }
        return result;
    }

    private GeneCopyNumberFactory(@NotNull final HmfGenomeRegion gene, @NotNull final List<PurpleCopyNumber> copyNumbers,
            final List<PurityAdjustedSomaticVariant> enrichedSomatics) {
        builder = new GeneCopyNumberBuilder(gene);
        RegionZipper.zip(copyNumbers, gene.exome(), this);
        addSomatics(enrichedSomatics);
    }

    private void addSomatics(@NotNull final List<PurityAdjustedSomaticVariant> somatics) {
        somatics.forEach(builder::addSomatic);
    }

    @NotNull
    private GeneCopyNumber geneCopyNumber() {
        return builder.build();
    }

    @Override
    public void enterChromosome(@NotNull final String chromosome) {
        // IGNORE
    }

    @Override
    public void primary(@NotNull final PurpleCopyNumber copyNumber) {
        builder.addCopyNumber(copyNumber);
    }

    @Override
    public void secondary(@NotNull final HmfExonRegion exon) {
        builder.addExon(exon);
    }
}
