package com.hartwig.hmftools.common.gene;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.region.hmfslicer.HmfExonRegion;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.zipper.RegionZipper;
import com.hartwig.hmftools.common.zipper.RegionZipperHandler;

import org.jetbrains.annotations.NotNull;

public class GeneCopyNumberFactory implements RegionZipperHandler<PurpleCopyNumber, HmfExonRegion> {

    @NotNull
    private final GeneCopyNumberBuilder builder;

    @NotNull
    public static List<GeneCopyNumber> geneCopyNumbers(List<HmfGenomeRegion> genes, List<PurpleCopyNumber> copyNumbers) {
        return genes.stream().map(x -> new GeneCopyNumberFactory(x, copyNumbers).geneCopyNumber()).collect(Collectors.toList());
    }

    private GeneCopyNumberFactory(@NotNull final HmfGenomeRegion gene, @NotNull final List<PurpleCopyNumber> copyNumbers) {
        builder = new GeneCopyNumberBuilder(gene);
        RegionZipper.zip(copyNumbers, gene.exome(), this);
    }

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
