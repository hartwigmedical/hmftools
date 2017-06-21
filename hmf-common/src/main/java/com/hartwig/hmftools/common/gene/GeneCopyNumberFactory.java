package com.hartwig.hmftools.common.gene;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.zipper.RegionZipper;
import com.hartwig.hmftools.common.zipper.RegionZipperHandler;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GeneCopyNumberFactory implements RegionZipperHandler<PurpleCopyNumber, HmfGenomeRegion> {

    public static List<GeneCopyNumber> geneCopyNumbers(List<HmfGenomeRegion> genes, List<PurpleCopyNumber> copyNumbers) {
        GeneCopyNumberFactory handler = new GeneCopyNumberFactory();
        RegionZipper.zip(copyNumbers, genes, handler);
        return handler.geneCopyNumbers();
    }

    @Nullable
    private PurpleCopyNumber current;
    @Nullable
    private GeneCopyNumberBuilder builder;
    @NotNull
    private final List<GeneCopyNumber> geneCopyNumbers;

    private GeneCopyNumberFactory() {
        geneCopyNumbers = Lists.newArrayList();
    }

    private List<GeneCopyNumber> geneCopyNumbers() {
        if (builder != null) {
            finialiseBuilder();
        }
        return geneCopyNumbers;
    }

    @Override
    public void chromosome(@NotNull final String chromosome) {
        finialiseBuilder();
        current = null;
    }

    @Override
    public void primary(@NotNull final PurpleCopyNumber region) {
        current = region;
        if (builder != null) {
            builder.addCopyNumber(region);
        }
    }

    @Override
    public void secondary(@NotNull final HmfGenomeRegion region) {
        finialiseBuilder();
        builder = new GeneCopyNumberBuilder(region);
        if (current != null) {
            builder.addCopyNumber(current);
        }
    }

    private void finialiseBuilder() {
        if (builder != null) {
            geneCopyNumbers.add(builder.build());
            builder = null;
        }
    }
}
