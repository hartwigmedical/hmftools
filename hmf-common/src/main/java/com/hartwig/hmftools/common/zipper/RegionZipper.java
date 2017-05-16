package com.hartwig.hmftools.common.zipper;

import java.util.List;

import com.hartwig.hmftools.common.chromosome.Chromosomes;
import com.hartwig.hmftools.common.region.GenomeRegion;

public class RegionZipper {

    public static <S extends GenomeRegion, T extends GenomeRegion> void zip(List<S> left, List<T> right,
            RegionZipperHandler<S, T> handler) {

        int i = 0, j = 0;
        while (i < left.size() || j < right.size()) {

            S myLeft = i < left.size() ? left.get(i) : null;
            T myRight = j < right.size() ? right.get(j) : null;

            if (myLeft == null || compare(myLeft, myRight) > 0) {
                handler.right(myRight);
                j++;
            } else {
                handler.left(myLeft);
                i++;
            }
        }
    }

    private static int compare(GenomeRegion position, GenomeRegion region) {
        int positionChromosome = Chromosomes.asInt(position.chromosome());
        int regionChromosome = Chromosomes.asInt(region.chromosome());
        if (positionChromosome < regionChromosome) {
            return -1;
        }
        if (positionChromosome > regionChromosome) {
            return 1;
        }

        if (position.start() < region.start()) {
            return -1;
        }

        if (position.start() > region.start()) {
            return 1;
        }

        return 0;
    }

}
