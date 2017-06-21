package com.hartwig.hmftools.common.zipper;

import java.util.List;

import com.hartwig.hmftools.common.chromosome.Chromosomes;
import com.hartwig.hmftools.common.region.GenomeRegion;

public class RegionZipper {

    public static <S extends GenomeRegion, T extends GenomeRegion> void zip(List<S> primary, List<T> secondary,
            RegionZipperHandler<S, T> handler) {

        String chromosome = "";

        int i = 0, j = 0;
        while (i < primary.size() || j < secondary.size()) {

            S leftRegion = i < primary.size() ? primary.get(i) : null;
            T rightRegion = j < secondary.size() ? secondary.get(j) : null;

            if (leftRegion == null || (rightRegion != null && compare(leftRegion, rightRegion) > 0)) {
                if (!rightRegion.chromosome().equals(chromosome)) {
                    chromosome = rightRegion.chromosome();
                    handler.chromosome(chromosome);
                }
                handler.secondary(rightRegion);
                j++;
            } else {
                if (!leftRegion.chromosome().equals(chromosome)) {
                    chromosome = leftRegion.chromosome();
                    handler.chromosome(chromosome);
                }
                handler.primary(leftRegion);
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
