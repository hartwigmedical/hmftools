package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.peach.PeachGenotype;

import org.jetbrains.annotations.NotNull;

public class Peach {

    private Peach(){

    }

    @NotNull
    public static List<PeachGenotype> sort(@NotNull List<PeachGenotype> peachGenotypes) {
        return peachGenotypes.stream().sorted((genotype1, genotype2) -> {
                if (genotype1.gene().equals(genotype2.gene())) {
                    return genotype1.haplotype().compareTo(genotype2.haplotype());
                } else {
                    return genotype1.gene().compareTo(genotype2.gene());
                }

        }).collect(Collectors.toList());
    }
}
