package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.sv.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;

import org.jetbrains.annotations.NotNull;

public final class GeneFusions {

    private GeneFusions() {
    }

    @NotNull
    public static List<LinxFusion> sort(@NotNull List<LinxFusion> fusions) {
        return fusions.stream().sorted((fusion1, fusion2) -> {
            if (fusion1.likelihood() == fusion2.likelihood()) {
                if (fusion1.geneStart().equals(fusion2.geneStart())) {
                    return fusion1.geneEnd().compareTo(fusion2.geneEnd());
                } else {
                    return fusion1.geneStart().compareTo(fusion2.geneStart());
                }
            } else {
                return fusion1.likelihood() == FusionLikelihoodType.HIGH ? -1 : 1;
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    public static Set<String> uniqueGeneFusions(@NotNull List<LinxFusion> fusions) {
        Set<String> genes = Sets.newHashSet();
        for (LinxFusion fusion : fusions) {
            genes.add(name(fusion));
        }
        return genes;
    }

    @NotNull
    public static String name(@NotNull LinxFusion fusion) {
        return fusion.geneStart() + " - " + fusion.geneEnd();
    }

    @NotNull
    public static String transcriptUrl(@NotNull String transcriptField) {
        return "http://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=" + transcriptField;
    }
}
