package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.fusion.ReportableGeneFusion;

import org.jetbrains.annotations.NotNull;

public final class GeneFusions {

    private GeneFusions() {
    }

    @NotNull
    public static List<ReportableGeneFusion> sort(@NotNull List<ReportableGeneFusion> fusions) {
        return fusions.stream().sorted((fusion1, fusion2) -> {
            if (fusion1.geneStart().equals(fusion2.geneStart())) {
                return fusion1.geneEnd().compareTo(fusion2.geneEnd());
            } else {
                return fusion1.geneStart().compareTo(fusion2.geneStart());
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    public static Set<String> uniqueGeneFusions(@NotNull List<ReportableGeneFusion> fusions) {
        Set<String> genes = Sets.newHashSet();
        for (ReportableGeneFusion fusion : fusions) {
            genes.add(name(fusion));
        }
        return genes;
    }

    @NotNull
    public static String name(@NotNull ReportableGeneFusion fusion) {
        return fusion.geneStart() + " - " + fusion.geneEnd();
    }

    @NotNull
    public static String transcriptUrl(@NotNull String transcriptField) {
        return "http://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=" + transcriptField;
    }
}
