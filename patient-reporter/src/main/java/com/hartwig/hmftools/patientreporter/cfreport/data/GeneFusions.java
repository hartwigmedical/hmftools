package com.hartwig.hmftools.patientreporter.cfreport.data;

import com.hartwig.hmftools.patientreporter.report.data.GeneFusionDataSource;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneFusion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import static com.hartwig.hmftools.common.fusions.KnownFusionsModel.*;

public final class GeneFusions {

    private GeneFusions() {
    }

    @NotNull
    public static String[] geneFusions(@NotNull final List<ReportableGeneFusion> fusions) {
        final List<String> returnVariants = new ArrayList<>();
        for (ReportableGeneFusion fusion : GeneFusionDataSource.sort(fusions)) {
            returnVariants.add(GeneFusionDataSource.name(fusion));
        }
        return returnVariants.toArray(new String[0]);
    }

    @NotNull
    public static List<ReportableGeneFusion> sort(@NotNull final List<ReportableGeneFusion> fusions) {
        return fusions.stream().sorted((fusion1, fusion2) -> {
            if (fusion1.geneStart().equals(fusion2.geneStart())) {
                return fusion1.geneEnd().compareTo(fusion2.geneEnd());
            } else {
                return fusion1.geneStart().compareTo(fusion2.geneStart());
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    public static String getName(@NotNull final ReportableGeneFusion fusion) {
        return fusion.geneStart() + " - " + fusion.geneEnd();
    }

    @NotNull
    public static String getSourceUrl(@NotNull String sourceName) {
        switch (sourceName) {
            case ONCOKB:
                return "http://oncokb.org/#/";
            case COSMIC:
                return "https://cancer.sanger.ac.uk/cosmic";
            case CGI:
                return "https://www.cancergenomeinterpreter.org/biomarkers";
            case CIVIC:
                return "https://civicdb.org/browse/somaticVariants";
            default:
                return Strings.EMPTY;
        }
    }
}
