package com.hartwig.hmftools.patientreporter.cfreport.data;

import com.hartwig.hmftools.patientreporter.report.data.GeneFusionDataSource;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneFusion;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import static com.hartwig.hmftools.common.fusions.KnownFusionsModel.*;
import static com.hartwig.hmftools.common.fusions.KnownFusionsModel.CIVIC;

public class GeneFusions {

    @NotNull
    public  static String[] geneFusions(@NotNull List<ReportableGeneFusion> fusions) {
        final List<String> returnVariants = new ArrayList<>();
        for (ReportableGeneFusion fusion : GeneFusionDataSource.sort(fusions)) {
            returnVariants.add(GeneFusionDataSource.name(fusion));
        }
        return returnVariants.toArray(new String[0]);
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
    public static String getName(@NotNull ReportableGeneFusion fusion) {
        return fusion.geneStart() + " - " + fusion.geneEnd();
    }

    @NotNull
    public static String getPloidyToCopiesString(@Nullable Double ploidy, boolean hasReliablePurityFit) {

        if (!hasReliablePurityFit) {
            return Util.NAString;
        } else {
            return ploidy != null ? String.format("%.1f", ploidy) : Strings.EMPTY;
        }

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
