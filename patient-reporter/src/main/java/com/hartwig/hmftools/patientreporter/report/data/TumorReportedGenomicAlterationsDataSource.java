package com.hartwig.hmftools.patientreporter.report.data;


import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneDisruption;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneFusion;
import com.hartwig.hmftools.patientreporter.variants.ReportableSomaticVariant;

import org.jetbrains.annotations.NotNull;

public class TumorReportedGenomicAlterationsDataSource {

    private TumorReportedGenomicAlterationsDataSource() {
    }

    @NotNull
    public static String countSomaticVariants(@NotNull List<ReportableSomaticVariant> variants) {
        int countSomaticVariants = 0;
        for (final ReportableSomaticVariant variant : SomaticVariantDataSource.sort(variants)) {
            countSomaticVariants++;
        }
        return Integer.toString(countSomaticVariants);
    }

    @NotNull
    public static String somaticVariantsWithDriver(@NotNull List<ReportableSomaticVariant> variants) {
        List<String> somaticVariants = Lists.newArrayList();
        for (final ReportableSomaticVariant variant : SomaticVariantDataSource.sort(variants)) {
            if (SomaticVariantDataSource.driverField(variant).equals("High")) {
                somaticVariants.add(variant.gene());
            }
        }
        return String.join(", ", somaticVariants);
    }

    @NotNull
    public static String amplificationGenes(@NotNull final List<GeneCopyNumber> copyNumbers) {
        List<String> geneCopyNumbersAmplification = Lists.newArrayList();
        for (final GeneCopyNumber copyNumber : copyNumbers) {
            if (GeneCopyNumberDataSource.type(copyNumber).equals("gain")) {
                geneCopyNumbersAmplification.add(copyNumber.gene());
            }
        }
        return String.join(", ", geneCopyNumbersAmplification);
    }

    @NotNull
    public static String lossGenes(@NotNull final List<GeneCopyNumber> copyNumbers) {
        List<String> geneCopyNumbersLoss = Lists.newArrayList();
        for (final GeneCopyNumber copyNumber : copyNumbers) {
            if (GeneCopyNumberDataSource.type(copyNumber).equals("full loss") || GeneCopyNumberDataSource.type(copyNumber)
                    .equals("partial loss")) {
                geneCopyNumbersLoss.add(copyNumber.gene());
            }
        }
        return String.join(", ", geneCopyNumbersLoss);
    }

    @NotNull
    public static String geneFusions(@NotNull List<ReportableGeneFusion> fusions) {
        List<String> geneFusions = Lists.newArrayList();
        for (ReportableGeneFusion fusion : GeneFusionDataSource.sort(fusions)) {
            geneFusions.add(GeneFusionDataSource.name(fusion));
        }
        return String.join(", ", geneFusions);
    }

    @NotNull
    public static String geneDisruptions(@NotNull List<ReportableGeneDisruption> disruptions) {
        List<String> geneDisruptions = Lists.newArrayList();
        for (ReportableGeneDisruption disruption : GeneDisruptionDataSource.sort(disruptions)) {
            geneDisruptions.add(disruption.gene());
        }
        return String.join(", ", geneDisruptions);
    }
}

