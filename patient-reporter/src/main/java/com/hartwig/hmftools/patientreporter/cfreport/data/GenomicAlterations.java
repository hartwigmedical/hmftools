package com.hartwig.hmftools.patientreporter.cfreport.data;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.patientreporter.report.data.GeneCopyNumberDataSource;
import com.hartwig.hmftools.patientreporter.report.data.GeneFusionDataSource;
import com.hartwig.hmftools.patientreporter.report.data.SomaticVariantDataSource;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneFusion;
import com.hartwig.hmftools.patientreporter.variants.ReportableSomaticVariant;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * All methods to handle Genomic Alteration data
 */
public class GenomicAlterations {

    @NotNull
    public static String[] somaticVariantsWithDriver(@NotNull List<ReportableSomaticVariant> variants) {
        final List<String> returnVariants = new ArrayList<>();
        for (final ReportableSomaticVariant variant : SomaticVariantDataSource.sort(variants)) {
            if (SomaticVariantDataSource.driverField(variant).equals("High")) {
                returnVariants.add(variant.gene());
            }
        }
        return returnVariants.toArray(new String[0]);
    }

    public static int countSomaticVariants(@NotNull List<ReportableSomaticVariant> variants) {
        return variants.size();
    }

    @NotNull
    public static String[] amplificationGenes(@NotNull List<GeneCopyNumber> copyNumbers) {
        final List<String> returnVariants = new ArrayList<>();
        for (GeneCopyNumber copyNumber : copyNumbers) {
            if (GeneCopyNumberDataSource.type(copyNumber).equals("gain")) {
                returnVariants.add(copyNumber.gene());
            }
        }
        return returnVariants.toArray(new String[0]);
    }

    @NotNull
    public static String[] lossGenes(@NotNull List<GeneCopyNumber> copyNumbers) {
        final List<String> returnVariants = new ArrayList<>();
        Set<String> geneCopyNumbersLoss = Sets.newHashSet();
        for (GeneCopyNumber copyNumber : copyNumbers) {
            if (GeneCopyNumberDataSource.type(copyNumber).equals("full loss") || GeneCopyNumberDataSource.type(copyNumber)
                    .equals("partial loss")) {
                returnVariants.add(copyNumber.gene());
            }
        }
        return returnVariants.toArray(new String[0]);
    }

    @NotNull
    public  static String[] geneFusions(@NotNull List<ReportableGeneFusion> fusions) {
        final List<String> returnVariants = new ArrayList<>();
        for (ReportableGeneFusion fusion : GeneFusionDataSource.sort(fusions)) {
            returnVariants.add(GeneFusionDataSource.name(fusion));
        }
        return returnVariants.toArray(new String[0]);
    }

}
