package com.hartwig.hmftools.patientreporter.cfreport;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.patientreporter.report.data.GeneCopyNumberDataSource;
import com.hartwig.hmftools.patientreporter.report.data.GeneFusionDataSource;
import com.hartwig.hmftools.patientreporter.report.data.SomaticVariantDataSource;
import com.hartwig.hmftools.patientreporter.structural.ReportableGeneFusion;
import com.hartwig.hmftools.patientreporter.variants.ReportableSomaticVariant;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.text.DecimalFormat;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

public class DataUtility {

   // Number formatting
    private static final DecimalFormat PERCENTAGE_FORMAT = new DecimalFormat("#'%'");

    // Missing/invalid data indicators
    public static final String NoneString = "NONE";
    public static final String NAString = "N/A";

    /**
     * Remap v from [inMin, inMax] to [0, 100]
     */
    public static double mapPercentage(final double v, final double inMin, final double inMax) {
        return map(v, inMin, inMax, 0, 100);
    }

    /**
     * Remap v from [inMin, inMax] to [outMin, outMax]
     */
    public static double map(final double v, final double inMin, final double inMax, final double outMin, final double outMax) {
        return (v - inMin) * (outMax - outMin) / (inMax - inMin) + outMin;
    }

    @NotNull
    public static String formatPercentage(final double percentage) {
        return PERCENTAGE_FORMAT.format(percentage);
    }

    @NotNull
    public static String formatDate(@Nullable final LocalDate date) {
        final DateTimeFormatter formatter = DateTimeFormatter.ofPattern(ReportResources.DATE_TIME_FORMAT);
        return date != null ? formatter.format(date) : "?";
    }

    /**
     * All members dealing with Tumor Purity
     */
    public static class TumorPurity {
        public static final int RANGE_MIN = 0;
        public static final int RANGE_MAX = 1;
    }

    /**
     * All members dealing with HR Deficiency data
     */
    public static class HrDeficiency {

        public static final double RANGE_MIN = 0;
        public static final double RANGE_MAX = 1;


        /**
         * Interpret HR Deficiency value. It's assumed the purity fit has been checked
         */
        @NotNull
        public static String interpretToString(final double chordHrdScore) {
            return new DecimalFormat("#.##").format(chordHrdScore);
        }

//        public static double graphValue(final double value) {
//
//        }

    }


    /**
     * All members dealing with Tumor Mutational Load
     */
    public static class MutationalLoad {

        public static final int RANGE_MIN = 1;
        public static final int RANGE_MAX = 1000;
        public static final int THRESHOLD = 140;

        /**
         * Interpret mutational load value to be either "High" or "Low". It's assumed the purity
         * fit has been checked
         */
        @NotNull
        public static String interpretToString(final int mutationalLoad) {
            if (mutationalLoad > THRESHOLD) {
                return "High";
            } else {
                return "Low";
            }
        }

        public static double scale(double value) {
            return Math.log10(value);
        }

    }

    public static class MutationalBurden {

        public static final double RANGE_MIN = 1E-2;
        public static final double RANGE_MAX = 120;

    }

    /**
     * All members dealing with Tumor Mutational Load
     */
    public static class MicroSatellite {

        public static final double RANGE_MIN = 1E-2;
        public static final double RANGE_MAX = 100;
        public static final double THRESHOLD = 4;

        /**
         * Interpret micro satellite indel value to be either "Stable" or "Instable". It's assumed the purity
         * fit has been checked
         */
        @NotNull
        public static String interpretToString(double microSatelliteIndelsPerMb) {
            if (microSatelliteIndelsPerMb > THRESHOLD) {
                return "Unstable";
            } else {
                return "Stable";
            }
        }

        public static double scale(double value) {
            return Math.log10(value);
        }

    }

    /**
     * All methods to handle Genomic Alteration data
     */
    public static class GenomicAlterations {

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

}
