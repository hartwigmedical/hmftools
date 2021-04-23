package com.hartwig.hmftools.purple.fitting;

import static java.lang.String.format;

import static com.google.common.collect.Lists.newArrayList;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleCommon.formatDbl;
import static com.hartwig.hmftools.purple.config.ConfigSupplier.getConfigValue;
import static com.hartwig.hmftools.purple.config.SomaticFitConfig.SOMATIC_MIN_PEAK_DEFAULT;
import static com.hartwig.hmftools.purple.config.SomaticFitConfig.SOMATIC_MIN_VARIANTS_DEFAULT;
import static com.hartwig.hmftools.purple.fitting.SomaticHistogramPeaks.calcProbabilityUpperBound;

import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurity;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

class SomaticPurityFitter
{
    // kernel density parameters
    private final int mKdMinPeak;
    private final int mKdMinSomatics;
    private final double mMinPurity;
    private final double mMaxPurity;

    // histogram parameters
    private final double mBinWidth;
    private final double mMinPeakWeight;
    private final double mMinPeakPerc;
    private final int mPeakWidth;

    // stand-alone analysis
    private final String mSomaticVariantFile;
    private final String mOututDir;

    // config string
    private static final String BIN_WIDTH = "som_fit_bin_width";
    private static final String PEAK_WIDTH = "som_fit_peak_width";
    private static final String MIN_PEAK_PERC = "som_fit_min_peak_perc";
    private static final String MIN_PEAK_WEIGHT = "som_fit_min_peak_weight";

    private static final String SOM_VARIANT_VCF = "somatic_vcf_file";

    private static final int MAX_HOTSPOT_SNV_COUNT = 1000;
    private static final double UPPER_DEPTH_BOUND_PERC = 1.25;
    private static final double LOWER_DEPTH_BOUND_PERC = 0.75;
    private static final double BIN_WIDTH_DEFAULT = 0.005;
    private static final int PEAK_WIDTH_DEFAULT = 2; // 2 either side of the peak's bucket
    private static final double MIN_PEAK_WEIGHT_DEFAULT = 20; // 2 either side of the peak's bucket
    private static final double MIN_PEAK_PERC_DEFAULT = 0.03;

    public SomaticPurityFitter(final CommandLine cmd)
    {
        mBinWidth = getConfigValue(cmd, BIN_WIDTH, BIN_WIDTH_DEFAULT);
        mMinPeakWeight = getConfigValue(cmd, MIN_PEAK_WEIGHT, MIN_PEAK_WEIGHT_DEFAULT);
        mMinPeakPerc = getConfigValue(cmd, MIN_PEAK_PERC, MIN_PEAK_PERC_DEFAULT);
        mPeakWidth = getConfigValue(cmd, PEAK_WIDTH, PEAK_WIDTH_DEFAULT);

        mKdMinPeak = SOMATIC_MIN_PEAK_DEFAULT;
        mKdMinSomatics = SOMATIC_MIN_VARIANTS_DEFAULT;
        mMinPurity = 0.08;
        mMaxPurity = 1.0;

        mOututDir = parseOutputDir(cmd);
        mSomaticVariantFile = cmd.getOptionValue(SOM_VARIANT_VCF);
    }

    public static void addOptions(@NotNull Options options)
    {
        options.addOption(SOM_VARIANT_VCF, true, "Somatic variant VCF");
        options.addOption(OUTPUT_DIR, true, "Output directory");

        options.addOption(BIN_WIDTH, true, "");
        options.addOption(PEAK_WIDTH, true, "");
        options.addOption(MIN_PEAK_PERC, true, "");
        options.addOption(MIN_PEAK_WEIGHT, true, "");
    }

    SomaticPurityFitter(final int minPeak, final int minSomatics, final double minPurity, final double maxPurity)
    {
        mBinWidth = BIN_WIDTH_DEFAULT;
        mMinPeakWeight = MIN_PEAK_WEIGHT_DEFAULT;
        mMinPeakPerc = MIN_PEAK_PERC_DEFAULT;
        mPeakWidth = PEAK_WIDTH_DEFAULT;

        mKdMinPeak = minPeak;
        mKdMinSomatics = minSomatics;
        mMinPurity = minPurity;
        mMaxPurity = maxPurity;

        mOututDir = null;
        mSomaticVariantFile = null;
    }

    public Optional<FittedPurity> fromSomatics(final List<SomaticVariant> variants, final List<FittedPurity> allCandidates)
    {
        if(variants.size() < mKdMinSomatics)
        {
            PPL_LOGGER.info("somatic variants count({}) too low for somatic fit", variants.size());
            return Optional.empty();
        }

        PPL_LOGGER.info("Looking for peak somatic allelic frequencies");

        Optional<FittedPurity> kdFit = SomaticKernelDensityPeaks.fitPurity(
                allCandidates, variants, mKdMinSomatics, mKdMinPeak, mMinPurity, mMaxPurity);

        double maxPurity = findPurityPeak(variants);

        if(maxPurity <= 0)
        {
            PPL_LOGGER.info("Somatic purity from VAFs failed");
            return Optional.empty();
        }

        maxPurity *= 2; // since a diploid sample

        FittedPurity fittedPurity = ImmutableFittedPurity.builder()
                .score(1)
                .diploidProportion(1)
                .normFactor(1)
                .purity(maxPurity)
                .somaticPenalty(1)
                .ploidy(2)
                .build();

        PPL_LOGGER.info("somatic purity fit: peakModel({}) kernelDensity({})",
                formatDbl.format(maxPurity), kdFit.isPresent() ? formatDbl.format(kdFit.get().purity()) : 0);

        return Optional.of(fittedPurity);
    }

    private double findPurityPeak(final List<SomaticVariant> variants)
    {
        final List<WeightedPloidy> weightedVAFs = newArrayList();

        final List<SomaticVariant> hotspotVariants = Lists.newArrayList();

        for(SomaticVariant variant : variants)
        {
            if(!variant.filter().equals(SomaticVariantFactory.PASS_FILTER))
                continue;

            if(variant.type() != VariantType.SNP)
                continue;

            if(!HumanChromosome.contains(variant.chromosome()) || !HumanChromosome.fromString(variant.chromosome()).isAutosome())
                continue;

            double varWeight = 1; // variant.isHotspot() ? HOTSPOT_WEIGHT : 1;

            weightedVAFs.add(ModifiableWeightedPloidy.create()
                    .setPloidy(variant.alleleFrequency())
                    .setAlleleReadCount(variant.alleleReadCount())
                    .setTotalReadCount(variant.totalReadCount())
                    .setWeight(varWeight));

            if(variant.isHotspot())
                hotspotVariants.add(variant);
        }

        double medianTotalReadCount = calcMedianDepth(weightedVAFs);
        double lowerReadCount = medianTotalReadCount * LOWER_DEPTH_BOUND_PERC;
        double upperReadCount = medianTotalReadCount * UPPER_DEPTH_BOUND_PERC;

        PPL_LOGGER.info(String.format("read count median(%.1f up=%.1f down=%.1f)", medianTotalReadCount, lowerReadCount, upperReadCount));

        final List<WeightedPloidy> validVAFs = weightedVAFs.stream()
                .filter(x -> x.totalReadCount() >= lowerReadCount)
                .filter(x -> x.totalReadCount() <= upperReadCount)
                .collect(Collectors.toList());

        SomaticHistogramPeaks somaticHistogram = new SomaticHistogramPeaks(2, mBinWidth, mPeakWidth, mMinPeakWeight, mMinPeakPerc);
        double maxPurity = somaticHistogram.findVafPeak(validVAFs);

        if(maxPurity <= 0)
            return 0;

        if(weightedVAFs.size() <= MAX_HOTSPOT_SNV_COUNT)
        {
            double upperPeak = calcProbabilityUpperBound(weightedVAFs.size(), maxPurity);

            // look for any hotspot with a VAF higher than this bound
            SomaticVariant topPurityVar = null;

            for(SomaticVariant var : hotspotVariants)
            {
                if(var.alleleFrequency() < upperPeak || var.alleleFrequency() > 0.5)
                    continue;

                if(topPurityVar == null || var.alleleFrequency() > topPurityVar.alleleFrequency())
                    topPurityVar = var;
            }

            if(topPurityVar != null)
            {
                PPL_LOGGER.info("purity({}) set from hotspot({}:{}) vs peak({} upper={})",
                        formatDbl.format(topPurityVar.alleleFrequency()), topPurityVar.chromosome(), topPurityVar.position(),
                        formatDbl.format(maxPurity), format("%.4f", upperPeak));

                maxPurity = topPurityVar.alleleFrequency();
            }
        }

        return maxPurity;
    }

    private double calcMedianDepth(final List<WeightedPloidy> weightedVariants)
    {
        final List<Integer> totalReadCounts = Lists.newArrayList();

        for(WeightedPloidy var : weightedVariants)
        {
            int index = 0;
            while(index < totalReadCounts.size())
            {
                if(var.totalReadCount() < totalReadCounts.get(index))
                    break;

                ++index;
            }

            totalReadCounts.add(index, var.totalReadCount());
        }

        int midIndex = totalReadCounts.size() / 2;
        if((totalReadCounts.size() % 2) == 0)
            return (totalReadCounts.get(midIndex - 1) + totalReadCounts.get(midIndex)) * 0.5;
        else
            return totalReadCounts.get(midIndex - 1);
    }

}
