package com.hartwig.hmftools.purple.fitting;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.google.common.collect.Lists.newArrayList;
import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleCommon.formatDbl;
import static com.hartwig.hmftools.purple.PurpleCommon.formatPurity;
import static com.hartwig.hmftools.purple.config.PurpleConstants.SNV_HOTSPOT_MAX_SNV_COUNT;
import static com.hartwig.hmftools.purple.config.PurpleConstants.SNV_HOTSPOT_VAF_PROBABILITY;
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
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.jetbrains.annotations.NotNull;

public class SomaticPurityFitter
{
    // kernel density parameters
    private final int mKdMinPeak;
    private final int mKdMinSomatics;
    private final double mMinPurity;
    private final double mMaxPurity;

    public SomaticPurityFitter(final int minPeak, final int minSomatics, final double minPurity, final double maxPurity)
    {
        mKdMinPeak = minPeak;
        mKdMinSomatics = minSomatics;
        mMinPurity = minPurity;
        mMaxPurity = maxPurity;
    }

    public Optional<FittedPurity> fromSomatics(
            final List<SomaticVariant> variants, final List<StructuralVariant> structuralVariants, final List<FittedPurity> allCandidates)
    {
        if(variants.size() < mKdMinSomatics)
        {
            PPL_LOGGER.info("somatic variants count({}) too low for somatic fit", variants.size());
            return Optional.empty();
        }

        PPL_LOGGER.info("looking for peak somatic allelic frequencies");

        Optional<FittedPurity> kdFit = SomaticKernelDensityPeaks.fitPurity(
                allCandidates, variants, mKdMinSomatics, mKdMinPeak, mMinPurity, mMaxPurity);

        if(!kdFit.isPresent())
            return kdFit;

        double peakPurity = kdFit.get().purity();

        PPL_LOGGER.info("peak somatic purity({})", formatPurity(peakPurity));

        // check for a hotspot variant with a higher VAF
        int snvCount = (int)variants.stream().filter(x -> !x.isFiltered()).count();

        if(snvCount > SNV_HOTSPOT_MAX_SNV_COUNT)
            return kdFit;

        double maxHotspotVaf = 0;

        for(SomaticVariant variant : variants)
        {
            if(variant.isFiltered() || !variant.isHotspot())
                continue;

            if(!HumanChromosome.contains(variant.chromosome()) || !HumanChromosome.fromString(variant.chromosome()).isAutosome())
                continue;

            if(variant.alleleFrequency() * 2 <= peakPurity || variant.alleleFrequency() > 0.5)
                continue;

            if(variant.alleleFrequency() < maxHotspotVaf)
                continue;

            // test this variants allele read count vs what's expected from the somatic peak
            if(!belowRequiredProbability(peakPurity, variant.totalReadCount(), variant.alleleReadCount()))
                continue;

            PPL_LOGGER.info(String.format("hotspot(%s:%d) vaf(%.3f %d/%d)",
                    variant.chromosome(), variant.position(),
                    variant.alleleFrequency(), variant.alleleReadCount(), variant.totalReadCount()));

            maxHotspotVaf = max(variant.alleleFrequency(), maxHotspotVaf);
        }

        for(StructuralVariant sv : structuralVariants)
        {
            if(!sv.hotspot() || sv.isFiltered() || sv.end() == null)
                continue;

            double alleleFrequency = min(sv.start().alleleFrequency(), sv.end().alleleFrequency());

            if(alleleFrequency < maxHotspotVaf)
                continue;

            /*
            sv.start().tumorVariantFragmentCount()
            if(!belowRequiredProbability(peakPurity, variant.totalReadCount(), variant.alleleReadCount()))
                continue;

            PPL_LOGGER.info(String.format("hotspot(%s:%d) vaf(%.3f %d/%d)",
                    variant.chromosome(), variant.position(),
                    variant.alleleFrequency(), variant.alleleReadCount(), variant.totalReadCount()));

            maxHotspotVaf = max(variant.alleleFrequency(), maxHotspotVaf);
            */
        }

        maxHotspotVaf *= 2;

        if(maxHotspotVaf > peakPurity)
        {
            return Optional.of(ImmutableFittedPurity.builder()
                    .score(1)
                    .diploidProportion(1)
                    .normFactor(1)
                    .purity(maxHotspotVaf)
                    .somaticPenalty(1)
                    .ploidy(2)
                    .build());
        }
        else
        {
            return kdFit;
        }
    }

    private boolean belowRequiredProbability(double peakPurity, int totalReadCount, int alleleReadCount)
    {
        double expectedAlleleReadCount = peakPurity * 0.5 * totalReadCount;

        PoissonDistribution poissonDist = new PoissonDistribution(expectedAlleleReadCount);
        double poissonProb = 1 - poissonDist.cumulativeProbability(alleleReadCount - 1);

        return poissonProb < SNV_HOTSPOT_VAF_PROBABILITY;
    }
}
