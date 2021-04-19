package com.hartwig.hmftools.purple.fitting;

import static com.google.common.collect.Lists.newArrayList;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleCommon.formatDbl;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.purity.SomaticPeak;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.clonality.ModifiableWeightedPloidy;
import com.hartwig.hmftools.common.variant.clonality.PeakModel;
import com.hartwig.hmftools.purple.config.ConfigSupplier;

import org.jetbrains.annotations.NotNull;

class SomaticPurityFitter
{
    private final ConfigSupplier mConfig;

    private final int minPeak;
    private final int minSomatics;
    private final double minPurity;
    private final double maxPurity;

    SomaticPurityFitter(final ConfigSupplier config,
            final int minPeak, final int minSomatics, final double minPurity, final double maxPurity)
    {
        mConfig = config;
        this.minPeak = minPeak;
        this.minSomatics = minSomatics;
        this.minPurity = minPurity;
        this.maxPurity = maxPurity;
    }

    public Optional<FittedPurity> fromSomatics(@NotNull final List<FittedPurity> allCandidates, @NotNull final List<SomaticVariant> variants)
    {
        if(variants.size() < minSomatics)
        {
            PPL_LOGGER.info("somatic variants count({}) too low for somatic fit", variants.size());
            return Optional.empty();
        }

        PPL_LOGGER.info("Looking for peak somatic allelic frequencies");

        Optional<FittedPurity> kdFit = fitPurityFromKernelDensity(allCandidates, variants);

        final List<SomaticPeak> modelPeaks = findPeaks(variants);

        if(modelPeaks.isEmpty())
        {
            PPL_LOGGER.info("Somatic purity from VAFs failed");
            return Optional.empty();
        }

        double maxPurity = modelPeaks.stream().mapToDouble(x -> x.alleleFrequency()).max().orElse(0);

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
                maxPurity, kdFit.isPresent() ? formatDbl.format(kdFit.get().purity()) : 0);

        return Optional.of(fittedPurity);
    }

    private Optional<FittedPurity> fitPurityFromKernelDensity(
            @NotNull final List<FittedPurity> allCandidates, @NotNull final List<SomaticVariant> variants)
    {
        final List<SomaticPeak> peaks = SomaticKernelDensityPeaks.findSomaticPeaks(variants);

        // First try and get the largest implied purity where peak count > minPeak
        int maxPeak = 0;
        for(int i = peaks.size() - 1; i >= 0; i--)
        {
            SomaticPeak peak = peaks.get(i);
            double impliedPurity = peak.alleleFrequency() * 2;
            if(inPurityRange(impliedPurity))
            {
                if(peak.count() >= minPeak)
                {
                    Optional<FittedPurity> diploid = diploid(impliedPurity, allCandidates);
                    if(diploid.isPresent())
                    {
                        PPL_LOGGER.info("Somatic implied purity: {}", impliedPurity);
                        return diploid;
                    }
                    else
                    {
                        PPL_LOGGER.warn("Unable to find diploid solution for implied purity: {}", impliedPurity);
                    }
                }
                maxPeak = Math.max(maxPeak, peak.count());
            }
        }

        // Failing that, get the implied purity with the largest peak
        if(maxPeak > 0)
        {
            for(int i = peaks.size() - 1; i >= 0; i--)
            {
                SomaticPeak peak = peaks.get(i);
                double impliedPurity = peak.alleleFrequency() * 2;
                if(peak.count() == maxPeak)
                {
                    Optional<FittedPurity> diploid = diploid(impliedPurity, allCandidates);
                    if(diploid.isPresent())
                    {
                        PPL_LOGGER.info("Somatic implied purity: {}", impliedPurity);
                        return diploid;
                    }
                    else
                    {
                        PPL_LOGGER.warn("Unable to find diploid solution for implied purity: {}", impliedPurity);
                    }
                }
            }
        }

        PPL_LOGGER.info("Unable to determine somatic implied purity.");
        return Optional.empty();
    }

    private boolean inPurityRange(double impliedPurity)
    {
        return Doubles.greaterOrEqual(impliedPurity, minPurity) && Doubles.lessOrEqual(impliedPurity, maxPurity);
    }

    private static Optional<FittedPurity> diploid(double purity, List<FittedPurity> diploidCandidates)
    {
        return diploidCandidates.stream().filter(x -> Doubles.equal(x.purity(), purity)).findFirst();
    }

    private List<SomaticPeak> findPeaks(final List<SomaticVariant> variants)
    {
        final List<ModifiableWeightedPloidy> weightedPloidies = newArrayList();

        double clonalityMaxPloidy = mConfig.somaticConfig().clonalityMaxPloidy();

        for(SomaticVariant variant : variants)
        {
            if(Doubles.lessThan(variant.variantCopyNumber(), clonalityMaxPloidy)
                    && variant.filter().equals(SomaticVariantFactory.PASS_FILTER)
                    && HumanChromosome.contains(variant.chromosome()) && HumanChromosome.fromString(variant.chromosome()).isAutosome())
            {
                // AllelicDepth depth = variant.adjustedVAF().allelicDepth(mConfig.commonConfig().tumorSample());

                weightedPloidies.add(ModifiableWeightedPloidy.create()
                        .setPloidy(variant.alleleFrequency())
                        .setAlleleReadCount(variant.alleleReadCount())
                        .setTotalReadCount(variant.totalReadCount())
                        .setWeight(1));
            }
        }

        final SomaticHistogramPeaks somaticHistogramPeaks = new SomaticHistogramPeaks(10, 0.05);

        return somaticHistogramPeaks.model(weightedPloidies);
    }
}
