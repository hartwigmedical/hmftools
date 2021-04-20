package com.hartwig.hmftools.purple.fitting;

import static java.lang.String.format;

import static com.google.common.collect.Lists.newArrayList;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleCommon.formatDbl;
import static com.hartwig.hmftools.purple.fitting.SomaticHistogramPeaks.calcProbabilityUpperBound;

import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.purity.SomaticPeak;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.clonality.ModifiableWeightedPloidy;
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
                        PPL_LOGGER.debug("Somatic implied purity: {}", impliedPurity);
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
                        PPL_LOGGER.debug("Somatic implied purity: {}", impliedPurity);
                        return diploid;
                    }
                    else
                    {
                        PPL_LOGGER.warn("Unable to find diploid solution for implied purity: {}", impliedPurity);
                    }
                }
            }
        }

        PPL_LOGGER.debug("Unable to determine somatic implied purity.");
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

    private static final int MAX_HOTSPOT_SNV_COUNT = 1000;

    private double findPurityPeak(final List<SomaticVariant> variants)
    {
        final List<ModifiableWeightedPloidy> weightedVAFs = newArrayList();

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

        double maxPurity = SomaticHistogramPeaks.findVafPeak(2, 0.005, 2, weightedVAFs);

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
}
