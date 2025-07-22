package com.hartwig.hmftools.purple.fittingsnv;

import static java.lang.Math.abs;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;
import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.PaveVcfTags.GNOMAD_FREQ;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.MAPPABILITY_TAG;
import static com.hartwig.hmftools.common.variant.VariantType.INDEL;
import static com.hartwig.hmftools.purple.PurpleConstants.HOTSPOT_GNOMAD_FREQ_THRESHOLD;
import static com.hartwig.hmftools.purple.PurpleConstants.PURITY_INCREMENT_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_TUMOR_ONLY_HOTSPOT_VAF_CUTOFF;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_TUMOR_ONLY_VAF_MAX;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleUtils.formatPurity;
import static com.hartwig.hmftools.purple.PurpleConstants.SNV_FITTING_MAPPABILITY;
import static com.hartwig.hmftools.purple.PurpleConstants.SNV_FITTING_MAX_REPEATS;
import static com.hartwig.hmftools.purple.PurpleConstants.SNV_HOTSPOT_MAX_SNV_COUNT;
import static com.hartwig.hmftools.purple.PurpleConstants.SNV_HOTSPOT_VAF_PROBABILITY;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_TUMOR_ONLY_VAF_MIN;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_TUMOR_ONLY_PLOIDY_MAX;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_TUMOR_ONLY_PLOIDY_MIN;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_FIT_TUMOR_ONLY_PURITY_MIN;
import static com.hartwig.hmftools.purple.fittingsnv.SomaticReadjustmentFit.calcReadjustmentPurity;
import static com.hartwig.hmftools.purple.region.ObservedRegionFactory.EXCLUDED_IMMUNE_REGIONS;

import java.util.Collections;
import java.util.List;
import java.util.Optional;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGenePanel;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.utils.collection.Multimaps;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.filter.HumanChromosomeFilter;
import com.hartwig.hmftools.common.variant.filter.NTFilter;
import com.hartwig.hmftools.common.variant.filter.SGTFilter;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

import com.google.common.collect.Lists;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.filter.CompoundFilter;

public class SomaticPurityFitter
{
    private final int mMinReadCount;
    private final int mMaxReadCount;
    private final double mMinPurity;
    private final double mMaxPurity;

    // kernel density parameters
    private final int mKdMinPeak;
    private final int mKdMinSomatics;

    public SomaticPurityFitter(
            int minPeak, int minSomatics, int minReadCount, int maxReadCount, double minPurity, double maxPurity)
    {
        mMinReadCount = minReadCount;
        mMaxReadCount = maxReadCount;
        mMinPurity = minPurity;
        mMaxPurity = maxPurity;
        mKdMinPeak = minPeak;
        mKdMinSomatics = minSomatics;
    }

    private enum FilterReason
    {
        FILTERED,
        NON_SNV,
        GERMLINE_DIPLOID,
        GERMLINE_ALLELE_COUNT,
        REPEAT_COUNT,
        GNOMAD_FREQ,
        TIER,
        MAX_REPEATS,
        MAPPABILITY;
    }

    public static List<SomaticVariant> findFittingVariants(
            boolean tumorOnlyMode, final List<SomaticVariant> variants, final List<ObservedRegion> observedRegions)
    {
        List<SomaticVariant> fittingVariants = Lists.newArrayList();

        CompoundFilter filter = new CompoundFilter(true);
        filter.add(new SGTFilter());
        filter.add(new HumanChromosomeFilter());
        filter.add(new NTFilter());

        GenomeRegionSelector<ObservedRegion> observedRegionSelector = GenomeRegionSelectorFactory.createImproved(
                Multimaps.fromRegions(observedRegions));

        for(SomaticVariant variant : variants)
        {
            if(!variant.isPass() || !filter.test(variant.context()))
            {
                logFilteredFittingCandidate(variant, "non-passing");
                continue;
            }

            if(variant.type() == VariantType.MNP)
            {
                logFilteredFittingCandidate(variant, "excluded MNV");
                continue;
            }
            else if(tumorOnlyMode && variant.type() == INDEL)
            {
                // only use in tumor-only mode and if not in a repeat context
                if(variant.decorator().repeatCount() > 0)
                {
                    logFilteredFittingCandidate(variant, "invalid indel");
                    continue;
                }
            }

            if(!isFittingCandidate(variant))
                continue;

            /*
            if(variant.type() == INDEL)
            {
                PPL_LOGGER.debug("variant({}) used for fitting with vaf({})", variant, format("%.2f", variant.alleleFrequency()));
            }
            */

            Optional<ObservedRegion> region = observedRegionSelector.select(variant);

            GermlineStatus germlineStatus = region.isPresent() ? region.get().germlineStatus() : GermlineStatus.UNKNOWN;

            if(!variant.isHotspotType() && germlineStatus != GermlineStatus.DIPLOID)
            {
                logFilteredFittingCandidate(variant, "germline not diploid");
                continue;
            }

            fittingVariants.add(variant);
        }

        // PPL_LOGGER.debug("variants({}) used for fitting({})", variants.size(), fittingVariants.size());

        return fittingVariants;
    }

    private static boolean isFittingCandidate(final SomaticVariant variant)
    {
        if(EXCLUDED_IMMUNE_REGIONS.stream().anyMatch(x -> x.containsPosition(variant.chromosome(), variant.position())))
        {
            logFilteredFittingCandidate(variant, "immune region");
            return false;
        }
        
        if(!variant.hasTumorAlleleDepth() || variant.tumorAlleleDepth().TotalReadCount == 0)
        {
            logFilteredFittingCandidate(variant, "zero tumor depth");
            return false;
        }

        VariantTier variantTier = variant.decorator().tier();

        boolean isHotspotType = variant.isHotspotType();
        
        double variantGnomadFreq = variant.context().getAttributeAsDouble(GNOMAD_FREQ, -1);

        if(variantGnomadFreq > 0 && variantGnomadFreq < HOTSPOT_GNOMAD_FREQ_THRESHOLD && isHotspotType)
            return true;

        if(!isHotspotType)
        {
            if(variant.context().hasAttribute(GNOMAD_FREQ))
            {
                logFilteredFittingCandidate(variant, "gnomad frequency");
                return false;
            }

            if(variantTier == VariantTier.LOW_CONFIDENCE || variantTier == VariantTier.UNKNOWN)
            {
                logFilteredFittingCandidate(variant, "low tier");
                return false;
            }

            if(variant.decorator().repeatCount() > SNV_FITTING_MAX_REPEATS)
            {
                logFilteredFittingCandidate(variant, "max repeats");
                return false;
            }

            if(variant.context().hasAttribute(MAPPABILITY_TAG) && variant.decorator().mappability() < SNV_FITTING_MAPPABILITY)
            {
                logFilteredFittingCandidate(variant, "mappability");
                return false;
            }

            if(variant.referenceAlleleReadCount() > 0)
            {
                logFilteredFittingCandidate(variant, "germline allele count");
                return false;
            }
        }

        return true;
    }

    private static void logFilteredFittingCandidate(final SomaticVariant variant, final String reason)
    {
        if(!PPL_LOGGER.isTraceEnabled())
            return;

        PPL_LOGGER.trace("variant({}) excluded from fitting: {}", variant, reason);
    }

    @Nullable
    public FittedPurity fitfromSomatics(
            final List<SomaticVariant> somaticVariants, final List<StructuralVariant> hotspotSVs,
            final List<FittedPurity> diploidCandidates, final List<PurpleCopyNumber> copyNumbers, final Gender gender)
    {
        List<SomaticVariant> filteredSomatics = somaticVariants.stream()
                .filter(x -> x.isHotspotType() || (x.totalReadCount() >= mMinReadCount && x.totalReadCount() <= mMaxReadCount))
                .collect(toList());

        double somaticPeakPurity = 0;
        FittedPurity somaticFitPurity = null;

        if(filteredSomatics.size() >= mKdMinSomatics)
        {
            PPL_LOGGER.info("looking for peak somatic allelic frequencies");

            FittedPurity kdFit = SomaticKernelDensityPeaks.fitPurity(
                    diploidCandidates, filteredSomatics, mKdMinSomatics, mKdMinPeak, mMinPurity, mMaxPurity);

            if(kdFit != null)
            {
                somaticFitPurity = kdFit;
                somaticPeakPurity = kdFit.purity();
                PPL_LOGGER.info("peak somatic purity({})", formatPurity(somaticPeakPurity));
            }
        }
        else
        {
            PPL_LOGGER.info("somatic variants count({}) too low for somatic fit", filteredSomatics.size());
        }

        double hotspotPurity = findHotspotPurity(filteredSomatics, hotspotSVs, somaticPeakPurity, gender);

        if(hotspotPurity > somaticPeakPurity)
        {
            FittedPurity matchedFittedPurity = findMatchedFittedPurity(hotspotPurity, diploidCandidates);

            if(matchedFittedPurity != null)
                return matchedFittedPurity;
        }

        if(somaticFitPurity != null)
        {
            double reassessmentPurity = calcReadjustmentPurity(filteredSomatics, somaticFitPurity, copyNumbers);

            if(reassessmentPurity > somaticPeakPurity)
            {
                FittedPurity matchedFittedPurity = findMatchedFittedPurity(reassessmentPurity, diploidCandidates);

                if(matchedFittedPurity != null)
                    return matchedFittedPurity;
            }
        }

        return somaticFitPurity;
    }

    public static boolean highlyDiploidSomaticOrPanel(final FittedPurity normalPurityFit, final boolean highlyDiploidByFitPurity)
    {
        return (normalPurityFit.purity() > SOMATIC_FIT_TUMOR_ONLY_PURITY_MIN
            && normalPurityFit.ploidy() > SOMATIC_FIT_TUMOR_ONLY_PLOIDY_MIN
            && normalPurityFit.ploidy() < SOMATIC_FIT_TUMOR_ONLY_PLOIDY_MAX) || highlyDiploidByFitPurity;
    }

    protected static FittedPurity findMatchedFittedPurity(double purity, final List<FittedPurity> allCandidates)
    {
        // find the closest purity with diploid ploidy
        FittedPurity closestPurity = null;
        double closestDiff = 0;
        double purityEpsilon = PURITY_INCREMENT_DEFAULT * 0.25;

        for(FittedPurity fittedPurity : allCandidates)
        {
            if(abs(fittedPurity.ploidy() - 2) > 0.005)
                continue;

            double diff = abs(fittedPurity.purity() - purity);

            if(closestPurity == null || diff < closestDiff)
            {
                if(diff < purityEpsilon)
                    return fittedPurity;

                closestDiff = diff;
                closestPurity = fittedPurity;
            }
        }

        return closestPurity;
    }

    @Nullable
    public FittedPurity fitFromSomaticsOnly(
            final DriverGenePanel driverGenes, final List<SomaticVariant> variants, final List<FittedPurity> allCandidates)
    {
        List<Double> variantVafs = Lists.newArrayList();

        for(SomaticVariant variant : variants)
        {
            if(!variant.isHotspotType())
            {
                if(variant.variantImpact() == null)
                    continue;

                CodingEffect codingEffect = variant.variantImpact().CanonicalCodingEffect;

                if(codingEffect != NONSENSE_OR_FRAMESHIFT && codingEffect != MISSENSE)
                    continue;

                DriverGene driverGene = driverGenes.driverGenes().stream()
                        .filter(x -> x.gene().equals(variant.variantImpact().GeneName)).findFirst().orElse(null);

                if(driverGene == null)
                    continue;

                if(codingEffect == NONSENSE_OR_FRAMESHIFT && !driverGene.reportNonsenseAndFrameshift())
                    continue;
                else if(codingEffect == MISSENSE && !driverGene.reportMissenseAndInframe())
                    continue;
            }

            double vaf = variant.alleleFrequency();

            if(variant.isHotspotType())
            {
                vaf = min(vaf, SOMATIC_FIT_TUMOR_ONLY_HOTSPOT_VAF_CUTOFF);
            }
            else
            {
                if(vaf < SOMATIC_FIT_TUMOR_ONLY_VAF_MIN || vaf > SOMATIC_FIT_TUMOR_ONLY_VAF_MAX)
                    continue;
            }

            variantVafs.add(vaf);
        }

        if(variantVafs.isEmpty())
            return null;

        double vaf75thPercentile = calc75thPercentileValue(variantVafs);

        double somaticPurity = vaf75thPercentile * 2;

        PPL_LOGGER.info("somatic VAF-based purity({}) from {} variants", formatPurity(somaticPurity), variantVafs.size());

        FittedPurity matchedFittedPurity = findMatchedFittedPurity(somaticPurity, allCandidates);

        if(matchedFittedPurity != null)
            return matchedFittedPurity;

        return null;
    }

    @VisibleForTesting
    protected static double calc75thPercentileValue(final List<Double> variantVafs)
    {
        if(variantVafs.isEmpty())
            return 0;

        if(variantVafs.size() == 1)
            return variantVafs.get(0);

        Collections.sort(variantVafs);

        if(variantVafs.size() <= 3)
        {
            // take the highest if not at the max value for hotspots
            int topIndex = variantVafs.size() - 1;

            if(variantVafs.get(topIndex) < SOMATIC_FIT_TUMOR_ONLY_HOTSPOT_VAF_CUTOFF)
            {
                return variantVafs.get(topIndex);
            }
            else
            {
                return variantVafs.get(max(topIndex - 1, 0));
            }
        }

        // double index75thPercentile = 0.75 * (variantVafs.size() + 1) - 1;
        double index75thPercentile = 0.75 * (variantVafs.size() - 1);

        int lowerIndex = min((int)floor(index75thPercentile), variantVafs.size() - 2);
        int upperIndex = lowerIndex + 1;

        double lowerVaf = variantVafs.get(lowerIndex);
        double upperVaf = variantVafs.get(upperIndex);

        double indexFraction = index75thPercentile - lowerIndex;
        return lowerVaf + indexFraction * (upperVaf - lowerVaf);
    }

    private double findHotspotPurity(
            final List<SomaticVariant> somaticVariants, final List<StructuralVariant> hotspotSVs,
            final double somaticPeakPurity, final Gender gender)
    {
        // check for a hotspot variant with a higher VAF
        if(somaticVariants.size() > SNV_HOTSPOT_MAX_SNV_COUNT)
            return 0;

        double maxHotspotVaf = 0;

        for(SomaticVariant variant : somaticVariants)
        {
            if(!variant.isHotspotType())
                continue;

            HumanChromosome chromosome = HumanChromosome.fromString(variant.chromosome());

            if(!chromosome.isAutosome())
            {
                if(!(gender == Gender.FEMALE && chromosome == HumanChromosome._X))
                    continue;
            }

            if(variant.alleleFrequency() * 2 <= somaticPeakPurity || variant.alleleFrequency() > 0.5)
                continue;

            if(variant.alleleFrequency() < maxHotspotVaf)
                continue;

            // test this variants allele read count vs what's expected from the somatic peak
            if(!belowRequiredProbability(somaticPeakPurity, variant.totalReadCount(), variant.alleleReadCount()))
                continue;

            PPL_LOGGER.info(format("hotspot(%s:%d) vaf(%.3f %d/%d)",
                    variant.chromosome(), variant.position(),
                    variant.alleleFrequency(), variant.alleleReadCount(), variant.totalReadCount()));

            maxHotspotVaf = max(variant.alleleFrequency(), maxHotspotVaf);
        }

        for(StructuralVariant sv : hotspotSVs)
        {
            if(!sv.hotspot() || sv.isFiltered() || sv.end() == null)
                continue;

            double alleleFrequency = min(sv.start().alleleFrequency(), sv.end().alleleFrequency());

            if(alleleFrequency < maxHotspotVaf)
                continue;

            PPL_LOGGER.info(String.format("hotspotSV(%s %s:%d-%s:%d) vaf(%.3f)",
                    sv.type(), sv.chromosome(true), sv.position(true), sv.chromosome(false), sv.position(false),
                    alleleFrequency));

            maxHotspotVaf = alleleFrequency;
        }

        return maxHotspotVaf * 2;
    }

    private boolean belowRequiredProbability(double peakPurity, int totalReadCount, int alleleReadCount)
    {
        if(peakPurity <= 0)
            return true;

        double expectedAlleleReadCount = peakPurity * 0.5 * totalReadCount;

        PoissonDistribution poissonDist = new PoissonDistribution(expectedAlleleReadCount);
        double probability = 1 - poissonDist.cumulativeProbability(alleleReadCount - 1);

        return probability < SNV_HOTSPOT_VAF_PROBABILITY;
    }
}
