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
import static com.hartwig.hmftools.purple.PurpleConstants.HOTSPOT_GNOMAD_FREQ_THRESHOLD;
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
import java.util.StringJoiner;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.pathogenic.PathogenicSummaryFactory;
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

import org.apache.commons.lang3.ObjectUtils;
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

    public static List<SomaticVariant> findFittingVariants(final List<SomaticVariant> variants, final List<ObservedRegion> observedRegions)
    {
        List<SomaticVariant> fittingVariants = Lists.newArrayList();

        CompoundFilter filter = new CompoundFilter(true);
        filter.add(new SGTFilter());
        filter.add(new HumanChromosomeFilter());
        filter.add(new NTFilter());

        GenomeRegionSelector<ObservedRegion> observedRegionSelector = GenomeRegionSelectorFactory.createImproved(
                Multimaps.fromRegions(observedRegions));

        final int[] filterCounts = new int[FilterReason.values().length];

        for(SomaticVariant variant : variants)
        {
            if(variant.type() != VariantType.SNP)
            {
                ++filterCounts[FilterReason.NON_SNV.ordinal()];
                continue;
            }

            if(!variant.isPass() || !filter.test(variant.context()))
            {
                ++filterCounts[FilterReason.FILTERED.ordinal()];
                continue;
            }

            if(!isFittingCandidate(variant, filterCounts))
                continue;

            Optional<ObservedRegion> region = observedRegionSelector.select(variant);

            GermlineStatus germlineStatus = region.isPresent() ? region.get().germlineStatus() : GermlineStatus.UNKNOWN;

            if(!variant.isHotspotType() && germlineStatus != GermlineStatus.DIPLOID)
            {
                ++filterCounts[FilterReason.GERMLINE_DIPLOID.ordinal()];
                continue;
            }

            fittingVariants.add(variant);
        }

        if(PPL_LOGGER.isDebugEnabled())
        {
            StringJoiner filterCountsStr = new StringJoiner(", ");
            for(FilterReason reason : FilterReason.values())
            {
                filterCountsStr.add(format("%s=%d", reason, filterCounts[reason.ordinal()]));
            }

            PPL_LOGGER.debug("variants({}) fitting({}) filters: {}", variants.size(), fittingVariants.size(), filterCountsStr);
        }

        return fittingVariants;
    }

    private static boolean isFittingCandidate(final SomaticVariant variant, final int[] filterCounts)
    {
        if(EXCLUDED_IMMUNE_REGIONS.stream().anyMatch(x -> x.containsPosition(variant.chromosome(), variant.position())))
            return false;
        
        if(!variant.hasTumorAlleleDepth() || variant.tumorAlleleDepth().TotalReadCount == 0)
            return false;

        VariantTier variantTier = variant.decorator().tier();

        boolean isHotspotType = variant.isHotspotType();
        
        double variantGnomadFreq = variant.context().getAttributeAsDouble(GNOMAD_FREQ, -1);

        if(variantGnomadFreq > 0 && variantGnomadFreq < HOTSPOT_GNOMAD_FREQ_THRESHOLD && isHotspotType)
            return true;

        if(!isHotspotType)
        {
            if(variant.context().hasAttribute(GNOMAD_FREQ))
            {
                ++filterCounts[FilterReason.GNOMAD_FREQ.ordinal()];
                return false;
            }

            if(variantTier == VariantTier.LOW_CONFIDENCE || variantTier == VariantTier.UNKNOWN)
            {
                ++filterCounts[FilterReason.TIER.ordinal()];
                return false;
            }

            if(variant.decorator().repeatCount() > SNV_FITTING_MAX_REPEATS)
            {
                ++filterCounts[FilterReason.MAX_REPEATS.ordinal()];
                return false;
            }

            if(variant.context().hasAttribute(MAPPABILITY_TAG) && variant.decorator().mappability() < SNV_FITTING_MAPPABILITY)
            {
                ++filterCounts[FilterReason.MAPPABILITY.ordinal()];
                return false;
            }

            if(variant.referenceAlleleReadCount() > 0)
            {
                ++filterCounts[FilterReason.GERMLINE_ALLELE_COUNT.ordinal()];
                return false;
            }
        }

        return true;
    }

    @Nullable
    public FittedPurity fromSomatics(
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

    public static boolean highlyDiploidSomaticOrPanel(final FittedPurity normalPurityFit)
    {
        return normalPurityFit.purity() > SOMATIC_FIT_TUMOR_ONLY_PURITY_MIN
            && normalPurityFit.ploidy() > SOMATIC_FIT_TUMOR_ONLY_PLOIDY_MIN
            && normalPurityFit.ploidy() < SOMATIC_FIT_TUMOR_ONLY_PLOIDY_MAX;
    }

    protected static FittedPurity findMatchedFittedPurity(double purity, final List<FittedPurity> allCandidates)
    {
        // find the closest purity with diploid ploidy
        FittedPurity closestPurity = null;
        double closestDiff = 0;

        for(FittedPurity fittedPurity : allCandidates)
        {
            if(abs(fittedPurity.ploidy() - 2) > 0.005)
                continue;

            double diff = abs(fittedPurity.purity() - purity);

            if(closestPurity == null || diff < closestDiff)
            {
                closestDiff = diff;
                closestPurity = fittedPurity;
            }
        }

        return closestPurity;
    }

    @Nullable
    public FittedPurity fitFromSomatics(
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
                if(vaf > SOMATIC_FIT_TUMOR_ONLY_HOTSPOT_VAF_CUTOFF)
                {
                    vaf = SOMATIC_FIT_TUMOR_ONLY_HOTSPOT_VAF_CUTOFF;
                    variantVafs.add(vaf);
                }
                else
                {
                    variantVafs.add(vaf);
                }
            }
            else
            {
                if(vaf >= SOMATIC_FIT_TUMOR_ONLY_VAF_MIN && vaf <= SOMATIC_FIT_TUMOR_ONLY_VAF_MAX)
                {
                    variantVafs.add(vaf);
                }
            }
        }

        if(variantVafs.isEmpty())
            return null;

        Collections.sort(variantVafs);

        double vaf75thPercentile;

        if(variantVafs.size() >= 4)
        {
            double index75thPercentile = 0.75 * (variantVafs.size() + 1);

            int lowerIndex = min((int)floor(index75thPercentile), variantVafs.size() - 2);
            int upperIndex = lowerIndex + 1;

            double lowerVaf = variantVafs.get(lowerIndex);
            double upperVaf = variantVafs.get(upperIndex);

            vaf75thPercentile = lowerVaf + (index75thPercentile - lowerIndex) * (upperVaf - lowerVaf);
        }
        else if(variantVafs.size() == 3)
        {
            if(Collections.max(variantVafs) < SOMATIC_FIT_TUMOR_ONLY_HOTSPOT_VAF_CUTOFF)
            {
                vaf75thPercentile = variantVafs.get(2);
            }
            else
            {
                vaf75thPercentile = variantVafs.get(1);
            }
        }
        else if(variantVafs.size() == 2)
        {
            vaf75thPercentile = (variantVafs.get(0) + variantVafs.get(1)) * 0.5;
        }
        else
        {
            vaf75thPercentile = variantVafs.get(0);
        }

        double somaticPurity = vaf75thPercentile * 2;

        PPL_LOGGER.info("somatic VAF-based purity({}) from {} variants", formatPurity(somaticPurity), variantVafs.size());

        FittedPurity matchedFittedPurity = findMatchedFittedPurity(somaticPurity, allCandidates);

        if(matchedFittedPurity != null)
            return matchedFittedPurity;

        return null;
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
