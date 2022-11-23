package com.hartwig.hmftools.purple.fitting;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.PaveVcfTags.GNOMAD_FREQ;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleUtils.formatPurity;
import static com.hartwig.hmftools.purple.config.PurpleConstants.SNV_FITTING_MAPPABILITY;
import static com.hartwig.hmftools.purple.config.PurpleConstants.SNV_FITTING_MAX_REPEATS;
import static com.hartwig.hmftools.purple.config.PurpleConstants.SNV_HOTSPOT_MAX_SNV_COUNT;
import static com.hartwig.hmftools.purple.config.PurpleConstants.SNV_HOTSPOT_VAF_PROBABILITY;

import java.util.List;
import java.util.Optional;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurity;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.utils.collection.Multimaps;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.filter.HumanChromosomeFilter;
import com.hartwig.hmftools.common.variant.filter.NTFilter;
import com.hartwig.hmftools.common.variant.filter.SGTFilter;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

import org.apache.commons.compress.utils.Lists;
import org.apache.commons.math3.distribution.PoissonDistribution;

import htsjdk.variant.variantcontext.filter.CompoundFilter;

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

            if(germlineStatus != GermlineStatus.DIPLOID)
            {
                ++filterCounts[FilterReason.GERMLINE_DIPLOID.ordinal()];
                continue;
            }

            PPL_LOGGER.trace("variant({}) germlineStatus({}) used for fitting", variant.toString(), germlineStatus);

            fittingVariants.add(variant);
        }

        StringJoiner filterCountsStr = new StringJoiner(", ");
        for(FilterReason reason : FilterReason.values())
        {
            filterCountsStr.add(format("%s=%d", reason, filterCounts[reason.ordinal()]));
        }

        PPL_LOGGER.debug("variants({}) fitting({}) filters: {}", variants.size(), fittingVariants.size(), filterCountsStr);

        return fittingVariants;
    }

    private static boolean isFittingCandidate(final SomaticVariant variant, final int[] filterCounts)
    {
        if(!variant.hasTumorAlleleDepth() || variant.tumorAlleleDepth().totalReadCount() == 0)
            return false;

        VariantTier variantTier = variant.decorator().tier();

        if(variantTier != VariantTier.HOTSPOT)
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

            if(variant.decorator().mappability() < SNV_FITTING_MAPPABILITY)
            {
                ++filterCounts[FilterReason.MAPPABILITY.ordinal()];
                return false;
            }
        }

        if(variant.referenceAlleleReadCount() > 0)
        {
            ++filterCounts[FilterReason.GERMLINE_ALLELE_COUNT.ordinal()];
            return false;
        }

        return true;
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
            if(!variant.isPass() || !variant.isHotspot())
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

            PPL_LOGGER.info(format("hotspot(%s:%d) vaf(%.3f %d/%d)",
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
