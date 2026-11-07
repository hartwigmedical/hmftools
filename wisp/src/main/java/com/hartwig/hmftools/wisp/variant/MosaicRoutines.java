package com.hartwig.hmftools.wisp.variant;

import static java.lang.Math.ceil;
import static java.lang.Math.max;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.wisp.ResultsWriter.formatPurityValue;
import static com.hartwig.hmftools.wisp.WispConstants.OUTLIER_MIN_ALLELE_FRAGS;
import static com.hartwig.hmftools.wisp.WispConstants.OUTLIER_MIN_ALLELE_FRAGS_WITH_DUAL;
import static com.hartwig.hmftools.wisp.WispConstants.OUTLIER_MIN_AVG_VAF_MULTIPLE;
import static com.hartwig.hmftools.wisp.WispConstants.OUTLIER_MIN_SAMPLE_PERC;
import static com.hartwig.hmftools.wisp.WispConstants.OUTLIER_MIN_SAMPLE_RETEST_PERC;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.variant.SomaticPurityCalcs.estimatedPurity;

import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.apache.commons.math3.distribution.BinomialDistribution;

public final class MosaicRoutines
{
    public static boolean onlyMosaicFilters(final SomaticVariant variant, final GenotypeFragments sampleFragData)
    {
        return variant.filterReasons().stream().allMatch(x -> FilterReason.isMosiacType(x))
            && sampleFragData.filterReasons().stream().allMatch(x -> FilterReason.isMosiacType(x));
    }

    public static List<SomaticVariant> findOutlierVariants(
            final String sampleId, final List<SomaticVariant> filteredVariants, double initialSampleTotalAD)
    {
        // check for CHIP or similar outlier variants and remove them from variants used for purity estimates
        // recompute the total sample AD if any outlier variants are found and repeat
        double sampleTotalAD = initialSampleTotalAD;

        // first compute an average VAF to test outliers against
        double vafTotal = 0;

        for(SomaticVariant variant : filteredVariants)
        {
            GenotypeFragments sampleFragData = variant.findSampleFragmentData(sampleId);

            vafTotal += sampleFragData.vaf();
        }

        double averageVaf = vafTotal / (double)filteredVariants.size();

        List<SomaticVariant> allOutlierVariants = Lists.newArrayList();
        List<SomaticVariant> outlierVariants = Lists.newArrayList();
        double minSamplePerc = OUTLIER_MIN_SAMPLE_PERC;

        while(true)
        {
            for(SomaticVariant variant : filteredVariants)
            {
                GenotypeFragments sampleFragData = variant.findSampleFragmentData(sampleId);

                boolean sufficientFrags = sampleFragData.AlleleCount >= OUTLIER_MIN_ALLELE_FRAGS
                        || (sampleFragData.AlleleCount >= OUTLIER_MIN_ALLELE_FRAGS_WITH_DUAL && sampleFragData.UmiCounts.AlleleDual >= 1);

                if(!sufficientFrags)
                    continue;

                if(sampleFragData.vaf() < OUTLIER_MIN_AVG_VAF_MULTIPLE * averageVaf)
                    continue;

                if((sampleFragData.AlleleCount / variant.variantCnFloored()) / sampleTotalAD > minSamplePerc)
                {
                    sampleFragData.markOutlier();
                    outlierVariants.add(variant);
                }
            }

            if(outlierVariants.isEmpty())
                break;

            allOutlierVariants.addAll(outlierVariants);

            for(SomaticVariant variant : outlierVariants)
            {
                CT_LOGGER.debug(format("sample(%s) outlier variant(%s) normalised ad(%.1f) vs sampleTotal(%.1f)",
                        sampleId, variant, variant.findSampleFragmentData(sampleId).AlleleCount / variant.variantCnFloored(), sampleTotalAD));

                filteredVariants.remove(variant);
            }

            outlierVariants.clear();

            // increase threshold and recompute sample allele total less the identified CHIP variants
            minSamplePerc = OUTLIER_MIN_SAMPLE_RETEST_PERC;

            sampleTotalAD = 0;
            for(SomaticVariant variant : filteredVariants)
            {
                GenotypeFragments sampleFragData = variant.findSampleFragmentData(sampleId);
                sampleTotalAD += sampleFragData.AlleleCount / max(variant.VariantCopyNumber, 1);
            }
        }

        return allOutlierVariants;
    }

    public static String formOutlierVariantsInfo(final String sampleId, final List<SomaticVariant> outlierVariants, final double noiseRate)
    {
        StringJoiner sjOutlier = new StringJoiner(ITEM_DELIM);

        for(SomaticVariant outlier : outlierVariants)
        {
            GenotypeFragments sampleFragData = outlier.findSampleFragmentData(sampleId);
            GenotypeFragments tumorFragData = outlier.tumorFragmentData();

            FragmentTotals variantFragTotals = new FragmentTotals();

            variantFragTotals.addVariantData(
                    outlier.CopyNumber, outlier.VariantCopyNumber, tumorFragData.AlleleCount, sampleFragData.AlleleCount,
                    tumorFragData.Depth, sampleFragData.Depth, false);

            double impliedTF = estimatedPurity(variantFragTotals.rawSampleVaf(), noiseRate, variantFragTotals);

            sjOutlier.add(format("%s %d/%d %.2f %s",
                    outlier, sampleFragData.AlleleCount, sampleFragData.Depth, sampleFragData.vaf(), formatPurityValue(impliedTF)));
        }

        return sjOutlier.toString();
    }

    public static MosaicData assessMosaics(
            final String sampleId, final List<SomaticVariant> filteredVariants, int outlierVariantCount,
            int mosaicEligibleVariantCount, double sampleTotalAD, double sampleAfTotal)
    {
        int preOutlierFilteredCount = filteredVariants.size() + outlierVariantCount;
        boolean checkOutlierFrags = preOutlierFilteredCount >= OUTLIER_MIN_AVG_VAF_MULTIPLE;

        // start at the lower percentage as in outlier identification
        double minOutlierSamplePerc = outlierVariantCount == 0 ? OUTLIER_MIN_SAMPLE_PERC : OUTLIER_MIN_SAMPLE_RETEST_PERC;

        double totalExpectedMosaicFrags = 0;
        double totalExpectedDualMosaicFrags = 0;

        for(SomaticVariant variant : filteredVariants)
        {
            GenotypeFragments sampleFragData = variant.findSampleFragmentData(sampleId);

            if(sampleFragData.Depth == 0)
                continue;

            GenotypeFragments referenceFragData = variant.referenceFragmentData();
            int referenceDp = referenceFragData.Depth;

            double dualDepthRatio = sampleFragData.UmiCounts.TotalDual / sampleFragData.Depth;
            double vaf = sampleFragData.vaf();
            int sampleAd = sampleFragData.AlleleCount;
            int sampleDp = sampleFragData.Depth;

            int minFragsAsOutlier = 0;

            if(checkOutlierFrags)
            {
                // ceiling(CNFloored * minOutlierSamplePerc  * (sampleTotalAD - sampleAD / CNFloored) / (1 - minOutlierSamplePerc))
                double cn = variant.variantCnFloored();

                int minOutlierFragsAD = (int)ceil(cn * minOutlierSamplePerc * (sampleTotalAD - sampleAd / cn) / (1 - minOutlierSamplePerc));

                int minOutlierFragsAF = (int)ceil(OUTLIER_MIN_AVG_VAF_MULTIPLE * sampleDp * (sampleAfTotal - vaf) / outlierVariantCount);

                // Then, starting at minFragsToOutlier = max(adMinFrags, afMinFrags), keep testing increasing minFrags until
                // minFragsToOutlier >= 4 || (minFragsToOutlier >= 2 && round(minFragsToOutlier * siteDualRatio) >= 1), and return that
                minFragsAsOutlier = max(minOutlierFragsAD, minOutlierFragsAF);

                while(minFragsAsOutlier < 4 && !(minFragsAsOutlier >= 2 && round(minFragsAsOutlier * dualDepthRatio) >= 1))
                {
                    ++minFragsAsOutlier;
                }
            }

            // test this for each modelled cell generation
            int cellGeneration = 1;

            // start with cell generation 1 and continue until expectedAD < 1e-3, which will depend on site sample depth
            while(true)
            {
                // for cell generation n, the expected AF of a mosaic variant is 2^-n. Expected AD is expected AF * sampleDP
                double expectedAf = pow(2, -cellGeneration);
                double expectedAd = sampleDp * expectedAf;

                if(expectedAd < 0.001)
                    break;

                // the chance that such a variant would have 0 germline AD in the original T/N is (1-expectedAF)^germlineDP
                double zeroGermlineADProb = pow(1 - expectedAf, referenceDp);

                double notOutlierFilteredProb;

                if(checkOutlierFrags)
                {
                    // The chance that such a variant would NOT be outlier filtered can be computed as:
                    // Binom(n=sampleDP-1, expectedAF).cdf(minFragsToOutlier - 2), or 1 if noPossibleOutlier=F.
                    BinomialDistribution binDistribution = new BinomialDistribution(sampleDp - 1, expectedAf);
                    notOutlierFilteredProb = binDistribution.cumulativeProbability(minFragsAsOutlier - 2);
                }
                else
                {
                    notOutlierFilteredProb = 1;
                }

                // the expected number of mosaic frags for that variant, given that cell generation, is:
                //   expectedAD * zeroGermlineADProb * notOutlierFilteredProb.
                //
                // since one prob is based on random sampling of germlineDP and the other is random sampling of ctDNA sampleDP,
                // they are largely orthogonal and can be multiplied in this way. Call this expectedFrags
                double expectedMosaicFags = expectedAd * zeroGermlineADProb * notOutlierFilteredProb;

                // the expected number of dual mosaic frags for that variant is expectedFrags * siteDualRatio, unless dual-filtered
                double expectedDualMosaicFags = !sampleFragData.isDualFiltered() ? expectedMosaicFags * dualDepthRatio : 0;

                totalExpectedMosaicFrags += expectedMosaicFags;
                totalExpectedDualMosaicFrags += expectedDualMosaicFags;

                ++cellGeneration;
            }
        }

        return new MosaicData(mosaicEligibleVariantCount, totalExpectedMosaicFrags, totalExpectedDualMosaicFrags);
    }

}
