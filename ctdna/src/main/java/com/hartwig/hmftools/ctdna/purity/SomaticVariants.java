package com.hartwig.hmftools.ctdna.purity;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LIST_SEPARATOR;
import static com.hartwig.hmftools.common.variant.SageVcfTags.RC_REALIGNED;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_QUALITY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.MAPPABILITY_TAG;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.LOW_QUAL_NOISE_CUTOFF;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.MAX_REPEAT_COUNT;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.MIN_QUAL_PER_AD;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.MAX_SUBCLONAL_LIKELIHOOD;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.PURPLE_CTDNA_SOMATIC_VCF_ID;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.SAMPLE_ALLELE_OUTLIER_PROBABIILTY;
import static com.hartwig.hmftools.ctdna.purity.SomaticPurityCalc.LOW_PROBABILITY;
import static com.hartwig.hmftools.ctdna.purity.SomaticPurityCalc.calcPoissonNoiseValue;
import static com.hartwig.hmftools.ctdna.purity.SomaticVariantResult.INVALID_RESULT;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import org.apache.commons.math3.distribution.PoissonDistribution;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticVariants
{
    private final PurityConfig mConfig;
    private final ResultsWriter mResultsWriter;

    private final SampleData mSample;
    private final List<SomaticVariant> mVariants;

    public SomaticVariants(final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;
        mSample = sample;

        mVariants = Lists.newArrayList();
    }

    public boolean loadVariants()
    {
        String somaticVcf = mConfig.SomaticVcf;

        if(somaticVcf.isEmpty())
        {
            somaticVcf = mConfig.SampleDataDir + mSample.TumorId + PURPLE_CTDNA_SOMATIC_VCF_ID;

            if(!mSample.VcfTag.isEmpty())
                somaticVcf += mSample.VcfTag + ".";

            somaticVcf += "vcf.gz";
        }
        else
        {
            if(!Files.exists(Paths.get(somaticVcf)))
                somaticVcf = mConfig.SampleDataDir + somaticVcf;
        }

        CT_LOGGER.debug("loading somatic variant VCF: {}", somaticVcf);

        List<String> targetSampleIds = Lists.newArrayList(mSample.TumorId);
        mSample.CtDnaSamples.forEach(x -> targetSampleIds.add(x));

        VcfFileReader vcfFileReader = new VcfFileReader(somaticVcf);

        if(!vcfFileReader.fileValid())
        {
            CT_LOGGER.error("failed to read somatic vcf({})", somaticVcf);
            return false;
        }

        VCFHeader vcfHeader = vcfFileReader.vcfHeader();

        for(String sampleId : targetSampleIds)
        {
            if(!vcfHeader.getGenotypeSamples().contains(sampleId))
            {
                CT_LOGGER.error("patient({}) missing sample({}) in vcf({})", mSample.PatientId, sampleId, somaticVcf);
                return false;
            }
        }

        int filteredCount = 0;
        int variantCount = 0;

        for(VariantContext variantContext : vcfFileReader.iterator())
        {
            if(variantContext.isFiltered())
            {
                ++filteredCount;
                continue;
            }

            ++variantCount;

            try
            {
                processVariant(targetSampleIds, variantContext);
            }
            catch(Exception e)
            {
                e.printStackTrace();
                CT_LOGGER.error("error processing VCF({}): {}", somaticVcf, e.toString());
                return false;
            }

            if(variantCount > 0 && (variantCount % 100000) == 0)
            {
                CT_LOGGER.info("processed {} variants", variantCount);
            }
        }

        CT_LOGGER.info("processed {} filtered({}) somatic variants from VCF({})", variantCount, filteredCount, somaticVcf);

        return true;
    }

    private void processVariant(final List<String> targetSampleIds, final VariantContext variantContext)
    {
        VariantContextDecorator variant = new VariantContextDecorator(variantContext);

        double subclonalLikelihood = variant.context().getAttributeAsDouble(SUBCLONAL_LIKELIHOOD_FLAG, 0);

        boolean isFiltered = filterVariant(variant, subclonalLikelihood);

        SomaticVariant somaticVariant = null;

        for(Genotype genotype : variantContext.getGenotypes())
        {
            if(!targetSampleIds.contains(genotype.getSampleName()))
                continue;

            if(somaticVariant == null)
            {
                somaticVariant = new SomaticVariant(
                        variant.chromosome(), variant.position(), variant.ref(), variant.alt(), variant.tier(),
                        variant.type(), variant.repeatCount(), variant.mappability(), subclonalLikelihood, !isFiltered);
                mVariants.add(somaticVariant);
            }

            if(genotype == null || genotype.getExtendedAttributes().isEmpty())
                continue;

            int depth = genotype.getDP();
            int alleleCount = genotype.getAD()[1];
            UmiTypeCounts umiTypeCounts = UmiTypeCounts.fromAttribute(genotype.getExtendedAttribute(UMI_TYPE_COUNTS, null));

            int qualTotal = 0;

            if(alleleCount > 0)
            {
                final String[] qualCounts = genotype.getExtendedAttribute(READ_CONTEXT_QUALITY, 0).toString()
                        .split(LIST_SEPARATOR, -1);

                for(int i = 0; i <= RC_REALIGNED; ++i)
                {
                    qualTotal += Integer.parseInt(qualCounts[i]);
                }
            }

            somaticVariant.Samples.add(new GenotypeFragments(genotype.getSampleName(), alleleCount, depth, qualTotal, umiTypeCounts));
        }
    }

    public SomaticVariantResult processSample(final String sampleId, final PurityContext purityContext)
    {
        // only include variants which satisfy the min avg qual check in the ctDNA sample
        SomaticVariantCounts tumorCounts = new SomaticVariantCounts();
        SomaticVariantCounts sampleCounts = new SomaticVariantCounts();
        SomaticVariantCounts sampleCountsDual = new SomaticVariantCounts();

        int totalVariants = 0;
        int calcVariants = 0;

        UmiTypeCounts umiTypeCounts = new UmiTypeCounts();

        double sampleTumorAlleleRatio = calcTumorSampleAlleleRatio(sampleId);
        double sampleTumorVafRatio = calcTumorSampleVafRatio(sampleId);

        for(SomaticVariant variant : mVariants)
        {
            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);
            GenotypeFragments tumorFragData = variant.findGenotypeData(mSample.TumorId);

            if(sampleFragData == null || tumorFragData == null)
                continue;

            ++totalVariants;

            boolean useForTotals = useVariantForPurityCalcs(variant, sampleFragData);

            if(useForTotals)
            {
                // useForTotals = !variantIsVafOutlier(sampleTumorVafRatio, variant, tumorFragData, sampleFragData);

                /*
                useForTotals = !variantIsVafOutlier(
                        sampleTumorAlleleRatio, tumorFragData.AlleleCount, sampleFragData.UmiCounts.alleleTotal());
                 */
            }

            if(mConfig.WriteFilteredSomatics || useForTotals)
            {
                String filter = variant.PassFilters && useForTotals ? "PASS" : (!variant.PassFilters ? "FILTERED" : "NO_FRAGS");
                mResultsWriter.writeVariant(mSample.PatientId, sampleId, variant, sampleFragData, filter);
            }

            if(!useForTotals)
                continue;

            ++calcVariants;

            // take the tumor values
            tumorCounts.addFragmentCount(tumorFragData.Depth);
            tumorCounts.addAlleleFragmentCount(tumorFragData.AlleleCount, tumorFragData.QualTotal);

            // take the sample values
            sampleCounts.addFragmentCount(sampleFragData.UmiCounts.total());
            sampleCountsDual.addFragmentCount(sampleFragData.UmiCounts.dualTotal());

            umiTypeCounts.add(sampleFragData.UmiCounts);

            sampleCounts.addAlleleFragmentCount(sampleFragData.UmiCounts.alleleTotal(), sampleFragData.QualTotal);
            sampleCountsDual.addAlleleFragmentCount(sampleFragData.UmiCounts.AlleleDual, 0);
        }

        if(totalVariants == 0)
            return INVALID_RESULT;

        int sampleDepthTotal = sampleCounts.totalFragments();
        if(sampleDepthTotal == 0)
            return INVALID_RESULT;

        double tumorDepthTotal = tumorCounts.totalFragments();
        if(tumorDepthTotal == 0)
            return INVALID_RESULT;

        double tumorPurity = purityContext.bestFit().purity();
        double tumorPloidy = purityContext.bestFit().ploidy();

        double tumorVaf = tumorCounts.alleleFragments() / tumorDepthTotal;
        double adjustedTumorVaf = tumorVaf * (tumorPloidy * tumorPurity + 2 * (1 - tumorPurity)) / tumorPurity / tumorPloidy;

        double qualPerAllele = sampleCounts.alleleFragments() > 0 ? sampleCounts.allelelQualTotal() / (double)sampleCounts.alleleFragments() : 0;

        double lowQualNoiseFactor = qualPerAllele < LOW_QUAL_NOISE_CUTOFF ?
                1.0 * (LOW_QUAL_NOISE_CUTOFF - qualPerAllele) / (LOW_QUAL_NOISE_CUTOFF - MIN_QUAL_PER_AD) : 0;

        double allFragsNoise = sampleDepthTotal / 1000000.0 * mConfig.NoiseReadsPerMillion + lowQualNoiseFactor;

        FragmentCalcResult allFragsResult = SomaticPurityCalc.calc(
                tumorPloidy, adjustedTumorVaf, sampleDepthTotal, sampleCounts.alleleFragments(), allFragsNoise);

        double dualFragsNoise = sampleCountsDual.totalFragments() / 1000000.0 * mConfig.NoiseReadsPerMillionDualStrand + lowQualNoiseFactor;

        FragmentCalcResult dualFragsResult = SomaticPurityCalc.calc(
                tumorPloidy, adjustedTumorVaf, sampleCountsDual.totalFragments(), sampleCountsDual.alleleFragments(), dualFragsNoise);

        // calculate a limit-of-detection (LOD), being the number of fragments that would return a 95% confidence of a tumor presence
        double lodFragments = calcPoissonNoiseValue((int)round(allFragsNoise), LOW_PROBABILITY);

        FragmentCalcResult lodFragsResult = SomaticPurityCalc.calc(
                tumorPloidy, adjustedTumorVaf, sampleDepthTotal, (int)round(lodFragments), allFragsNoise);

        return new SomaticVariantResult(
                true, totalVariants, calcVariants, sampleCounts, umiTypeCounts, qualPerAllele,
                tumorVaf, adjustedTumorVaf, allFragsResult, dualFragsResult, lodFragsResult);
    }

    private double calcTumorSampleVafRatio(final String sampleId)
    {
        double tumorVafTotal = 0;
        double sampleVafTotal = 0;

        for(SomaticVariant variant : mVariants)
        {
            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);
            GenotypeFragments tumorFragData = variant.findGenotypeData(mSample.TumorId);

            if(sampleFragData == null || tumorFragData == null)
                continue;

            if(useVariantForPurityCalcs(variant, sampleFragData))
            {
                if(tumorFragData.Depth == 0 || sampleFragData.UmiCounts.total() == 0)
                    continue;

                tumorVafTotal += tumorFragData.AlleleCount / (double)tumorFragData.Depth;
                sampleVafTotal += sampleFragData.AlleleCount / (double)sampleFragData.UmiCounts.total();
            }
        }

        return tumorVafTotal > 0 ? sampleVafTotal / tumorVafTotal : 0;
    }

    private double calcTumorSampleAlleleRatio(final String sampleId)
    {
        int tumorFragments = 0;
        int sampleFragments = 0;

        for(SomaticVariant variant : mVariants)
        {
            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);
            GenotypeFragments tumorFragData = variant.findGenotypeData(mSample.TumorId);

            if(sampleFragData == null || tumorFragData == null)
                continue;

            if(useVariantForPurityCalcs(variant, sampleFragData))
            {
                // tumorFragments += tumorFragData.AlleleCount;
                // sampleFragments += sampleFragData.AlleleCount;

                if(tumorFragData.Depth == 0 || sampleFragData.UmiCounts.total() == 0)
                    continue;

                double tumorVaf = tumorFragData.AlleleCount / (double)tumorFragData.Depth;;
                double sampleVaf = sampleFragData.AlleleCount / (double)sampleFragData.UmiCounts.total();

                tumorFragments += tumorVaf;
                sampleFragments += sampleVaf;
            }
        }

        return tumorFragments > 0 ? sampleFragments / (double)tumorFragments : 0;
    }
    private boolean variantIsVafOutlier(
            double sampleTumorVafRatio, final SomaticVariant variant,
            final GenotypeFragments tumorFragData, final GenotypeFragments sampleFragData)
    {
        int sampleAlleleCount = sampleFragData.AlleleCount;

        if(sampleAlleleCount == 1)
            return false;

        double tumorVaf = tumorFragData.AlleleCount / (double)tumorFragData.Depth;
        double expectedVaf = tumorVaf / sampleTumorVafRatio;
        double expectedCount = expectedVaf * sampleFragData.UmiCounts.total();

        // poison ( AD cfDNA,sum(AD cfDNA) / sum(AD tumor)*AD tumor, true)
        // double expectedCount = sampleTumorVafRatio * tumorAlleleCount;

        if(sampleAlleleCount <= expectedCount)
            return false;

        if(expectedCount == 0)
            return sampleAlleleCount > 10;

        PoissonDistribution poisson = new PoissonDistribution(expectedCount);
        double probability = 1 - poisson.cumulativeProbability(sampleAlleleCount - 1);

        if(probability > SAMPLE_ALLELE_OUTLIER_PROBABIILTY)
            return false;

        CT_LOGGER.info(format("sample(%s) variant(%s) tumor(%d vaf=%.3f) sample(%d/%d expVaf=%.3f expCount=%.0f) ratio(%.3g) is outlier",
                sampleFragData.SampleName, variant.toString(), tumorFragData.AlleleCount, tumorVaf,
                sampleFragData.UmiCounts.alleleTotal(), sampleFragData.UmiCounts.total(),
                expectedVaf, expectedCount, sampleTumorVafRatio));

        return true;
    }

    private boolean variantIsVafOutlier(double sampleTumorAlleleRatio, int tumorAlleleCount, int sampleAlleleCount)
    {
        if(sampleAlleleCount == 1)
            return false;

        // poison ( AD cfDNA,sum(AD cfDNA) / sum(AD tumor)*AD tumor, true)
        double expectedCount = sampleTumorAlleleRatio * tumorAlleleCount;

        if(sampleAlleleCount <= expectedCount)
            return false;

        PoissonDistribution poisson = new PoissonDistribution(expectedCount);
        double probability = 1 - poisson.cumulativeProbability(sampleAlleleCount - 1);
        return probability < SAMPLE_ALLELE_OUTLIER_PROBABIILTY;
    }

    private boolean useVariantForPurityCalcs(final SomaticVariant variant, final GenotypeFragments sampleFragData)
    {
        return variant.PassFilters
                && (sampleFragData.qualPerAlleleFragment() > MIN_QUAL_PER_AD || sampleFragData.UmiCounts.alleleTotal() == 0);
    }

    private boolean filterVariant(final VariantContextDecorator variant, double subclonalLikelihood)
    {
        if(variant.context().isFiltered())
            return true;

        if(variant.type() != VariantType.SNP)
            return true;

        if(variant.context().hasAttribute(MAPPABILITY_TAG) && variant.mappability() < 1)
            return true;

        if(variant.repeatCount() > MAX_REPEAT_COUNT)
            return true;

        if(variant.tier() == VariantTier.LOW_CONFIDENCE)
            return true;

        if(subclonalLikelihood > MAX_SUBCLONAL_LIKELIHOOD)
            return true;

        return false;
    }

}
