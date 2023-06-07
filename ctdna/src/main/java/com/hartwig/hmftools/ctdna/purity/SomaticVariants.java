package com.hartwig.hmftools.ctdna.purity;

import static com.hartwig.hmftools.common.stats.PoissonCalcs.calcPoissonNoiseValue;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LIST_SEPARATOR;
import static com.hartwig.hmftools.common.variant.SageVcfTags.RC_REALIGNED;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_QUALITY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.MAPPABILITY_TAG;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.MAX_REPEAT_COUNT;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.MIN_QUAL_PER_AD;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.MAX_SUBCLONAL_LIKELIHOOD;
import static com.hartwig.hmftools.ctdna.purity.PurityConstants.PURPLE_CTDNA_SOMATIC_VCF_ID;
import static com.hartwig.hmftools.ctdna.purity.SomaticPurityCalc.LOW_PROBABILITY;
import static com.hartwig.hmftools.ctdna.purity.SomaticVariantResult.INVALID_RESULT;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.VcfFileReader;

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

        for(SomaticVariant variant : mVariants)
        {
            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);
            GenotypeFragments tumorFragData = variant.findGenotypeData(mSample.TumorId);

            if(sampleFragData == null || tumorFragData == null)
                continue;

            ++totalVariants;

            boolean useForTotals = variant.PassFilters
                    && (sampleFragData.qualPerAlleleFragment() >= MIN_QUAL_PER_AD || sampleFragData.UmiCounts.alleleTotal() == 0);

            if(mConfig.WriteFilteredSomatics || useForTotals)
            {
                String filter = variant.PassFilters && useForTotals ? "PASS" : (!variant.PassFilters ? "FILTERED" : "NO_FRAGS");
                mResultsWriter.writeVariant(mSample.PatientId, sampleId, variant, sampleFragData, filter);
            }

            if(!useForTotals)
                continue;

            ++calcVariants;

            // take the tumor values
            tumorCounts.VariantDepths.add(tumorFragData.Depth);

            if(tumorFragData.Depth > 0)
                tumorCounts.NonZeroVariantDepths.add(tumorFragData.Depth);

            tumorCounts.AllelelQualTotal += tumorFragData.QualTotal;
            tumorCounts.AlleleFragments += tumorFragData.AlleleCount;

            // take the sample values
            sampleCounts.VariantDepths.add(sampleFragData.UmiCounts.total());
            sampleCountsDual.VariantDepths.add(sampleFragData.UmiCounts.dualTotal());

            if(sampleFragData.Depth > 0)
                sampleCounts.NonZeroVariantDepths.add(sampleFragData.Depth);

            umiTypeCounts.add(sampleFragData.UmiCounts);

            sampleCounts.AllelelQualTotal += sampleFragData.QualTotal;
            sampleCounts.AlleleFragments += sampleFragData.UmiCounts.alleleTotal();
            sampleCountsDual.AlleleFragments += sampleFragData.UmiCounts.AlleleDual;
        }

        if(totalVariants == 0)
            return INVALID_RESULT;

        int sampleDepthTotal = sampleCounts.depthTotal();
        if(sampleDepthTotal == 0)
            return INVALID_RESULT;

        double tumorDepthTotal = tumorCounts.depthTotal();
        if(tumorDepthTotal == 0)
            return INVALID_RESULT;

        double tumorPurity = purityContext.bestFit().purity();
        double tumorPloidy = purityContext.bestFit().ploidy();

        double tumorVaf = tumorCounts.AlleleFragments / tumorDepthTotal;
        double adjustedTumorVaf = tumorVaf * (tumorPloidy * tumorPurity + 2 * (1 - tumorPurity)) / tumorPurity / tumorPloidy;

        // ctDNA_TF = 2 * cfDNA_VAF / [ PLOIDY * ADJ_PRIMARY_VAF + cfDNA_VAF * ( 2 - PLOIDY)]
        // ADJ_PRIMARY_VAF= PRIMARY_VAF * [ PURITY*PLOIDY - 2*(1-PURITY)]/PURITY/PLOIDY

        /*
        double noise = sampleDepthTotal / 1000000.0 * mConfig.NoiseReadsPerMillion;
        double sampleVaf = max(sampleCounts.AlleleFragments - noise, 0) / (double)sampleDepthTotal;
        double samplePurity = 2 * sampleVaf / (tumorPloidy * adjustedTumorVaf + sampleVaf * (2 - tumorPloidy));

        double probability = 1;

        if(sampleCounts.AlleleFragments > noise && noise > 0)
        {
            PoissonDistribution poisson = new PoissonDistribution(noise);
            probability = 1 - poisson.cumulativeProbability(sampleCounts.AlleleFragments - 1);
        }
        */

        FragmentCalcResult allFragsResult = SomaticPurityCalc.calc(
                tumorPloidy, adjustedTumorVaf, sampleCounts.depthTotal(), sampleCounts.AlleleFragments, mConfig.NoiseReadsPerMillion);

        FragmentCalcResult dualFragsResult = SomaticPurityCalc.calc(
                tumorPloidy, adjustedTumorVaf, sampleCountsDual.depthTotal(), sampleCountsDual.AlleleFragments, mConfig.NoiseReadsPerMillionDualStrand);

        // calculate a limit-of-detection (LOD), being the number of fragments that would return a 95% confidence of a tumor presence
        double lodFragments = calcPoissonNoiseValue((int)allFragsResult.Noise, LOW_PROBABILITY);

        FragmentCalcResult lodFragsResult = SomaticPurityCalc.calc(
                tumorPloidy, adjustedTumorVaf, sampleCountsDual.depthTotal(), (int)lodFragments, mConfig.NoiseReadsPerMillion);

        double qualPerAllele = sampleCounts.AlleleFragments > 0 ? sampleCounts.AllelelQualTotal / (double)sampleCounts.AlleleFragments : 0;

        return new SomaticVariantResult(
                true, totalVariants, calcVariants, sampleCounts, umiTypeCounts, qualPerAllele,
                tumorVaf, adjustedTumorVaf, allFragsResult, dualFragsResult, lodFragsResult);

        /*
        return new SomaticVariantResult(
                true, totalVariants, calcVariants, sampleDepthTotal,
                umiTypeCounts[umiIndex++], umiTypeCounts[umiIndex++], umiTypeCounts[umiIndex++],
                sampleCounts.AlleleFragments, umiTypeCounts[umiIndex++], umiTypeCounts[umiIndex++], umiTypeCounts[umiIndex], qualPerAllele,
                sampleCounts.medianDepth(false), sampleCounts.NonZeroVariantDepths.size(), sampleCounts.medianDepth(true),
                tumorVaf, adjustedTumorVaf, allFragsResult.VAF, allFragsResult.EstimatedPurity, allFragsResult.PurityProbability);
        */
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
