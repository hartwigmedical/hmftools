package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LIST_SEPARATOR;
import static com.hartwig.hmftools.common.variant.SageVcfTags.RC_REALIGNED;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_QUALITY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.MAPPABILITY_TAG;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.common.CommonUtils.DEFAULT_PROBE_LENGTH;
import static com.hartwig.hmftools.wisp.common.CommonUtils.generateMutationSequence;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.MAX_SUBCLONAL_LIKELIHOOD;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SUBCLONAL_VCN_THRESHOLD;
import static com.hartwig.hmftools.wisp.purity.variant.FilterReason.GC_RATIO;
import static com.hartwig.hmftools.wisp.purity.variant.FilterReason.LOW_CONFIDENCE;
import static com.hartwig.hmftools.wisp.purity.variant.FilterReason.LOW_QUAL_PER_AD;
import static com.hartwig.hmftools.wisp.purity.variant.FilterReason.MAPPABILITY;
import static com.hartwig.hmftools.wisp.purity.variant.FilterReason.NON_SNV;
import static com.hartwig.hmftools.wisp.purity.variant.FilterReason.NO_FILTER;
import static com.hartwig.hmftools.wisp.purity.variant.FilterReason.NO_PASS;
import static com.hartwig.hmftools.wisp.purity.variant.FilterReason.REPEAT_COUNT;
import static com.hartwig.hmftools.wisp.purity.variant.FilterReason.SUBCLONAL;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityCalc.LOW_PROBABILITY;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityCalc.calcPoissonNoiseValue;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityCalc.estimatedPurity;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticVariantResult.INVALID_RESULT;
import static com.hartwig.hmftools.wisp.purity.variant.UmiTypeCounts.NO_UMI_COUNTS;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.wisp.purity.PurityConfig;
import com.hartwig.hmftools.wisp.purity.WriteType;
import com.hartwig.hmftools.wisp.purity.ResultsWriter;
import com.hartwig.hmftools.wisp.common.SampleData;
import com.hartwig.hmftools.wisp.purity.PurityConstants;

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
            somaticVcf = mConfig.SampleDataDir + mSample.TumorId + PurityConstants.PURPLE_CTDNA_SOMATIC_VCF_ID;

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
        double sequenceGcRatio = -1;

        if(mConfig.RefGenome != null)
        {
            String variantRefContext = generateMutationSequence(
                    mConfig.RefGenome, DEFAULT_PROBE_LENGTH, variant.chromosome(), variant.position(), variant.ref(), variant.alt());
            sequenceGcRatio = calcGcPercent(variantRefContext);
        }

        List<FilterReason> filterReasons = checkFilters(variant, subclonalLikelihood, sequenceGcRatio);

        SomaticVariant somaticVariant = null;

        for(Genotype genotype : variantContext.getGenotypes())
        {
            if(!targetSampleIds.contains(genotype.getSampleName()))
                continue;

            if(somaticVariant == null)
            {
                somaticVariant = new SomaticVariant(variant, subclonalLikelihood, filterReasons);

                somaticVariant.setSequenceGcRatio(sequenceGcRatio);

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

            if(umiTypeCounts != NO_UMI_COUNTS)
            {
                // override basic AD and DP if umit totals are set
                depth = umiTypeCounts.total();
                alleleCount = umiTypeCounts.alleleTotal();
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
        int frag1Variants = 0;
        int frag2PlusVariants = 0;
        int fragTotal = 0;
        double depthFragTotal = 0;

        UmiTypeCounts umiTypeCounts = new UmiTypeCounts();

        for(SomaticVariant variant : mVariants)
        {
            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);
            GenotypeFragments tumorFragData = variant.findGenotypeData(mSample.TumorId);

            if(sampleFragData == null || tumorFragData == null)
                continue;

            ++totalVariants;

            List<FilterReason> filterReasons = Lists.newArrayList(variant.filterReasons());

            boolean useForTotals = false;

            if(filterReasons.isEmpty())
            {
                if(sampleFragData.isLowQual())
                {
                    filterReasons.add(LOW_QUAL_PER_AD);
                }
                else
                {
                    useForTotals = true;
                }
            }

            if(mConfig.writeType(WriteType.SOMATICS_ALL) || useForTotals)
            {
                // String filter = variant.PassFilters && useForTotals ? "PASS" : (!variant.PassFilters ? "FILTERED" : "NO_FRAGS");
                mResultsWriter.writeVariant(mSample.PatientId, sampleId, variant, sampleFragData, tumorFragData, filterReasons);
            }

            if(!useForTotals)
                continue;

            ++calcVariants;

            if(sampleFragData.AlleleCount >= 2)
                ++frag2PlusVariants;
            else if(sampleFragData.AlleleCount == 1)
                ++frag1Variants;

            fragTotal += sampleFragData.AlleleCount;
            depthFragTotal += sampleFragData.AlleleCount * sampleFragData.Depth;

            // take the tumor values
            tumorCounts.addFragmentCount(tumorFragData.Depth);
            tumorCounts.addAlleleFragmentCount(tumorFragData.AlleleCount, tumorFragData.QualTotal);

            // take the sample values
            sampleCounts.addFragmentCount(sampleFragData.Depth);
            sampleCountsDual.addFragmentCount(sampleFragData.UmiCounts.dualTotal());

            umiTypeCounts.add(sampleFragData.UmiCounts);

            sampleCounts.addAlleleFragmentCount(sampleFragData.AlleleCount, sampleFragData.QualTotal);
            sampleCountsDual.addAlleleFragmentCount(sampleFragData.UmiCounts.AlleleDual, 0);
        }

        if(totalVariants == 0)
            return INVALID_RESULT;

        int sampleDepthTotal = sampleCounts.totalFragments();
        if(sampleDepthTotal == 0)
            return INVALID_RESULT;

        double tumorVaf;

        if(!mConfig.hasSyntheticTumor())
        {
            double tumorDepthTotal = tumorCounts.totalFragments();
            if(tumorDepthTotal == 0)
                return INVALID_RESULT;

            tumorVaf = mConfig.hasSyntheticTumor() ? 0.5 : tumorCounts.alleleFragments() / tumorDepthTotal;
        }
        else
        {
            tumorVaf = 0.5;
        }

        double tumorPurity = purityContext.bestFit().purity();
        double tumorPloidy = purityContext.bestFit().ploidy();

        double adjustedTumorVaf = tumorVaf * (tumorPloidy * tumorPurity + 2 * (1 - tumorPurity)) / tumorPurity / tumorPloidy;

        double qualPerAllele = sampleCounts.alleleFragments() > 0 ? sampleCounts.allelelQualTotal() / (double)sampleCounts.alleleFragments() : 0;

        double lowQualNoiseFactor = qualPerAllele < PurityConstants.LOW_QUAL_NOISE_CUTOFF ?
                1.0 * (PurityConstants.LOW_QUAL_NOISE_CUTOFF - qualPerAllele) / (
                        PurityConstants.LOW_QUAL_NOISE_CUTOFF - PurityConstants.MIN_QUAL_PER_AD) : 0;

        double allFragsNoise = sampleDepthTotal / 1000000.0 * mConfig.NoiseReadsPerMillion + lowQualNoiseFactor;

        FragmentCalcResult allFragsResult = SomaticPurityCalc.calc(
                tumorPloidy, adjustedTumorVaf, sampleDepthTotal, sampleCounts.alleleFragments(), allFragsNoise);

        double rawSamplePurity = allFragsResult.EstimatedPurity;

        double weightedAvgDepth = fragTotal > 0 ? depthFragTotal / fragTotal : sampleDepthTotal / (double)calcVariants;

        ClonalityResult modelResult = ClonalityResult.INVALID_RESULT;
        ClonalityMethod clonalityMethod = ClonalityMethod.NONE;

        if(frag2PlusVariants >= PurityConstants.VAF_PEAK_MODEL_MIN_FRAG_VARIANTS && weightedAvgDepth > PurityConstants.VAF_PEAK_MODEL_MIN_AVG_DEPTH)
        {
            ClonalityModel model = new VafPeakModel(mConfig, mResultsWriter, mSample, mVariants);
            modelResult = model.calculate(sampleId, allFragsResult);

            if(modelResult == ClonalityResult.INVALID_RESULT)
                clonalityMethod = ClonalityMethod.NO_PEAK; // record the attempt
        }
        else if(frag1Variants + frag2PlusVariants >= PurityConstants.LOW_COUNT_MODEL_MIN_FRAG_VARIANTS && weightedAvgDepth < PurityConstants.LOW_COUNT_MODEL_MIN_AVG_DEPTH)
        {
            ClonalityModel model = new LowCountModel(mConfig, mResultsWriter, mSample, mVariants);
            modelResult = model.calculate(sampleId, allFragsResult);
        }

        int clonalityVarCount = calcVariants;
        double clonalityDropout = 0;

        if(modelResult != ClonalityResult.INVALID_RESULT)
        {
            clonalityMethod = modelResult.Method;
            clonalityVarCount = modelResult.VariantCount;
            clonalityDropout = modelResult.DropoutRate;

            double modelPurity = estimatedPurity(modelResult.Vaf, tumorPloidy, adjustedTumorVaf);
            double modelPurityLow = estimatedPurity(modelResult.VafLow, tumorPloidy, adjustedTumorVaf);
            double modelPurityHigh = estimatedPurity(modelResult.VafHigh, tumorPloidy, adjustedTumorVaf);

            // CT_LOGGER.debug(format("sample(%s) tumor(purity=%.4f ploidy=%.1f adjVaf=%.4f) vaf(basic=%.6f peak=%.6f) newPurity(%.6f)",
            //        sampleId, tumorPurity, tumorPloidy, adjustedTumorVaf, allFragsResult.VAF, modelResult.Vaf, modelPurity));

            allFragsResult = new FragmentCalcResult(
                    modelResult.Vaf, modelPurity, allFragsResult.PurityProbability, modelPurityLow, modelPurityHigh);
        }

        double dualFragsNoise = sampleCountsDual.totalFragments() / 1000000.0 * mConfig.NoiseReadsPerMillionDualStrand + lowQualNoiseFactor;

        FragmentCalcResult dualFragsResult = SomaticPurityCalc.calc(
                tumorPloidy, adjustedTumorVaf, sampleCountsDual.totalFragments(), sampleCountsDual.alleleFragments(), dualFragsNoise);

        // calculate a limit-of-detection (LOD), being the number of fragments that would return a 95% confidence of a tumor presence
        double lodFragments = calcPoissonNoiseValue((int)round(allFragsNoise), LOW_PROBABILITY);

        FragmentCalcResult lodFragsResult = SomaticPurityCalc.calc(
                tumorPloidy, adjustedTumorVaf, sampleDepthTotal, (int)round(lodFragments), allFragsNoise);

        // CT_LOGGER.info(format("patient(%s) sample(%s) sampleTotalFrags(%d) noise(%.1f) LOD(%.6f)",
        //        mSample.PatientId, sampleId, sampleDepthTotal, allFragsNoise, lodFragsResult.EstimatedPurity));

        return new SomaticVariantResult(
                true, totalVariants, calcVariants,  frag1Variants, frag2PlusVariants,
                clonalityVarCount, clonalityMethod, clonalityDropout, weightedAvgDepth, sampleCounts, umiTypeCounts,
                tumorVaf, adjustedTumorVaf, rawSamplePurity, allFragsResult, dualFragsResult, lodFragsResult);
    }

    /*
    private boolean useVariantForPurityCalcs(final SomaticVariant variant, final GenotypeFragments sampleFragData)
    {
        return !variant.isFiltered()
                && (sampleFragData.qualPerAlleleFragment() > PurityConstants.MIN_QUAL_PER_AD || sampleFragData.UmiCounts.alleleTotal() == 0);
    }
    */

    private List<FilterReason> checkFilters(final VariantContextDecorator variant, double subclonalLikelihood, double sequenceGcRatio)
    {
        List<FilterReason> filters = Lists.newArrayList();

        if(variant.context().isFiltered())
            filters.add(NO_PASS);

        if(variant.type() != VariantType.SNP)
            filters.add(NON_SNV);

        if(variant.context().hasAttribute(MAPPABILITY_TAG) && variant.mappability() < 1)
            filters.add(MAPPABILITY);

        if(variant.repeatCount() > PurityConstants.MAX_REPEAT_COUNT)
            filters.add(REPEAT_COUNT);

        if(variant.tier() == VariantTier.LOW_CONFIDENCE)
            filters.add(LOW_CONFIDENCE);

        if(subclonalLikelihood > MAX_SUBCLONAL_LIKELIHOOD && variant.variantCopyNumber() < SUBCLONAL_VCN_THRESHOLD)
            filters.add(SUBCLONAL);

        // check GC content
        if(mConfig.GcRatioMin > 0 && sequenceGcRatio >= 0 && sequenceGcRatio < mConfig.GcRatioMin)
            filters.add(GC_RATIO);

        return filters;
    }
}
