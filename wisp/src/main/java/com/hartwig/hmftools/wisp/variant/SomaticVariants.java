package com.hartwig.hmftools.wisp.variant;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;
import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.probabilityToPhredQual;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_RECALIBRATED_BASE_QUAL;
import static com.hartwig.hmftools.common.variant.SageVcfTags.NEARBY_INDEL_FLAG;
import static com.hartwig.hmftools.wisp.WispConstants.MAX_GERMLINE_AF_MAX;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.common.CommonUtils.DEFAULT_PROBE_LENGTH;
import static com.hartwig.hmftools.wisp.common.CommonUtils.generateMutationSequence;
import static com.hartwig.hmftools.wisp.WispConfig.isUltima;
import static com.hartwig.hmftools.wisp.WispConstants.MIN_AVG_EDGE_DISTANCE;
import static com.hartwig.hmftools.wisp.WispConstants.MIN_AVG_EDGE_DISTANCE_DUAL_ULTIMA;
import static com.hartwig.hmftools.wisp.WispConstants.MIN_DUAL_QUAL_FILTER_THRESHOLD;
import static com.hartwig.hmftools.wisp.WispConstants.MAX_SUBCLONAL_LIKELIHOOD;
import static com.hartwig.hmftools.wisp.WispConstants.SUBCLONAL_VCN_THRESHOLD;
import static com.hartwig.hmftools.wisp.WispConstants.HIGH_GERMLINE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.wisp.WispConstants.MAX_GERMLINE_AF;
import static com.hartwig.hmftools.wisp.WriteType.FRAG_LENGTHS;
import static com.hartwig.hmftools.wisp.variant.FilterReason.AVG_EDGE_DIST;
import static com.hartwig.hmftools.wisp.variant.FilterReason.DUAL_ERROR_RATE;
import static com.hartwig.hmftools.wisp.variant.FilterReason.GC_RATIO;
import static com.hartwig.hmftools.wisp.variant.FilterReason.LOW_CONFIDENCE;
import static com.hartwig.hmftools.wisp.variant.FilterReason.LOW_QUAL_PER_AD;
import static com.hartwig.hmftools.wisp.variant.FilterReason.NON_SNV;
import static com.hartwig.hmftools.wisp.variant.FilterReason.NO_PASS;
import static com.hartwig.hmftools.wisp.variant.FilterReason.GERMLINE_AF;
import static com.hartwig.hmftools.wisp.variant.FilterReason.NEARBY_INDEL;
import static com.hartwig.hmftools.wisp.variant.FilterReason.REPEAT_COUNT;
import static com.hartwig.hmftools.wisp.variant.FilterReason.SUBCLONAL;
import static com.hartwig.hmftools.wisp.variant.GenotypeFragments.AVG_BASE_QUAL_OLD_TAG;
import static com.hartwig.hmftools.wisp.variant.GenotypeFragments.REQUIRED_GENOTYPE_FIELDS;
import static com.hartwig.hmftools.wisp.variant.GenotypeFragments.buildFragments;
import static com.hartwig.hmftools.wisp.variant.MosaicRoutines.assessMosaics;
import static com.hartwig.hmftools.wisp.variant.MosaicRoutines.onlyMosaicFilters;
import static com.hartwig.hmftools.wisp.variant.SnvFitStatus.UNCERTAIN;
import static com.hartwig.hmftools.wisp.variant.SomaticPurityResult.INVALID_RESULT;

import java.io.BufferedWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.PaveVcfTags;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.wisp.WispConfig;
import com.hartwig.hmftools.wisp.WriteType;
import com.hartwig.hmftools.wisp.ResultsWriter;
import com.hartwig.hmftools.wisp.SampleData;
import com.hartwig.hmftools.wisp.WispConstants;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticVariants
{
    private final WispConfig mConfig;
    private final ResultsWriter mResultsWriter;

    private final SampleData mSample;
    private final List<SomaticVariant> mVariants;
    private final List<SimpleVariant> mProbeVariants;
    private final SomaticPurityEstimator mEstimator;
    private final BqrAdjustment mBqrAdjustment;
    private final BufferedWriter mSomaticWriter;
    private final SampleFragmentLengths mFragmentLengths;

    public SomaticVariants(final WispConfig config, final ResultsWriter resultsWriter, final SampleData sampleData)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;
        mSample = sampleData;
        mBqrAdjustment = new BqrAdjustment(mConfig);
        mEstimator = new SomaticPurityEstimator(mConfig, resultsWriter, sampleData, mBqrAdjustment);

        mSomaticWriter = mResultsWriter.getSomaticWriter();
        mFragmentLengths = new SampleFragmentLengths(config, resultsWriter, sampleData);

        mVariants = Lists.newArrayList();

        mProbeVariants = mConfig.ProbeVariants.getSampleVariants(mSample.TumorId);
    }

    public boolean loadVariants()
    {
        String somaticVcf;

        if(mConfig.SomaticVcf.isEmpty())
        {
            String vcfFilename = mConfig.SomaticDir + mSample.TumorId + WispConstants.PURPLE_APPENDED_SOMATIC_VCF_ID;

            if(!mSample.VcfTag.isEmpty())
                vcfFilename += mSample.VcfTag + ".";

            vcfFilename += "vcf.gz";

            if(Files.exists(Paths.get(vcfFilename)))
                somaticVcf = vcfFilename;
            else
                somaticVcf = mConfig.SomaticDir + mSample.TumorId + WispConstants.PURPLE_APPENDED_SOMATIC_VCF_ID + "vcf.gz";
        }
        else
        {
            somaticVcf = mConfig.getSomaticVcf(mSample.TumorId);

            if(!Files.exists(Paths.get(somaticVcf)))
                somaticVcf = mConfig.SampleDataDir + somaticVcf;
        }

        CT_LOGGER.debug("loading somatic variant VCF: {}", somaticVcf);

        List<String> targetSampleIds = Lists.newArrayList(mSample.TumorId);

        if(mSample.hasReference())
            targetSampleIds.add(mSample.ReferenceId);

        mSample.SampleIds.forEach(x -> targetSampleIds.add(x));

        VcfFileReader vcfFileReader = new VcfFileReader(somaticVcf);

        if(!vcfFileReader.fileValid())
        {
            CT_LOGGER.error("failed to read somatic vcf({})", somaticVcf);
            return false;
        }

        VCFHeader vcfHeader = vcfFileReader.vcfHeader();

        // check header has required fields
        if(REQUIRED_GENOTYPE_FIELDS.stream().anyMatch(x -> !vcfHeader.hasFormatLine(x)))
        {
            String requiredVcfTags = REQUIRED_GENOTYPE_FIELDS.stream().collect(Collectors.joining(";"));
            CT_LOGGER.error("patient({}) vcf({}) missing required fields({})",
                    mSample.PatientId, somaticVcf, requiredVcfTags);
            System.exit(1);
        }

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
                CT_LOGGER.info("loaded {} variants", variantCount);
            }
        }

        int matchedProbeCount = 0;
        if(mProbeVariants != null)
        {
            for(SomaticVariant variant : mVariants)
            {
                if(mProbeVariants.stream().anyMatch(x -> variant.matches(x)))
                {
                    ++matchedProbeCount;
                    variant.markProbeVariant();
                }
            }
        }

        CT_LOGGER.info("processed {} somatic variants from VCF({}), filtered({}) probeMatched({})",
                variantCount, filenamePart(somaticVcf), filteredCount, matchedProbeCount);

        if(mConfig.writeType(FRAG_LENGTHS))
            mFragmentLengths.processSample(somaticVcf, mVariants);

        return true;
    }

    private static final double NO_GC_RATIO = -1;

    protected void processVariant(final List<String> targetSampleIds, final VariantContext variantContext)
    {
        VariantContextDecorator variant = new VariantContextDecorator(variantContext);

        // check if the variant has been excluded in config
        if(mConfig.ExcludedSomatics.stream().anyMatch(x -> x.matches(variant)))
            return;

        double subclonalLikelihood = variant.context().getAttributeAsDouble(SUBCLONAL_LIKELIHOOD_FLAG, 0);
        boolean hasSyntheticTumor = mConfig.hasSyntheticTumor();
        double sequenceGcRatio = NO_GC_RATIO;

        if(mConfig.RefGenome != null && mSample.IsPanel)
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
                somaticVariant = new SomaticVariant(variant, subclonalLikelihood, filterReasons, hasSyntheticTumor);

                somaticVariant.setSequenceGcRatio(sequenceGcRatio);
                mVariants.add(somaticVariant);
            }

            if(genotype == null || genotype.getExtendedAttributes().isEmpty())
                continue;

            GenotypeFragments genotypeFragments = buildFragments(
                    genotype, mConfig.DisableDualFragments, genotype.getSampleName().equals(mSample.TumorId));

            if(genotype.getSampleName().equals(mSample.TumorId))
            {
                somaticVariant.setTumorFragmentData(genotypeFragments);
            }
            else if(genotype.getSampleName().equals(mSample.ReferenceId))
            {
                somaticVariant.setReferenceFragmentData(genotypeFragments);
            }
            else
            {
                somaticVariant.addSampleFragmentData(genotypeFragments);
            }
        }
    }

    public SomaticPurityResult processSample(final String sampleId)
    {
        // BQR values are used both in filters and purity calcs
        if(!mConfig.SkipBqr && mConfig.OutputDir != null)
            mBqrAdjustment.loadBqrData(sampleId);

        List<SomaticVariant> filteredVariants = Lists.newArrayList();

        double sampleTotalAD = 0; // these values will be normalised by copy number
        double sampleTotalAF = 0;

        int mosaicEligibleVariantCount = 0;

        // compile a set of filtered variants, and compute a total AD just for these
        for(SomaticVariant variant : mVariants)
        {
            GenotypeFragments sampleFragData = variant.findSampleFragmentData(sampleId);
            GenotypeFragments tumorFragData = variant.tumorFragmentData();

            if(sampleFragData == null || tumorFragData == null)
                continue;

            checkSampleDataFilters(variant, sampleFragData);

            if(onlyMosaicFilters(variant, sampleFragData))
                ++mosaicEligibleVariantCount;

            if(variant.isFiltered() || sampleFragData.isFiltered())
                continue;

            filteredVariants.add(variant);

            sampleTotalAD += sampleFragData.AlleleCount / variant.variantCnFloored();

            sampleTotalAF += sampleFragData.vaf();
        }

        SomaticPurityResult purityResult = INVALID_RESULT;

        if(!filteredVariants.isEmpty())
        {
            // check for outlier / CHIP variants and remove them from variants used for purity estimates
            List<SomaticVariant> outlierVariants = MosaicRoutines.findOutlierVariants(sampleId, filteredVariants, sampleTotalAD);

            MosaicData mosaicData = assessMosaics(
                    sampleId, filteredVariants, outlierVariants.size(), mosaicEligibleVariantCount, sampleTotalAD, sampleTotalAF);

            purityResult = mEstimator.calculatePurity(sampleId, filteredVariants, mVariants.size(), mosaicData, true);

            if(purityResult.Status == UNCERTAIN || !purityResult.mrdDetected())
            {
                SomaticPurityResult newPurityResult = mEstimator.calculatePurity(
                        sampleId, filteredVariants, mVariants.size(), mosaicData, false);

                if(newPurityResult.mrdDetected())
                {
                    CT_LOGGER.warn("sample({}) recovered MRD detection using high copy number variants, status({})",
                            sampleId, purityResult.Status);
                    purityResult = newPurityResult;
                }
            }

            String outlierInfo = MosaicRoutines.formOutlierVariantsInfo(sampleId, outlierVariants, purityResult.PurityCalcs.ErrorRate);
            purityResult.setOutlierVariantInfo(outlierInfo);
        }

        if(mConfig.writeType(WriteType.SOMATIC_DATA))
        {
            for(SomaticVariant variant : mVariants)
            {
                GenotypeFragments sampleFragData = variant.findSampleFragmentData(sampleId);
                GenotypeFragments tumorFragData = variant.tumorFragmentData();

                if(sampleFragData != null && tumorFragData != null)
                {
                    SomaticWriter.writeVariant(mSomaticWriter, mConfig, mSample, sampleId, variant, sampleFragData, tumorFragData);
                }
            }
        }

        return purityResult;
    }

    @VisibleForTesting
    protected List<FilterReason> checkFilters(final VariantContextDecorator variant, double subclonalLikelihood, double sequenceGcRatio)
    {
        List<FilterReason> filters = Lists.newArrayList();

        if(mSample.hasReference())
        {
            Genotype refGenotype = variant.context().getGenotype(mSample.ReferenceId);

            AllelicDepth refAllelicDepth = AllelicDepth.fromGenotype(refGenotype);
            double germlineAF = refAllelicDepth.alleleFrequency();

            // revert to previous tag if only that is present - this is permitted for the germline which can be run on v4.2+ of Sage
            String germlineABQStr = refGenotype.getExtendedAttribute(AVG_RECALIBRATED_BASE_QUAL, "").toString();
            String germlineABQOldStr = refGenotype.getExtendedAttribute(AVG_BASE_QUAL_OLD_TAG, "").toString();

            if(germlineABQStr.isEmpty() && !germlineABQOldStr.isEmpty())
                germlineABQStr = germlineABQOldStr;

            double germlineABQ = Double.parseDouble(germlineABQStr.toString().split(CSV_DELIM)[1]);

            if((germlineAF >= MAX_GERMLINE_AF && germlineABQ >= HIGH_GERMLINE_QUAL_THRESHOLD) || germlineAF >= MAX_GERMLINE_AF_MAX)
            {
                filters.add(GERMLINE_AF);
            }
        }

        if(variant.context().getCommonInfo().hasAttribute(NEARBY_INDEL_FLAG))
            filters.add(NEARBY_INDEL);

        if(variant.context().isFiltered())
            filters.add(NO_PASS);

        int maxRepeatCount = max(variant.repeatCount(), variant.altRepeatCount());
        if(maxRepeatCount > WispConstants.MAX_REPEAT_COUNT)
            filters.add(REPEAT_COUNT);

        if(mConfig.ApplyRefVariantFilters) // no others are applied
            return filters;

        if(variant.type() != VariantType.SNP)
            filters.add(NON_SNV);

        if(variant.context().hasAttribute(PaveVcfTags.MAPPABILITY) && variant.mappability() < 1)
            filters.add(FilterReason.MAPPABILITY);

        if(variant.tier() == VariantTier.LOW_CONFIDENCE)
            filters.add(LOW_CONFIDENCE);

        if(!mConfig.SkipSubclonalFilter)
        {
            if(subclonalLikelihood > MAX_SUBCLONAL_LIKELIHOOD && variant.variantCopyNumber() < SUBCLONAL_VCN_THRESHOLD)
                filters.add(SUBCLONAL);
        }

        // check GC content
        if(sequenceGcRatio != NO_GC_RATIO && mConfig.GcRatioMin > 0 && sequenceGcRatio < mConfig.GcRatioMin)
            filters.add(GC_RATIO);

        return filters;
    }

    public void checkSampleDataFilters(final SomaticVariant variant, final GenotypeFragments sampleFragData)
    {
        if(sampleFragData.AlleleCount > 0)
        {
            if(sampleFragData.qualPerAlleleFragment() <= WispConstants.MIN_QUAL_PER_AD)
                sampleFragData.addFilterReason(LOW_QUAL_PER_AD);

            if(isUltima() && sampleFragData.qualPerAlleleFragment() <= WispConstants.MIN_QUAL_PER_AD_DUAL_ULTIMA)
                sampleFragData.addDualFilterReason(LOW_QUAL_PER_AD);

            if(sampleFragData.averageReadDistance() >= 0)
            {
                if(sampleFragData.averageReadDistance() < MIN_AVG_EDGE_DISTANCE)
                    sampleFragData.addFilterReason(AVG_EDGE_DIST);

                if(isUltima() && sampleFragData.averageReadDistance() < MIN_AVG_EDGE_DISTANCE_DUAL_ULTIMA)
                    sampleFragData.addDualFilterReason(AVG_EDGE_DIST);
            }

            if(variant.Type == VariantType.SNP && !mConfig.SkipBqr && mBqrAdjustment.hasDualEntries())
            {
                if(!mBqrAdjustment.hasErrorRate(variant.TriNucContext, variant.Alt, ConsensusType.DUAL))
                {
                    CT_LOGGER.error("variant({}) missing dual BQR data", variant);
                    System.exit(1);
                }

                double dualBqrErrorRate = mEstimator.getBqrErrorRate(variant, sampleFragData, ConsensusType.DUAL);
                sampleFragData.setDualBqrErrorRate(dualBqrErrorRate);

                double dualBqrErrorQual = probabilityToPhredQual(dualBqrErrorRate);

                if(dualBqrErrorQual < MIN_DUAL_QUAL_FILTER_THRESHOLD)
                    sampleFragData.addDualFilterReason(DUAL_ERROR_RATE);
            }
        }
    }

    public BqrAdjustment bqrAdjustment() { return mBqrAdjustment; }
}
