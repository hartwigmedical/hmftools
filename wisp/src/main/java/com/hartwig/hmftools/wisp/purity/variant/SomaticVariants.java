package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LIST_SEPARATOR;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_QUALITY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.MAPPABILITY_TAG;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.common.CommonUtils.DEFAULT_PROBE_LENGTH;
import static com.hartwig.hmftools.wisp.common.CommonUtils.generateMutationSequence;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.CHIP_MIN_ALLELE_FRAGS;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.CHIP_MIN_SAMPLE_PERC;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.MAX_SUBCLONAL_LIKELIHOOD;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SUBCLONAL_VCN_THRESHOLD;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.SOMATICS_FILE_ID;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.addCommonFields;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.addCommonHeaderFields;
import static com.hartwig.hmftools.wisp.purity.variant.FilterReason.CHIP;
import static com.hartwig.hmftools.wisp.purity.variant.FilterReason.GC_RATIO;
import static com.hartwig.hmftools.wisp.purity.variant.FilterReason.LOW_CONFIDENCE;
import static com.hartwig.hmftools.wisp.purity.variant.FilterReason.LOW_QUAL_PER_AD;
import static com.hartwig.hmftools.wisp.purity.variant.FilterReason.MAPPABILITY;
import static com.hartwig.hmftools.wisp.purity.variant.FilterReason.NON_SNV;
import static com.hartwig.hmftools.wisp.purity.variant.FilterReason.NO_PASS;
import static com.hartwig.hmftools.wisp.purity.variant.FilterReason.REPEAT_COUNT;
import static com.hartwig.hmftools.wisp.purity.variant.FilterReason.SUBCLONAL;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityResult.INVALID_RESULT;
import static com.hartwig.hmftools.wisp.purity.variant.UmiTypeCounts.NO_UMI_COUNTS;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.utils.r.RExecutor;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantReadSupport;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.wisp.purity.PurityConfig;
import com.hartwig.hmftools.wisp.purity.WriteType;
import com.hartwig.hmftools.wisp.purity.ResultsWriter;
import com.hartwig.hmftools.wisp.purity.SampleData;
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
    private final List<ProbeVariant> mProbeVariants;
    private final SomaticPurityEstimator mEstimator;
    private final BufferedWriter mSomaticWriter;

    public SomaticVariants(final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;
        mSample = sample;
        mEstimator = new SomaticPurityEstimator(mConfig, resultsWriter, sample);
        mSomaticWriter = mResultsWriter.getSomaticWriter();

        mVariants = Lists.newArrayList();

        mProbeVariants = mConfig.ProbeVariants.getSampleVariants(mSample.TumorId);
    }

    public boolean loadVariants()
    {
        String somaticVcf;

        if(mConfig.SomaticVcf.isEmpty())
        {
            somaticVcf = mConfig.SomaticDir + mSample.TumorId + PurityConstants.PURPLE_CTDNA_SOMATIC_VCF_ID;

            if(!mSample.VcfTag.isEmpty())
                somaticVcf += mSample.VcfTag + ".";

            somaticVcf += "vcf.gz";
        }
        else
        {
            somaticVcf = mConfig.getSomaticVcf(mSample.TumorId);

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

        int matchedProbeCount = 0;
        if(mProbeVariants != null)
        {
            for(SomaticVariant variant : mVariants)
            {
                if(mProbeVariants.stream().anyMatch(x -> x.matches(variant)))
                {
                    ++matchedProbeCount;
                    variant.markProbeVariant();
                }
            }
        }

        CT_LOGGER.info("processed {} somatic variants from VCF({}), filtered({}) probeMatched({})",
                variantCount, filenamePart(somaticVcf), filteredCount, matchedProbeCount);

        return true;
    }

    private static final double NO_GC_RATIO = -1;

    private void processVariant(final List<String> targetSampleIds, final VariantContext variantContext)
    {
        VariantContextDecorator variant = new VariantContextDecorator(variantContext);

        double subclonalLikelihood = variant.context().getAttributeAsDouble(SUBCLONAL_LIKELIHOOD_FLAG, 0);
        double sequenceGcRatio = NO_GC_RATIO;

        if(mConfig.RefGenome != null && mSample.ApplyGcRatio)
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

                for(int i = 0; i <= VariantReadSupport.REALIGNED.ordinal(); ++i)
                {
                    qualTotal += Integer.parseInt(qualCounts[i]);
                }
            }

            if(umiTypeCounts != NO_UMI_COUNTS)
            {
                // override basic AD and DP if umit totals are set
                depth = umiTypeCounts.totalCount();
                alleleCount = umiTypeCounts.alleleCount();
            }

            somaticVariant.Samples.add(new GenotypeFragments(genotype.getSampleName(), alleleCount, depth, qualTotal, umiTypeCounts));
        }
    }

    public SomaticPurityResult processSample(final String sampleId, final PurityContext purityContext)
    {
        List<SomaticVariant> filteredVariants = Lists.newArrayList();

        int sampleTotalAD = 0;

        for(SomaticVariant variant : mVariants)
        {
            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);
            GenotypeFragments tumorFragData = variant.findGenotypeData(mSample.TumorId);

            if(sampleFragData == null || tumorFragData == null)
                continue;

            List<FilterReason> filterReasons = Lists.newArrayList(variant.filterReasons());

            boolean useForTotals = false;

            if(filterReasons.isEmpty())
            {
                // only include variants which satisfy the min avg qual check in the ctDNA sample
                if(sampleFragData.isLowQual())
                {
                    filterReasons.add(LOW_QUAL_PER_AD);
                }
                else
                {
                    useForTotals = true;
                }
            }

            if(mConfig.writeType(WriteType.SOMATIC_ALL) || useForTotals)
            {
                writeVariant(mSomaticWriter, mConfig, mSample, sampleId, variant, sampleFragData, tumorFragData, filterReasons);
            }

            if(!useForTotals)
                continue;

            filteredVariants.add(variant);

            sampleTotalAD += sampleFragData.AlleleCount;
        }

        if(filteredVariants.isEmpty())
            return INVALID_RESULT;

        // check for CHIP variants and remove them from variants used for purity estimates
        final double sampleAlleleTotal = sampleTotalAD;

        List<SomaticVariant> chipVariants = filteredVariants.stream()
                .filter(x -> isLikelyChipVariant(x, x.findGenotypeData(sampleId), sampleAlleleTotal)).collect(Collectors.toList());

        for(SomaticVariant variant : chipVariants)
        {
            CT_LOGGER.debug("sample({}) chip variant({}) ad({}) vs sampleTotal({})",
                    variant, variant.findGenotypeData(sampleId).AlleleCount, sampleTotalAD);

            filteredVariants.remove(variant);
            variant.addFilterReason(CHIP);
        }

        return mEstimator.calculatePurity(sampleId, purityContext, filteredVariants, mVariants.size(), chipVariants.size());
    }

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
        if(sequenceGcRatio != NO_GC_RATIO && mConfig.GcRatioMin > 0 && sequenceGcRatio < mConfig.GcRatioMin)
            filters.add(GC_RATIO);

        return filters;
    }

    private static boolean isLikelyChipVariant(
            final SomaticVariant variant, final GenotypeFragments sampleFragData, final double sampleAlleleTotal)
    {
        return sampleFragData.AlleleCount > CHIP_MIN_ALLELE_FRAGS
            && sampleFragData.AlleleCount / sampleAlleleTotal > CHIP_MIN_SAMPLE_PERC;
    }

    public static BufferedWriter initialiseVariantWriter(final PurityConfig config)
    {
        try
        {
            String fileName = config.formFilename(SOMATICS_FILE_ID);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            addCommonHeaderFields(sj, config);

            sj.add("Chromosome").add("Position").add("Ref").add("Alt").add("IsProbe");
            sj.add("Filter").add("Tier").add("Type").add("RepeatCount").add("Mappability").add("SubclonalPerc");
            sj.add("Gene").add("CodingEffect").add("Hotspot").add("Reported").add("VCN").add("CopyNumber");
            sj.add("TumorDP").add("TumorAD");
            sj.add("SampleDP").add("SampleAD").add("SampleDualDP").add("SampleDualAD").add("SampleQualPerAD").add("SeqGcRatio");

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to initialise variant output file: {}", e.toString());
            return null;
        }
    }

    private static synchronized void writeVariant(
            final BufferedWriter writer, final PurityConfig config,
            final SampleData sampleData, final String sampleId, final SomaticVariant variant,
            final GenotypeFragments sampleFragData, final GenotypeFragments tumorData, final List<FilterReason> filterReasons)
    {
        if(writer == null)
            return;

        try
        {
            VariantContextDecorator decorator = variant.decorator();
            VariantImpact variantImpact = decorator.variantImpact();

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            addCommonFields(sj, config, sampleData, sampleId);

            sj.add(variant.Chromosome).add(String.valueOf(variant.Position)).add(variant.Ref).add(variant.Alt);
            sj.add(String.valueOf(variant.isProbeVariant()));

            String filtersStr = filterReasons.stream().map(x -> x.toString()).collect(Collectors.joining(";"));

            if(filtersStr.isEmpty())
                filtersStr = PASS;

            sj.add(filtersStr).add(decorator.tier().toString()).add(variant.Type.toString()).add(String.valueOf(decorator.repeatCount()))
                    .add(format("%.2f", decorator.mappability())).add(format("%.2f", variant.SubclonalPerc));

            sj.add(variantImpact.CanonicalGeneName).add(variantImpact.CanonicalCodingEffect.toString())
                    .add(Hotspot.fromVariant(decorator.context()).toString()).add(String.valueOf(decorator.reported()))
                    .add(format("%.2f", decorator.variantCopyNumber())).add(format("%.2f", decorator.adjustedCopyNumber()));

            sj.add(String.valueOf(tumorData.Depth)).add(String.valueOf(tumorData.AlleleCount));
            sj.add(String.valueOf(sampleFragData.Depth)).add(String.valueOf(sampleFragData.AlleleCount));
            sj.add(String.valueOf(sampleFragData.UmiCounts.TotalDual)).add(String.valueOf(sampleFragData.UmiCounts.AlleleDual));
            sj.add(format("%.1f", sampleFragData.qualPerAlleleFragment()));
            sj.add(format("%.3f", variant.sequenceGcRatio()));

            writer.write(sj.toString());

            writer.newLine();
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write output file: {}", e.toString());
            System.exit(1);
        }
    }

    public static boolean plotSomaticVafs(final String patientId, final String sampleId, final PurityConfig config)
    {
        try
        {
            String summaryFile = config.formFilename(ResultsWriter.SUMMARY_FILE_ID);
            String somaticPeaksFile = config.formFilename(ResultsWriter.SOMATIC_PEAK_FILE_ID);
            //String somaticsFile = config.formFilename(SOMATICS_FILE_ID);

            if(!Files.exists(Paths.get(summaryFile)) || !Files.exists(Paths.get(somaticPeaksFile)))
            {
                CT_LOGGER.warn("plots missing required files: summary({}) somatics({})", summaryFile, somaticPeaksFile);
                return false;
            }

            int runCode = RExecutor.executeFromClasspath(
                    "plots/SomaticVafPlot.R", true, patientId, sampleId, summaryFile, somaticPeaksFile, config.PlotDir);

            return runCode == 0;
        }
        catch(Exception e)
        {
            CT_LOGGER.error("failed to generate CN plot with R script: {}", e.toString());
            return false;
        }
    }

}
