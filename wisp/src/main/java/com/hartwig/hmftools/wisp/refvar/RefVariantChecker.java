package com.hartwig.hmftools.wisp.refvar;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS;
import static com.hartwig.hmftools.wisp.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.wisp.purity.SampleData;
import com.hartwig.hmftools.wisp.purity.variant.UmiTypeCounts;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class RefVariantChecker
{
    private final String mOutputDir;
    private final String mOutputId;
    private final String mSampleVcfs;
    private final String mPurpleDir;

    private final List<SampleData> mSamples;

    private final BufferedWriter mWriter;

    private static final String SAMPLE_VCFS = "sample_vcf";

    public RefVariantChecker(final ConfigBuilder configBuilder)
    {
        mSampleVcfs = configBuilder.getValue(SAMPLE_VCFS);
        mPurpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
        mOutputDir = parseOutputDir(configBuilder);
        mOutputId = configBuilder.getValue(OUTPUT_ID);

        mSamples = SampleData.loadSampleDataFile(configBuilder.getValue(SAMPLE_ID_FILE));

        mWriter = initialiseWriter();
    }

    public void run()
    {
        if(mSamples.isEmpty() || mWriter == null)
            System.exit(1);

        CT_LOGGER.info("running probe variant selection for {} patients and {} samples",
                mSamples.size(), mSamples.stream().mapToInt(x -> x.SampleIds.size()).sum());

        for(SampleData sample : mSamples)
        {
            processSample(sample);
        }

        closeBufferedWriter(mWriter);

        CT_LOGGER.info("Reference variant checking complete");
    }

    private void processSample(final SampleData sample)
    {
        // read in Purple variants

        String purpleDir = convertWildcardSamplePath(mPurpleDir, sample.TumorId);
        String purpleSomaticVcf = PurpleCommon.purpleSomaticVcfFile(purpleDir, sample.TumorId);
        String purpleGermlineVcf = PurpleCommon.purpleGermlineVcfFile(purpleDir, sample.TumorId);

        Map<String,List<VariantContextDecorator>> tumorChrVariantMap = Maps.newHashMap();
        Map<String,List<VariantContextDecorator>> germlineChrVariantMap = Maps.newHashMap();

        loadPurpleVariants(purpleSomaticVcf, tumorChrVariantMap, sample.PatientId);
        loadPurpleVariants(purpleGermlineVcf, germlineChrVariantMap, sample.PatientId);

        // now process each sample in turn
        for(String sampleId : sample.SampleIds)
        {
            String sampleVcf = convertWildcardSamplePath(mSampleVcfs, sampleId);

            VcfFileReader sampleFileReader = new VcfFileReader(sampleVcf);

            if(!sampleFileReader.fileValid())
            {
                CT_LOGGER.error("failed to read sample vcf({})", sampleVcf);
                continue;
            }

            int sampleVarCount = 0;

            for(VariantContext variantContext : sampleFileReader.iterator())
            {
                if(variantContext.isFiltered())
                    continue;

                VariantContextDecorator variant = new VariantContextDecorator(variantContext);

                // if(variant.tier() != VariantTier.HOTSPOT && variant.tier() != VariantTier.PANEL)
                //    continue;

                ++sampleVarCount;

                VariantContextDecorator tumorVariant = findExistingVariant(variant, tumorChrVariantMap);
                VariantContextDecorator germlineVariant = findExistingVariant(variant, germlineChrVariantMap);

                writeVariant(sample.PatientId, sampleId, variant, tumorVariant, germlineVariant);

                CT_LOGGER.info("patient({}) sample({}) processed {} variants", sample.PatientId, sampleId, sampleVarCount);
            }
        }
    }

    private void loadPurpleVariants(
            final String purpleVcf, final Map<String,List<VariantContextDecorator>> chrVariantMap, final String patientId)
    {
        VcfFileReader vcfFileReader = new VcfFileReader(purpleVcf);

        if(!vcfFileReader.fileValid())
        {
            CT_LOGGER.warn("failed to read Purple vcf({})", purpleVcf);
            return;
        }

        int variantCount = 0;
        for(VariantContext variantContext : vcfFileReader.iterator())
        {
            if(variantContext.isFiltered())
                continue;

            List<VariantContextDecorator> chrVariants = chrVariantMap.get(variantContext.getContig());

            if(chrVariants == null)
            {
                chrVariants = Lists.newArrayList();
                chrVariantMap.put(variantContext.getContig(), chrVariants);
            }

            chrVariants.add(new VariantContextDecorator(variantContext));
            ++variantCount;
        }

        CT_LOGGER.debug("patient({}) loaded {} variants from Purple VCF({})", patientId, variantCount, purpleVcf);
    }

    private VariantContextDecorator findExistingVariant(
            final VariantContextDecorator refVariant, final Map<String,List<VariantContextDecorator>> chrVariantMap)
    {
        List<VariantContextDecorator> tumorVariants = chrVariantMap.get(refVariant.chromosome());

        if(tumorVariants == null)
            return null;

        return tumorVariants.stream()
                .filter(x -> x.position() == refVariant.position())
                .filter(x -> x.ref().equals(refVariant.ref()))
                .filter(x -> x.alt().equals(refVariant.alt()))
                .findFirst().orElse(null);
    }

    private BufferedWriter initialiseWriter()
    {
        if(mOutputDir == null)
            return null;

        try
        {
            String filename = mOutputDir + "cohort_ref_variant_matching";

            if(mOutputId != null)
                filename += "." + mOutputId;

            filename += TSV_EXTENSION;

            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("PatientId").add("SampleId");
            sj.add("Chromosome").add("Position").add("Ref").add("Alt").add("Type").add("Tier").add("Filter");
            sj.add("InTumor").add("TumorReported").add("TumorQual").add("TumorVaf").add("TumorHotspot");
            sj.add("InGermline").add("GermlineReported").add("GermlineQual").add("GermlineVaf").add("GermlineHotspot");
            sj.add("Gene").add("CodingEffect");
            sj.add("SampleDP").add("SampleAD").add("SampleDual").add("SampleAlleleDual");

            writer.write(sj.toString());
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    private synchronized void writeVariant(
            final String patientId, final String sampleId, final VariantContextDecorator variant,
            final VariantContextDecorator tumorVariant, final VariantContextDecorator germlineVariant)
    {
        if(mWriter == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(patientId).add(sampleId);
            sj.add(variant.chromosome()).add(String.valueOf(variant.position())).add(variant.ref()).add(variant.alt());
            sj.add(variant.type().toString()).add(variant.tier().toString()).add(variant.filter());

            for(int i = 0; i <= 1; ++i)
            {
                VariantContextDecorator existingVariant = (i == 0) ? tumorVariant : germlineVariant;

                if(existingVariant != null)
                {
                    sj.add("true");
                    sj.add(String.valueOf(existingVariant.reported()));
                    sj.add(format("%.0f", existingVariant.qual()));
                    sj.add(format("%.2f", existingVariant.adjustedVaf()));
                    sj.add(existingVariant.hotspot().toString());
                }
                else
                {
                    sj.add("false").add("false").add("-1").add("-1").add("N/A");
                }
            }

            sj.add(variant.gene());
            sj.add(variant.canonicalCodingEffect().toString());

            Genotype genotype = variant.context().getGenotype(sampleId);
            UmiTypeCounts umiTypeCounts = UmiTypeCounts.fromAttribute(genotype.getExtendedAttribute(UMI_TYPE_COUNTS, null));
            sj.add(String.valueOf(umiTypeCounts.totalCount()));
            sj.add(String.valueOf(umiTypeCounts.alleleCount()));
            sj.add(String.valueOf(umiTypeCounts.TotalDual));
            sj.add(String.valueOf(umiTypeCounts.AlleleDual));

            mWriter.write(sj.toString());
            mWriter.newLine();
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write variants: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        addSampleIdFile(configBuilder, true);
        configBuilder.addPath(SAMPLE_VCFS, true, "Sample VCFs");
        configBuilder.addPath(PURPLE_DIR_CFG, true, PURPLE_DIR_DESC);
        addOutputOptions(configBuilder, true);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);
        setLogLevel(configBuilder);

        RefVariantChecker refVariantChecker = new RefVariantChecker(configBuilder);
        refVariantChecker.run();
    }
}
