package com.hartwig.hmftools.ctdna.refvar;

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
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.ctdna.common.SampleData;
import com.hartwig.hmftools.ctdna.purity.variant.UmiTypeCounts;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class RefVariantChecker
{
    private final String mOutputDir;
    private final String mOutputId;
    private final String mCtDnaVcfs;
    private final String mPurpleDir;
    private final List<SampleData> mSamples;

    private final BufferedWriter mWriter;

    private static final String CTDNA_VCFS = "ctdna_vcf";

    public RefVariantChecker(final ConfigBuilder configBuilder)
    {
        mCtDnaVcfs = configBuilder.getValue(CTDNA_VCFS);
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
                mSamples.size(), mSamples.stream().mapToInt(x -> x.CtDnaSamples.size()).sum());

        for(SampleData sample : mSamples)
        {
            processSample(sample);
        }

        closeBufferedWriter(mWriter);

        CT_LOGGER.info("Probe variation selection complete");
    }

    private void processSample(final SampleData sample)
    {
        // read in Purple variants

        String purpleDir = convertWildcardSamplePath(mPurpleDir, sample.TumorId);
        String purpleVcf = PurpleCommon.purpleSomaticVcfFile(purpleDir, sample.TumorId);

        VcfFileReader vcfFileReader = new VcfFileReader(purpleVcf);

        if(!vcfFileReader.fileValid())
        {
            CT_LOGGER.error("failed to read Purple vcf({})", purpleVcf);
            return;
        }

        Map<String,List<VariantContextDecorator>> tumorChrVariants = Maps.newHashMap();

        int tumorVarCount = 0;
        for(VariantContext variantContext : vcfFileReader.iterator())
        {
            if(variantContext.isFiltered())
                continue;

            List<VariantContextDecorator> chrVariants = tumorChrVariants.get(variantContext.getContig());

            if(chrVariants == null)
            {
                chrVariants = Lists.newArrayList();
                tumorChrVariants.put(variantContext.getContig(), chrVariants);
            }

            chrVariants.add(new VariantContextDecorator(variantContext));
            ++tumorVarCount;
        }

        CT_LOGGER.info("patient({}) read {} Purple variants from vcf({})", sample.PatientId, tumorVarCount, purpleVcf);

        // now process each ctDNA sample in turn
        for(String ctDnaSampleId : sample.CtDnaSamples)
        {
            String ctDnaSampleVcf = convertWildcardSamplePath(mCtDnaVcfs, ctDnaSampleId);

            VcfFileReader sampleFileReader = new VcfFileReader(ctDnaSampleVcf);

            if(!sampleFileReader.fileValid())
            {
                CT_LOGGER.error("failed to read ctDna vcf({})", ctDnaSampleVcf);
                continue;
            }

            //List<VariantContextDecorator> ctDnaVariants = Lists.newArrayList();

            int ctdnaVarCount = 0;

            for(VariantContext variantContext : sampleFileReader.iterator())
            {
                if(variantContext.isFiltered())
                    continue;

                VariantContextDecorator variant = new VariantContextDecorator(variantContext);

                if(variant.tier() != VariantTier.HOTSPOT && variant.tier() != VariantTier.PANEL)
                    continue;

                // ctDnaVariants.add(variant);
                ++ctdnaVarCount;

                VariantContextDecorator tumorVariant = null;

                List<VariantContextDecorator> tumorVariants = tumorChrVariants.get(variant.chromosome());

                if(tumorVariants != null)
                {
                    tumorVariant = tumorVariants.stream()
                            .filter(x -> x.position() == variant.position())
                            .filter(x -> x.ref().equals(variant.ref()))
                            .filter(x -> x.alt().equals(variant.alt()))
                            .findFirst().orElse(null);
                }

                writeVariant(sample.PatientId, ctDnaSampleId, variant, tumorVariant);

                CT_LOGGER.info("patient({}) sample({}) processed {} variants", sample.PatientId, ctDnaSampleId, ctdnaVarCount);
            }
        }
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
            sj.add("InTumor").add("Reported").add("TumorQual").add("TumorVaf").add("Hotspot");
            sj.add("Gene").add("CodingEffect");
            sj.add("SampleDP").add("SampleAD").add("SampleRefDual").add("SampleAlleleDual");

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
            final String patientId, final String sampleId,
            final VariantContextDecorator variant, final VariantContextDecorator tumorVariant)
    {
        if(mWriter == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(patientId).add(sampleId);
            sj.add(variant.chromosome()).add(String.valueOf(variant.position())).add(variant.ref()).add(variant.alt());
            sj.add(variant.type().toString()).add(variant.tier().toString()).add(variant.filter());

            if(tumorVariant != null)
            {
                sj.add("true");
                sj.add(String.valueOf(tumorVariant.reported()));
                sj.add(format("%.0f", tumorVariant.qual()));
                sj.add(format("%.2f", tumorVariant.adjustedVaf()));
                sj.add(tumorVariant.hotspot().toString());
            }
            else
            {
                sj.add("false").add("false").add("-1").add("-1").add("N/A");
            }

            sj.add(variant.gene());
            sj.add(variant.canonicalCodingEffect().toString());

            Genotype genotype = variant.context().getGenotype(sampleId);
            UmiTypeCounts umiTypeCounts = UmiTypeCounts.fromAttribute(genotype.getExtendedAttribute(UMI_TYPE_COUNTS, null));
            sj.add(String.valueOf(umiTypeCounts.refTotal()));
            sj.add(String.valueOf(umiTypeCounts.alleleTotal()));
            sj.add(String.valueOf(umiTypeCounts.dualTotal()));
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
        final VersionInfo version = new VersionInfo("ctdna.version");
        CT_LOGGER.info("RefVariantChecker version: {}", version.version());

        ConfigBuilder configBuilder = new ConfigBuilder();

        addSampleIdFile(configBuilder, true);
        configBuilder.addPath(CTDNA_VCFS, true, "CtDNA VCFs");
        configBuilder.addPath(PURPLE_DIR_CFG, true, PURPLE_DIR_DESC);
        addOutputOptions(configBuilder, true);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);
        setLogLevel(configBuilder);

        RefVariantChecker refVariantChecker = new RefVariantChecker(configBuilder);
        refVariantChecker.run();
    }
}
