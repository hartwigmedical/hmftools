package com.hartwig.hmftools.purple.tools;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.BufferedWriter;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class SomaticVariantCompiler
{
    private final List<String> mSampleIds;
    private final List<String> mCommonFields;
    private final List<String> mGenotypeFields;
    private final String mOutputFile;
    private final String mPurpleDir;

    private static final String OUTPUT_FILE = "output_file";
    private static final String COMMON_FIELDS = "common_fields";
    private static final String GENOTYPE_FIELDS = "genotype_fields";

    public SomaticVariantCompiler(final ConfigBuilder configBuilder)
    {
        mSampleIds = loadSampleIdsFile(configBuilder.getValue(SAMPLE_ID_FILE));
        mOutputFile = configBuilder.getValue(OUTPUT_FILE);
        mPurpleDir = configBuilder.getValue(PURPLE_DIR_CFG);

        mCommonFields = Arrays.stream(configBuilder.getValue(COMMON_FIELDS).split(",", -1)).collect(Collectors.toList());
        mGenotypeFields = Arrays.stream(configBuilder.getValue(GENOTYPE_FIELDS).split(",", -1)).collect(Collectors.toList());
    }

    public void run()
    {
        PPL_LOGGER.info("compiling somatic variants for {} samples", mSampleIds.size());

        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add("SampleId").add("VariantInfo").add("Type").add("Tier").add("Qual").add("Filters");

            for(String field : mCommonFields)
            {
                sj.add(field);
            }

            for(String field : mGenotypeFields)
            {
                sj.add(field);
            }

            writer.write(sj.toString());
            writer.newLine();

            for(String sampleId : mSampleIds)
            {
                processSampleVariants(sampleId, writer);
            }

            writer.close();
        }
        catch(Exception e)
        {
            PPL_LOGGER.error("failed to write variant output file: {}", e.toString());
            e.printStackTrace();
            System.exit(1);
        }

        PPL_LOGGER.info("compilation complete");
    }

    private void processSampleVariants(final String sampleId, final BufferedWriter writer) throws Exception
    {
        String purpleDir = mPurpleDir.replaceAll("\\*", sampleId);
        String purpleVcf = PurpleCommon.purpleSomaticVcfFile(purpleDir, sampleId);

        VcfFileReader vcfFileReader = new VcfFileReader(purpleVcf);

        if(!vcfFileReader.fileValid())
        {
            return;
        }

        int variantCount = 0;

        for(VariantContext variantContext : vcfFileReader.iterator())
        {
            VariantContextDecorator variant = new VariantContextDecorator(variantContext);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(sampleId);
            sj.add(format("%s:%d %s>%s", variant.chromosome(), variant.position(), variant.ref(), variant.alt()));
            sj.add(variant.type().toString());
            sj.add(variant.tier().toString());
            sj.add(String.valueOf(variant.qual()));
            sj.add(variant.filter());

            for(String field : mCommonFields)
            {
                String fieldValue = variantContext.getAttributeAsString(field, "");
                sj.add(fieldValue);
            }

            for(String field : mGenotypeFields)
            {
                StringJoiner gtSj = new StringJoiner(ITEM_DELIM);

                for(Genotype genotypeContext : variantContext.getGenotypes())
                {
                    if(field.equals("AD"))
                    {
                        gtSj.add(String.valueOf(genotypeContext.getAD()[1]));
                    }
                    else if(field.equals("DP"))
                    {
                        gtSj.add(String.valueOf(genotypeContext.getDP()));
                    }
                    else
                    {
                        Object fieldValue = genotypeContext.getExtendedAttribute(field, null);
                        String fieldStr = fieldValue != null ? fieldValue.toString() : "";
                        gtSj.add(fieldStr);
                    }
                }

                sj.add(gtSj.toString());
            }

            writer.write(sj.toString());
            writer.newLine();

            ++variantCount;
        }

        PPL_LOGGER.debug("sample({}) loaded {} variants", sampleId, variantCount);
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        addSampleIdFile(configBuilder, true);
        configBuilder.addPathItem(OUTPUT_FILE, true, "Output filename");
        configBuilder.addConfigItem(COMMON_FIELDS, false, "Required VCF fields separated by ','");
        configBuilder.addConfigItem(GENOTYPE_FIELDS, false, "Required VCF genotype fields separated by ','");
        configBuilder.addConfigItem(PURPLE_DIR_CFG, true, PURPLE_DIR_DESC);
        addLoggingOptions(configBuilder);

        if(!configBuilder.parseCommandLine(args))
        {
            configBuilder.logItems();
            System.exit(1);
        }

        setLogLevel(configBuilder);

        SomaticVariantCompiler somaticVariantCompiler = new SomaticVariantCompiler(configBuilder);
        somaticVariantCompiler.run();
    }
}
