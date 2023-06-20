package com.hartwig.hmftools.purple.tools;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadDelimitedIdFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.PaveVcfTags.GNOMAD_FREQ;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.BufferedWriter;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.variant.variantcontext.VariantContext;

public class SomaticVariantCompiler
{
    private final List<String> mSampleIds;
    private final List<String> mRequiredFields;
    private final String mOutputFile;
    private final String mPurpleDir;

    private static final String OUTPUT_FILE = "output_file";
    private static final String PURPLE_DIR = "purple_dir";
    private static final String REQUIRED_FIELDS = "required_fields";

    public SomaticVariantCompiler(final CommandLine cmd)
    {
        mSampleIds = loadSampleIdsFile(cmd);
        mOutputFile = cmd.getOptionValue(OUTPUT_FILE);
        mPurpleDir = cmd.getOptionValue(PURPLE_DIR);

        mRequiredFields = Arrays.stream(cmd.getOptionValue(REQUIRED_FIELDS).split(",", -1)).collect(Collectors.toList());
    }

    public void run()
    {
        PPL_LOGGER.info("compiling somatic variants for {} samples", mSampleIds.size());

        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add("SampleId").add("VariantInfo").add("Type").add("Tier").add("Qual").add("Filters");

            for(String field : mRequiredFields)
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

            for(String field : mRequiredFields)
            {
                String fieldValue = variantContext.getAttributeAsString(field, "");
                sj.add(fieldValue);
            }
            
            writer.write(sj.toString());
            writer.newLine();

            ++variantCount;
        }

        PPL_LOGGER.debug("sample({}) loaded {} variants", sampleId, variantCount);
    }

    public static void main(final String[] args) throws ParseException
    {
        final Options options = new Options();
        addLoggingOptions(options);
        addSampleIdFile(options);
        options.addOption(OUTPUT_FILE, true, "Output filename");
        options.addOption(PURPLE_DIR, true, "Path to purple files");
        options.addOption(REQUIRED_FIELDS, true, "Required VCF fields separated by ','");

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        SomaticVariantCompiler somaticVariantCompiler = new SomaticVariantCompiler(cmd);
        somaticVariantCompiler.run();
    }

    private static CommandLine createCommandLine(final String[] args, final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
