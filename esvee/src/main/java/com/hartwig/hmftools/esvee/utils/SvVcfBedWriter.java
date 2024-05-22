package com.hartwig.hmftools.esvee.utils;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;

import java.io.BufferedWriter;
import java.nio.file.Files;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class SvVcfBedWriter
{
    private static final String INPUT_VCF = "input_vcf";
    private static final String OUTPUT_BED = "output_bed";

    private static final int BREAKEND_PROXIMITY = 10;

    public static void writeBedFromSvVcf(final ConfigBuilder configBuilder)
    {
        String inputVcf = configBuilder.getValue(INPUT_VCF);
        String outputBed = configBuilder.getValue(OUTPUT_BED);

        if(inputVcf == null || outputBed == null || !Files.exists(Paths.get(inputVcf)))
        {
            SV_LOGGER.error("invalid input VCF or missing output BED config");
            System.exit(1);
        }

        // SV_LOGGER.info("SV Append for {} SV breakends", mChrBreakendMap.values().stream().mapToInt(x -> x.size()).sum());
        int breakendCount = 0;

        try
        {
            VcfFileReader vcfFileReader = new VcfFileReader(inputVcf);

            if(!vcfFileReader.fileValid())
                System.exit(1);

            SV_LOGGER.info("reading input VCF {}", inputVcf);

            BufferedWriter writer = createBufferedWriter(outputBed, false);

            int startPosition = 0;
            int lastPosition = 0;
            String currentChromosome = "";

            for(VariantContext variant : vcfFileReader.iterator())
            {
                ++breakendCount;

                String chromosome = variant.getContig();
                int breakendPosition = variant.getStart();

                if(chromosome.equals(currentChromosome) && lastPosition > 0 && breakendPosition - lastPosition <= BREAKEND_PROXIMITY)
                {
                    lastPosition = breakendPosition;
                    continue;
                }

                if(lastPosition > 0)
                {
                    writer.write(format("%s\t%d\t%d",
                            currentChromosome, startPosition - BREAKEND_PROXIMITY, lastPosition + BREAKEND_PROXIMITY));
                    writer.newLine();
                }

                currentChromosome = chromosome;
                lastPosition = breakendPosition;
                startPosition = breakendPosition;
            }

            if(lastPosition > 0)
            {
                writer.write(format("%s\t%d\t%d",
                        currentChromosome, startPosition - BREAKEND_PROXIMITY, lastPosition + BREAKEND_PROXIMITY));
                writer.newLine();
            }

            vcfFileReader.close();
            writer.close();
        }
        catch(Exception e)
        {
            SV_LOGGER.error("error reading / writing: {}", e.toString());
            System.exit(1);
        }

        SV_LOGGER.info("wrote {} SV breakends entries to {}", breakendCount, outputBed);
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addPath(INPUT_VCF, true, "Input VCF");
        configBuilder.addConfigItem(OUTPUT_BED, true, "Output BED");
        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        writeBedFromSvVcf(configBuilder);
    }
}
