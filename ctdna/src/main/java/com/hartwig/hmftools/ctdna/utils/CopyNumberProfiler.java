package com.hartwig.hmftools.ctdna.utils;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.ctdna.purity.CopyNumberGcData;
import com.hartwig.hmftools.ctdna.purity.CopyNumberProfile;
import com.hartwig.hmftools.ctdna.purity.PurityConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class CopyNumberProfiler
{
    private static final String SAMPLE = "sample";
    private static final String COBALT_SAMPLE = "cobalt_sample";

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Sample ID");
        options.addOption(COBALT_SAMPLE, true, "Cobalt sample ID (eg for the ctDNA sample)");
        PurityConfig.addCommandLineOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        PurityConfig purityConfig = new PurityConfig(cmd);
        String sampleId = cmd.getOptionValue(SAMPLE);
        String cobaltSampleId = cmd.getOptionValue(COBALT_SAMPLE, sampleId);

        CopyNumberProfile cnProfile = new CopyNumberProfile(purityConfig);
        cnProfile.processSample(sampleId, cobaltSampleId);

        try
        {
            String fileName = purityConfig.OutputDir + sampleId + ".cn_segment_data.csv";

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("Chromosome,SegmentStart,SegmentEnd,CopyNumber,GcRatioCount,GcRatioMedian,GcRatioMean");
            writer.newLine();

            for(CopyNumberGcData cnSegment : cnProfile.copyNumberGcRatios())
            {
                writer.write(format("%s,%d,%d,%.2f,%d,%.4f,%.4f",
                        cnSegment.Chromosome, cnSegment.SegmentStart, cnSegment.SegmentEnd, cnSegment.CopyNumber,
                        cnSegment.count(), cnSegment.median(), cnSegment.mean()));
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write copy number segment file: {}", e.toString());
        }


        CT_LOGGER.info("Sample VCF analyser complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
