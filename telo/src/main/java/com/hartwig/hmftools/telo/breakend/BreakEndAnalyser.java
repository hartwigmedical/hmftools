package com.hartwig.hmftools.telo.breakend;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.telo.TeloConfig.TE_LOGGER;

import java.io.File;
import java.time.Instant;
import java.time.Duration;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BreakEndAnalyser
{
    private final BreakEndAnalyserConfig mConfig;

    public BreakEndAnalyser(final Options options, final String... args) throws ParseException
    {
        VersionInfo versionInfo = new VersionInfo("telo.version");
        TE_LOGGER.info("Telo version: {}", versionInfo.version());

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        mConfig = new BreakEndAnalyserConfig(cmd);
    }

    private void run()
    {
        if(!mConfig.isValid())
        {
            TE_LOGGER.error(" invalid config, exiting");
            System.exit(1);
        }

        TE_LOGGER.info("starting telomeric analysis");
        Instant start = Instant.now();

        processBam(mConfig);

        Instant finish = Instant.now();
        long seconds = Duration.between(start, finish).getSeconds();
        TE_LOGGER.info("Telo run complete, time taken: {}m {}s",  seconds / 60, seconds % 60);
    }

    public static void main(final String... args)
    {
        final Options options = BreakEndAnalyserConfig.createOptions();
        BreakEndAnalyser application = null;

        try
        {
            application = new BreakEndAnalyser(options, args);
        }
        catch(ParseException e)
        {
            TE_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("TeloApplication", options);
            System.exit(1);
        }

        application.run();
    }

    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private void processBam(@NotNull BreakEndAnalyserConfig config)
    {
        SamReaderFactory factory = SamReaderFactory.makeDefault();
        SamReader samReader = factory.open(new File(config.TelbamFile));

        TelomericSplitReadAnalyser breakEndAnalyser = new TelomericSplitReadAnalyser();

        //TE_LOGGER.info("processing region({})", baseRegion.toString());

        // do not change the follow line to use functions other than query without testing that unmapped reads are returned
        try (final SAMRecordIterator iterator = samReader.iterator())
        {
            while (iterator.hasNext())
            {
                final SAMRecord record = iterator.next();
                processReadRecord(record, breakEndAnalyser);
            }
        }

        breakEndAnalyser.consolidatePotentialBreakEnds();

        //TE_LOGGER.info("processed region({}) read count({})", baseRegion.toString(), mReadCount);
    }

    // we go through the records twice.
    // first parse we identify candidate soft clips, i.e. positions that are one side telomeric and the
    // other side mapped to genome
    // then we tidy up all these candidates and produce a report
    private void processReadRecord(@NotNull SAMRecord record, @NotNull TelomericSplitReadAnalyser analyser)
    {
        // work out if the record is a telomere
        analyser.processRead(record);
    }
}
