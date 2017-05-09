package com.hartwig.hmftools.breakpointinspector;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Set;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class BreakPointInspectorApplication {

    private static final String BAM_PATH = "bam";
    private static final String BREAK_POINT = "break";
    private static final String RANGE = "range";
    private static final String SV_LEN = "svlen";

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(BAM_PATH, true, "input BAM");
        options.addOption(BREAK_POINT, true, "position of break point in chrX:123456 format");
        options.addOption(RANGE, true, "base distance around breakpoint");
        options.addOption(SV_LEN, true, "length of the SV to inspect");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static String getFlagString(final Set<SAMFlag> flags) {
        ArrayList<String> names = new ArrayList<String>();
        for(final SAMFlag flag : flags)
        {
            names.add(flag.name());
        }
        return String.join("|", names);
    }

    public static void main(final String... args)
            throws ParseException, IOException, InvocationTargetException, NoSuchMethodException,
            InstantiationException, IllegalAccessException {

        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        // grab arguments
        final String bamPath = cmd.getOptionValue(BAM_PATH);
        final String breakPoint = cmd.getOptionValue(BREAK_POINT);
        final int range = Integer.parseInt(cmd.getOptionValue(RANGE, "500"));
        final int svLen = Integer.parseInt(cmd.getOptionValue(SV_LEN, "0"));

        if(bamPath == null || breakPoint == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Ecrf-Analyser", options);
            System.exit(1);
        }

        // parse breakpoint location
        final String[] split = breakPoint.split(":");
        if(split.length != 2)
        {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Ecrf-Analyser", options);
            System.exit(1);
        }
        final String chromosome = split[0];
        final int location = Integer.parseInt(split[1]);

        // load the file
        final File bamFile = new File(bamPath);
        final SamReader reader = SamReaderFactory.makeDefault().open(bamFile);

        // query the position
        final int index = reader.getFileHeader().getSequenceDictionary().getSequenceIndex(chromosome);
        final QueryInterval[] queryIntervals = {
                new QueryInterval(index, Math.max(0, location - range), location + range),
                new QueryInterval(index, Math.max(0, location + svLen - range), location + svLen + range)
        };

        // print  header
        System.out.println(String.join("\t",
                "Read Name",
                "Chromosome",
                "Position",
                "Mapping Quality",
                "Flags",
                "CIGAR",
                "Mate Chromosome",
                "Mate Position",
                "Read Length"
        ));

        // execute and parse the results
        CloseableIterator<SAMRecord> results = reader.query(queryIntervals, false);
        while(results.hasNext()) {
            final SAMRecord record = results.next();

            System.out.println(String.join("\t",
                    record.getReadName(),
                    record.getReferenceName(),
                    Integer.toString(record.getAlignmentStart()),
                    Integer.toString(record.getMappingQuality()),
                    getFlagString(record.getSAMFlags()),
                    record.getCigarString(),
                    record.getMateReferenceName(),
                    Integer.toString(record.getMateAlignmentStart()),
                    Integer.toString(record.getInferredInsertSize())
            ));
        }

    }
}
