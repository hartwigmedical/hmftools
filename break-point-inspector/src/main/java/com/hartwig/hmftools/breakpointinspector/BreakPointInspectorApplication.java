package com.hartwig.hmftools.breakpointinspector;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;

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

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(BAM_PATH, true, "input BAM");
        options.addOption(BREAK_POINT, true, "position of break point");
        options.addOption(RANGE, true, "base distance around breakpoint");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static BAMFileReader makeReader(@NotNull final String bamPath)
            throws IllegalAccessException, InvocationTargetException, InstantiationException, NoSuchMethodException {

        // extremely nasty work around for BAMFileReader constructors not being public
        final Class<BAMFileReader> klass =  BAMFileReader.class;
        final Class[] signature = {File.class, File.class, boolean.class, boolean.class, ValidationStringency.class, SAMRecordFactory.class};
        final Constructor<BAMFileReader> constructor = klass.getDeclaredConstructor(signature);
        constructor.setAccessible(true);

        // load the file
        final File bamFile = new File(bamPath);
        final File indexFile = new File(bamPath + ".bai"); // TODO: better java path handling?
        final BAMFileReader reader = constructor.newInstance(
                bamFile,
                indexFile,
                true,
                true,
                ValidationStringency.DEFAULT_STRINGENCY,
                DefaultSAMRecordFactory.getInstance()
        );

        return reader;
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

        // parse breakpoint location
        final String[] split = breakPoint.split(":");
        final String chromosome = split[0];
        final int location = Integer.parseInt(split[1]);

        final BAMFileReader reader = makeReader(bamPath);

        // query the position
        final int index = reader.getFileHeader().getSequenceDictionary().getSequenceIndex(chromosome);
        final QueryInterval[] queryIntervals = {
                new QueryInterval(index, Math.max(0, location - range), location + range)
        };
        CloseableIterator<SAMRecord> results = reader.query(queryIntervals, false);

        while(results.hasNext())
        {
            final SAMRecord record = results.next();
            System.out.println(record.getReadName());
        }

    }
}
