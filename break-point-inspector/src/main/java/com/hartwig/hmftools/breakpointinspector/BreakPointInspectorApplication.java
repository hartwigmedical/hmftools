package com.hartwig.hmftools.breakpointinspector;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.*;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import org.jetbrains.annotations.NotNull;

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
            throws ParseException, IOException {

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
        QueryInterval.optimizeIntervals(queryIntervals); // handle overlaps

        // print  header
        System.out.println(String.join("\t",
                "READ_NAME",
                "TEMPLATE_LENGTH",
                "ALIGNMENT",
                "CHROMOSOME",
                "POS",
                "MAPQ",
                "FLAGS",
                "CIGAR",
                "SEQ",
                "QUAL",
                "MATE_CHROMOSOME",
                "MATE_POS",
                "MATE_MAPQ",
                "MATE_FLAGS",
                "MATE_CIGAR",
                "MATE_SEQ",
                "MATE_QUAL"
        ));

        final Map<String, SAMRecord> firstReadMap = new Hashtable<>();

        int properPairReads = 0;
        int zeroQualityMappings = 0;

        // execute query and parse the results
        final SAMRecordIterator results = reader.query(queryIntervals, false);
        while(results.hasNext()) {
            final SAMRecord read = results.next();

            // TODO: should separate pairs from the start of the SV to the end of the SV
            if(read.getProperPairFlag()) {
                properPairReads++;
                continue;
            }
            if(read.getMappingQuality() == 0) {
                zeroQualityMappings++;
                continue;
            }

            final SAMRecord firstRead = firstReadMap.get(read.getReadName());
            if(firstRead != null) {

                assert(firstRead.getReadName().equals(read.getReadName()));
                assert(firstRead.getMateAlignmentStart() == read.getAlignmentStart());
                assert(firstRead.getMateReferenceName().equals(read.getReferenceName()));

                // TODO: sort output by firstRead position
                // TODO: or just print record.getSAMString() ??
                System.out.println(String.join("\t",
                        // first read
                        firstRead.getReadName(),
                        Integer.toString(firstRead.getInferredInsertSize()),
                        String.join(",", SamPairUtil.getPairOrientation(firstRead).toString(), SamPairUtil.getPairOrientation(read).toString()),
                        firstRead.getReferenceName(),
                        Integer.toString(firstRead.getAlignmentStart()),
                        Integer.toString(firstRead.getMappingQuality()),
                        getFlagString(firstRead.getSAMFlags()),
                        firstRead.getCigarString(),
                        firstRead.getReadString(),
                        firstRead.getBaseQualityString(),
                        // second read
                        read.getReferenceName(),
                        Integer.toString(read.getAlignmentStart()),
                        Integer.toString(read.getMappingQuality()),
                        getFlagString(read.getSAMFlags()),
                        read.getCigarString(),
                        read.getReadString(),
                        read.getBaseQualityString()
                ));
            } else {
                firstReadMap.put(read.getReadName(), read);
            }
        }

        System.out.println();
        System.out.println("STATS");
        System.out.println("PROPER_PAIRS\t" + properPairReads / 2);
        System.out.println("ZERO_MAPQ\t" + zeroQualityMappings);

    }
}
