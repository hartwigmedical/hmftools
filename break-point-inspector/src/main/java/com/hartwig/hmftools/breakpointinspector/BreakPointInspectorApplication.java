package com.hartwig.hmftools.breakpointinspector;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

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
    private static final String INCLUDE_PROPER = "proper";

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(BAM_PATH, true, "input BAM");
        options.addOption(BREAK_POINT, true, "position of break point in chrX:123456 format");
        options.addOption(RANGE, true, "base distance around breakpoint");
        options.addOption(SV_LEN, true, "length of the SV to inspect");
        options.addOption(INCLUDE_PROPER, false, "include proper pairs in output");
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
        for (final SAMFlag flag : flags) {
            names.add(flag.name());
        }
        return String.join("|", names);
    }

    private static String getOrientationString(final SAMRecord read) {
        return read.getReadPairedFlag() && !(read.getReadUnmappedFlag() || read.getMateUnmappedFlag()) ?
                SamPairUtil.getPairOrientation(read).toString() :
                "unaligned";
    }

    private static void printHelpAndExit(final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Break-Point-Inspector", "header", options, "footer", true);
        System.exit(1);
    }

    private static boolean doesSpanLocation(final SAMRecord first, final SAMRecord second, final int location) {
        assert first.getInferredInsertSize() == -second.getInferredInsertSize();
        return (first.getAlignmentStart() + first.getInferredInsertSize()) > location;
    }

    private static boolean doesSpanLocation(final SAMRecord read, final int location) {
        return (read.getAlignmentStart() - location) * (read.getAlignmentEnd() - location) < 0;
    }

    private static boolean doesSpanLocationUnclipped(final SAMRecord read, final int location) {
        return (read.getUnclippedStart() - location) * (read.getUnclippedEnd() - location) < 0;
    }

    public static void main(final String... args) throws ParseException, IOException {

        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        // grab arguments
        final String bamPath = cmd.getOptionValue(BAM_PATH);
        final String breakPoint = cmd.getOptionValue(BREAK_POINT);
        final int range = Integer.parseInt(cmd.getOptionValue(RANGE, "500"));
        final int svLen = Integer.parseInt(cmd.getOptionValue(SV_LEN, "0"));
        final boolean outputProperPairs = cmd.hasOption(INCLUDE_PROPER);

        if (bamPath == null || breakPoint == null) {
            printHelpAndExit(options);
        }

        // parse breakpoint location
        final String[] split = breakPoint.split(":");
        if (split.length != 2) {
            printHelpAndExit(options);
        }
        final String chromosome = split[0];
        final int location1 = Integer.parseInt(split[1]);
        final int location2 = location1 + svLen;

        // load the file
        final File bamFile = new File(bamPath);
        final SamReader reader = SamReaderFactory.makeDefault().open(bamFile);

        // query the position
        final int index = reader.getFileHeader().getSequenceDictionary().getSequenceIndex(chromosome);
        final QueryInterval[] queryIntervals = {
                new QueryInterval(index, Math.max(0, location1 - range), location1 + range),
                new QueryInterval(index, Math.max(0, location2 - range), location2 + range) };
        QueryInterval.optimizeIntervals(queryIntervals); // TODO: this doesn't do what I think it does

        // print  header
        System.out.println(
                String.join("\t", "READ_NAME", "TEMPLATE_LENGTH", "ALIGNMENT", "CHROMOSOME", "POS", "MAPQ", "FLAGS",
                        "CIGAR", "SEQ", "QUAL", "MATE_CHROMOSOME", "MATE_POS", "MATE_MAPQ", "MATE_FLAGS", "MATE_CIGAR",
                        "MATE_SEQ", "MATE_QUAL"));

        final Map<String, SAMRecord> readMap = new Hashtable<>();
        final Map<Integer, String> outputMap = new TreeMap<>();

        int normalReadsInLocation1 = 0;
        int normalReadsInLocation2 = 0;

        int interestingReads = 0;
        int zeroQualityMappings = 0;

        // execute query and parse the results
        final SAMRecordIterator results = reader.query(queryIntervals, false);
        while (results.hasNext()) {
            final SAMRecord read = results.next();

            if (read.getMappingQuality() == 0) {
                zeroQualityMappings++;
                // continue;
            }

            final boolean spanLocation1 = doesSpanLocation(read, location1);
            final boolean spanLocation2 = doesSpanLocation(read, location2);
            if (read.getProperPairFlag()) {
                if (spanLocation1)
                    normalReadsInLocation1++;
                if (spanLocation2)
                    normalReadsInLocation2++;
            }

            final SAMRecord firstRead = readMap.get(read.getReadName());
            if (firstRead == null) {
                // this is the first read we've seen of the pair
                readMap.put(read.getReadName(), read);
            } else {
                // we've found the mate
                readMap.remove(read.getReadName());

                assert firstRead.getReadName().equals(read.getReadName());
                assert firstRead.getMateReferenceName().equals(read.getReferenceName());
                assert firstRead.getMateAlignmentStart() == read.getAlignmentStart();
                assert firstRead.getProperPairFlag() == read.getProperPairFlag();

                // filter for the reads we want
                final boolean differentChromosomes = !read.getReferenceName().equals(read.getMateReferenceName());
                final boolean oneOfPairUnmapped = read.getReadUnmappedFlag() ^ read.getMateUnmappedFlag();
                final boolean isSecondaryOrSupplementary = read.getSupplementaryAlignmentFlag();
                final boolean isProper = read.getProperPairFlag();

                boolean output = false;
                if (!isProper) {
                    interestingReads += 2;
                    output = true;
                }

                if (!output)
                    continue;

                // TODO: sort output by firstRead position
                // @formatter:off
                outputMap.put(firstRead.getAlignmentStart(),
                        String.join("\t",
                            // first read
                            firstRead.getReadName(),
                            Integer.toString(firstRead.getInferredInsertSize()),
                            String.join(",", getOrientationString(firstRead), getOrientationString(read)), // TODO: I believe these should always be the same value
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
                // @formatter:on
            }
        }

        // print unpaired reads - may still have mate data
        for (Map.Entry<String, SAMRecord> entry : readMap.entrySet()) {
            final SAMRecord read = entry.getValue();

            // output single reads which are paired
            boolean output = false;
            final boolean spanLocation1 = doesSpanLocation(read, location1);
            final boolean spanLocation2 = doesSpanLocation(read, location2);
            if (read.getProperPairFlag()) {
                if (spanLocation1 || spanLocation2) {
                    output = outputProperPairs;
                }
            } else {
                if (spanLocation1 || spanLocation2) {
                    interestingReads++;
                    output = true;
                } else if (doesSpanLocationUnclipped(read, location1) || doesSpanLocationUnclipped(read, location2)) {
                    interestingReads++;
                    output = true;
                }
            }

            if (!output)
                continue;

            // @formatter:off
            outputMap.put(read.getAlignmentStart(),
                    String.join("\t",
                        read.getReadName(),
                        Integer.toString(read.getInferredInsertSize()),
                        getOrientationString(read),
                        read.getReferenceName(),
                        Integer.toString(read.getAlignmentStart()),
                        Integer.toString(read.getMappingQuality()),
                        getFlagString(read.getSAMFlags()),
                        read.getCigarString(),
                        read.getReadString(),
                        read.getBaseQualityString(),
                        read.getMateReferenceName(),
                        Integer.toString(read.getMateAlignmentStart())
            ));
            // @formatter:on
        }

        for (String entry : outputMap.values()) {
            System.out.println(entry);
        }

        System.out.println();
        System.out.println("-STATS-");
        System.out.println("INTERESTING_READS\t" + interestingReads);
        System.out.println("NORMAL_READS_BP1\t" + normalReadsInLocation1);
        System.out.println("NORMAL_READS_BP2\t" + normalReadsInLocation2);
        System.out.println("ZERO_MAPQ\t" + zeroQualityMappings);
        System.out.println("PAIRS_MISSING_READ\t" + readMap.size());

    }
}
