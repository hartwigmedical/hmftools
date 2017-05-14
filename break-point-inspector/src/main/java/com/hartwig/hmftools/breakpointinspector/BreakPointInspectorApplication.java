package com.hartwig.hmftools.breakpointinspector;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
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
    private static final String PROXIMITY = "proximity";
    private static final String SV_LEN = "svlen";
    private static final String INCLUDE_PROPER = "proper";

    private static class BreakPointStats {
        public int NormalReads = 0;
        public int InterestingReadIntersecting = 0;
        public int InterestingReadProximity = 0;
        public int SoftClippedIntersecting = 0;
        public int MatedToDifferentChromosome = 0;
        public int UnmappedMate = 0;
    }

    private static class OrientationStats {
        public int InnieCount = 0;
        public int OutieCount = 0;
        public int TandemCount = 0;
        public int UnorientedCount = 0;
    }

    private static class PairedRead {
        public SAMRecord First = null;
        public SAMRecord Second = null;
        public boolean FirstInteresting = false;
        public boolean SecondInteresting = false;
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(BAM_PATH, true, "input BAM");
        options.addOption(BREAK_POINT, true, "position of break point in chrX:123456 format");
        options.addOption(PROXIMITY, true, "base distance around breakpoint");
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

    private static void printStats(final BreakPointStats stats) {
        System.out.println("INTERESTING_INTERSECTING\t" + stats.InterestingReadIntersecting);
        System.out.println("INTERESTING_PROXIMITY\t" + stats.InterestingReadProximity);
        System.out.println("SOFT_CLIPPED_READS_INTERSECTING\t" + stats.SoftClippedIntersecting);
        System.out.println("MATED_TO_DIFF_CHROMOSOME\t" + stats.MatedToDifferentChromosome);
        System.out.println("UNMAPPED_MATE\t" + stats.UnmappedMate);
        System.out.println("NORMAL_INTERSECTING\t" + stats.NormalReads);
    }

    public static void main(final String... args) throws ParseException, IOException {

        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        // grab arguments
        final String bamPath = cmd.getOptionValue(BAM_PATH);
        final String breakPoint = cmd.getOptionValue(BREAK_POINT);
        final int range = Integer.parseInt(cmd.getOptionValue(PROXIMITY, "500"));
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

        final Map<String, PairedRead> readMap = new Hashtable<>();
        final Map<Integer, String> outputMap = new TreeMap<>();
        final List<SAMRecord> evidenceReads = new ArrayList<>();

        BreakPointStats statsLocation1 = new BreakPointStats();
        BreakPointStats statsLocation2 = new BreakPointStats();
        OrientationStats orientationStats = new OrientationStats();

        // execute query and parse the results
        final SAMRecordIterator results = reader.query(queryIntervals, false);
        while (results.hasNext()) {
            final SAMRecord read = results.next();

            // we classify all reads into location 1 or 2
            final boolean spanLocation1 = doesSpanLocation(read, location1);
            final boolean spanLocation2 = doesSpanLocation(read, location2);

            final boolean isProper = read.getProperPairFlag(); // i.e. normal insert size, right orientation etc.

            boolean evidence = false;
            if (isProper) {
                if (spanLocation1) {
                    statsLocation1.NormalReads++;
                }
                if (spanLocation2) {
                    statsLocation2.NormalReads++;
                }
            } else {
                // filter for the reads we want
                final boolean differentChromosomes = !read.getReferenceName().equals(read.getMateReferenceName());
                final boolean oneOfPairUnmapped = read.getReadUnmappedFlag() ^ read.getMateUnmappedFlag();
                final boolean isSecondaryOrSupplementary =
                        read.getNotPrimaryAlignmentFlag() || read.getSupplementaryAlignmentFlag(); // TODO: ???
                final boolean closerToLocation1 = Math.abs(read.getAlignmentStart() - location1) < Math.abs(
                        read.getAlignmentStart() - location2);

                if (differentChromosomes) {
                    if (closerToLocation1)
                        statsLocation1.MatedToDifferentChromosome++;
                    else
                        statsLocation2.MatedToDifferentChromosome++;
                    evidence = true;
                } else if (oneOfPairUnmapped) {
                    if (closerToLocation1)
                        statsLocation1.UnmappedMate++;
                    else
                        statsLocation2.UnmappedMate++;
                    evidence = true;
                } else if (isSecondaryOrSupplementary) {
                    // TODO: classify
                    evidence = true;
                } else if (spanLocation1) {
                    statsLocation1.InterestingReadIntersecting++;
                    evidence = true;
                } else if (spanLocation2) {
                    statsLocation2.InterestingReadIntersecting++;
                    evidence = true;
                } else if (doesSpanLocationUnclipped(read, location1)) {
                    statsLocation1.SoftClippedIntersecting++;
                    evidence = true;
                } else if (doesSpanLocationUnclipped(read, location2)) {
                    statsLocation2.SoftClippedIntersecting++;
                    evidence = true;
                } else if (closerToLocation1) {
                    statsLocation1.InterestingReadProximity++;
                    evidence = true;
                } else if (!closerToLocation1) {
                    statsLocation2.InterestingReadProximity++;
                    evidence = true;
                }

                if (evidence) {
                    evidenceReads.add(read);
                    if (read.getReadPairedFlag() && !read.getReadUnmappedFlag() && !read.getMateUnmappedFlag()) {
                        final SamPairUtil.PairOrientation orientation = SamPairUtil.getPairOrientation(read);
                        switch (orientation) {
                            case FR:
                                orientationStats.InnieCount++;
                                break;
                            case RF:
                                orientationStats.OutieCount++;
                                break;
                            case TANDEM:
                                orientationStats.TandemCount++;
                                break;
                        }
                    } else {
                        orientationStats.UnorientedCount++;
                    }
                }
            }

            PairedRead pairedRead = readMap.get(read.getReadName());
            if (pairedRead != null) {
                // we've found the mate
                final SAMRecord firstRead = pairedRead.First;
                assert firstRead.getReadName().equals(read.getReadName());
                assert firstRead.getMateReferenceName().equals(read.getReferenceName());
                assert firstRead.getMateAlignmentStart() == read.getAlignmentStart();
                assert firstRead.getProperPairFlag() == read.getProperPairFlag();

                pairedRead.Second = read;
                pairedRead.SecondInteresting = evidence;
            } else {
                // this is the first read we've seen of the pair
                pairedRead = new PairedRead();
                pairedRead.First = read;
                pairedRead.FirstInteresting = evidence;
                readMap.put(read.getReadName(), pairedRead);
            }
        }

        for (final PairedRead pair : readMap.values()) {

            assert pair.First != null;

            if (!outputProperPairs && !pair.FirstInteresting && !pair.SecondInteresting)
                continue;

            // @formatter:off
            String output = String.join("\t",
                    pair.First.getReadName(),
                    Integer.toString(pair.First.getInferredInsertSize()),
                    getOrientationString(pair.First),
                    pair.First.getReferenceName(),
                    Integer.toString(pair.First.getAlignmentStart()),
                    Integer.toString(pair.First.getMappingQuality()),
                    getFlagString(pair.First.getSAMFlags()),
                    pair.First.getCigarString(),
                    pair.First.getReadString(),
                    pair.First.getBaseQualityString()
            );
            if(pair.Second != null) {
                output += "\t" + String.join("\t",
                        pair.Second.getReferenceName(),
                        Integer.toString(pair.Second.getAlignmentStart()),
                        Integer.toString(pair.Second.getMappingQuality()),
                        getFlagString(pair.Second.getSAMFlags()),
                        pair.Second.getCigarString(),
                        pair.Second.getReadString(),
                        pair.Second.getBaseQualityString()
                );
            } else {
                output += "\t" + String.join("\t",
                        pair.First.getMateReferenceName(),
                        Integer.toString(pair.First.getMateAlignmentStart())
                        );
            }
            // @formatter:on

            outputMap.put(pair.First.getAlignmentStart(), output);
        }

        // print  header
        // @formatter:off
        System.out.println(
                String.join("\t",
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
        // @formatter:on

        // print entries sorted by chromosome position
        for (final String entry : outputMap.values())
            System.out.println(entry);

        System.out.println();
        System.out.println("-BP1_STATS-");
        printStats(statsLocation1);
        System.out.println("-BP2_STATS-");
        printStats(statsLocation2);

        System.out.println("-ORIENTATION-");
        System.out.println("INNIE_COUNT\t" + orientationStats.InnieCount);
        System.out.println("OUTIE_COUNT\t" + orientationStats.OutieCount);
        System.out.println("TANDEM_COUNT\t" + orientationStats.TandemCount);
        System.out.println("UNORIENTED_COUNTED\t" + orientationStats.UnorientedCount);
    }
}
