package com.hartwig.hmftools.breakpointinspector;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.TreeMultimap;

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

    private enum ReadType {
        UNSET,
        NORM,
        SPAN,
        SINGLE,
        DIFF_CHROMOSOME,
        UNMAPPED,
        SECONDARY,
        CHIMERIC
    }

    private enum ReadClassification {
        FILTERED,
        PROXIMITY,
        INTERSECT,
        STRADDLE,
        CLIPPED
    }

    private static class ReadInfo {
        public SAMRecord Read = null;
        public ReadType Type = ReadType.UNSET;
        public ReadClassification Classification = ReadClassification.FILTERED;
    }

    private static class ReadCollection {
        public ArrayList<ReadInfo> Reads = new ArrayList<>();
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
        ArrayList<String> names = new ArrayList<>();
        for (final SAMFlag flag : flags) {
            names.add(flag.name());
        }
        return String.join("|", names);
    }

    private static boolean isOrientable(final SAMRecord read) {
        return read.getReadPairedFlag() && !read.getReadUnmappedFlag() && !read.getMateUnmappedFlag();
    }

    @NotNull
    private static String getOrientationString(final SAMRecord read) {
        if (isOrientable(read)) {
            switch (SamPairUtil.getPairOrientation(read)) {
                case FR:
                    return "NORM";
                case RF:
                    return "OUTIE";
                case TANDEM:
                    return "TANDEM";
            }
        }
        return "UNALIGNED";
    }

    private static void printHelpAndExit(final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Break-Point-Inspector", "Retrieve reads from an indexed BAM", options, "", true);
        System.exit(1);
    }

    private static boolean readIntersectsLocation(final SAMRecord read, final int location) {
        assert !read.getReadUnmappedFlag();
        return read.getAlignmentStart() <= location && read.getAlignmentEnd() >= location;
    }

    private static boolean pairStraddlesLocation(final SAMRecord read, final int location) {
        assert !read.getReadUnmappedFlag() && !read.getMateUnmappedFlag();
        if (read.getInferredInsertSize() > 0)
            return read.getAlignmentStart() <= location
                    && (read.getAlignmentStart() + read.getInferredInsertSize()) >= location;
        else
            return read.getMateAlignmentStart() <= location
                    && (read.getMateAlignmentStart() - read.getInferredInsertSize()) >= location;
    }

    public static void main(final String... args) throws ParseException, IOException, RuntimeException {

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
        final int refIndex = reader.getFileHeader().getSequenceDictionary().getSequenceIndex(chromosome);
        QueryInterval[] queryIntervals = {
                new QueryInterval(refIndex, Math.max(0, location1 - range), location1 + range),
                new QueryInterval(refIndex, Math.max(0, location2 - range), location2 + range) };
        queryIntervals = QueryInterval.optimizeIntervals(queryIntervals);

        final Map<String, ReadCollection> readMap = new Hashtable<>();
        final TreeMultimap<Integer, String> outputMap = TreeMultimap.create();

        // execute query and parse the results
        final SAMRecordIterator results = reader.query(queryIntervals, false);
        while (results.hasNext()) {
            final SAMRecord read = results.next();
            final ReadCollection collection = readMap.computeIfAbsent(read.getReadName(), k -> new ReadCollection());

            final ReadInfo readInfo = new ReadInfo();
            readInfo.Read = read;
            collection.Reads.add(readInfo);

            final boolean isProper = read.getProperPairFlag(); // i.e. normal insert size, right orientation etc.

            if (isProper) {
                readInfo.Type = ReadType.NORM;
                // only care about normals in Location1
                if (readIntersectsLocation(read, location1)) {
                    readInfo.Classification = ReadClassification.INTERSECT;
                } else if (pairStraddlesLocation(read, location1)) {
                    readInfo.Classification = ReadClassification.STRADDLE;
                } else {
                    readInfo.Classification = ReadClassification.FILTERED;
                }
            } else if (read.getReadUnmappedFlag()) {
                readInfo.Type = ReadType.UNMAPPED;
                readInfo.Classification = ReadClassification.FILTERED;
            } else if (read.getMateUnmappedFlag()) {
                readInfo.Type = ReadType.SINGLE;
                readInfo.Classification = ReadClassification.FILTERED;
            } else if (read.getReferenceIndex() != read.getMateReferenceIndex()) {
                readInfo.Type = ReadType.DIFF_CHROMOSOME;
                readInfo.Classification = ReadClassification.FILTERED;
            } else {
                // determine type
                if (read.getSupplementaryAlignmentFlag())
                    readInfo.Type = ReadType.CHIMERIC;
                else if (read.getNotPrimaryAlignmentFlag())
                    readInfo.Type = ReadType.SECONDARY;
                else if (read.getReadPairedFlag())
                    readInfo.Type = ReadType.SPAN;

                // determine classification
                final boolean clipped = read.getUnclippedStart() != read.getAlignmentStart()
                        || read.getUnclippedEnd() != read.getAlignmentEnd();
                if (clipped && read.getAlignmentStart() - 1 == location1) {
                    readInfo.Classification = ReadClassification.CLIPPED;
                } else if (clipped && read.getAlignmentEnd() == location1) {
                    readInfo.Classification = ReadClassification.CLIPPED;
                } else if (clipped && read.getAlignmentStart() - 1 == location2) {
                    readInfo.Classification = ReadClassification.CLIPPED;
                } else if (clipped && read.getAlignmentEnd() == location2) {
                    readInfo.Classification = ReadClassification.CLIPPED;
                } else if (readIntersectsLocation(read, location1) || readIntersectsLocation(read, location2)) {
                    readInfo.Classification = ReadClassification.INTERSECT;
                } else {
                    readInfo.Classification = ReadClassification.PROXIMITY;
                }
            }
        }

        // print  header
        // @formatter:off
        System.out.println(
                String.join("\t",
                        "READ_NAME",
                        "CLASSIFICATION",
                        "TLEN",
                        "ORIENTATION",
                        "CHROMOSOME",
                        "ALIGNMENT_START",
                        "ALIGNMENT_END",
                        "MAPQ",
                        "FLAGS",
                        "CIGAR",
                        "SEQ",
                        "QUAL",
                        "MATE_CHROMOSOME",
                        "MATE_POS"
                ));
        // @formatter:on

        int negativePairs = 0;
        int positivePairs = 0;
        int negativeSplitReads = 0;
        int positiveSplitReads = 0;

        for (final ReadCollection collection : readMap.values()) {

            for (final ReadInfo info : collection.Reads) {

                // filter reads when missing mate
                final SAMRecord r = info.Read;
                if (collection.Reads.stream().noneMatch(s -> r.getAlignmentStart() == s.Read.getMateAlignmentStart()
                        && r.getReferenceIndex() == s.Read.getMateReferenceIndex()))
                    continue;

                if (info.Type == ReadType.NORM) {
                    switch (info.Classification) {
                        case INTERSECT:
                            negativeSplitReads++;
                            break;
                        case STRADDLE:
                            if (r.getInferredInsertSize() > 0) // we want PAIRS not READS
                                negativePairs++;
                            break;
                        default:
                            continue;
                    }
                    if (!outputProperPairs)
                        continue;
                } else if (info.Classification == ReadClassification.FILTERED) {
                    continue;
                } else {
                    switch (info.Classification) {
                        case CLIPPED:
                            positiveSplitReads++;
                            if (r.getInferredInsertSize() > 0)
                                positivePairs++;
                            break;
                        case INTERSECT:
                        case PROXIMITY:
                            if (r.getInferredInsertSize() > 0) // we want PAIRS not reads
                                positivePairs++;
                            break;
                        default:
                            break;
                    }
                }

                // @formatter:off
                String output = String.join("\t",
                        r.getReadName(),
                        String.join(",", info.Type.toString(), info.Classification.toString()),
                        Integer.toString(r.getInferredInsertSize()),
                        getOrientationString(r),
                        r.getReferenceName(),
                        Integer.toString(r.getAlignmentStart()),
                        Integer.toString(r.getAlignmentEnd()),
                        Integer.toString(r.getMappingQuality()),
                        getFlagString(r.getSAMFlags()),
                        r.getCigarString(),
                        r.getReadString(),
                        r.getBaseQualityString(),
                        r.getMateReferenceName(),
                        Integer.toString(r.getMateAlignmentStart())
                );
                outputMap.put(r.getAlignmentStart(), output);
                // @formatter:on
            }
        }

        for (String s : outputMap.values())
            System.out.println(s);

        System.out.println();
        System.out.println("-EVIDENCE_READS-");
        System.out.println("PR_NEGATIVE\t" + negativePairs);
        System.out.println("PR_POSITIVE\t" + positivePairs);
        System.out.println("SR_NEGATIVE\t" + negativeSplitReads);
        System.out.println("SR_POSITIVE\t" + positiveSplitReads);
    }
}
