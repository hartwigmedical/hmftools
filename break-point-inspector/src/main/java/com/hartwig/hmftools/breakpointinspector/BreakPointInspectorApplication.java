package com.hartwig.hmftools.breakpointinspector;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

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
    private static final String INCLUDE_FILTERED = "filtered";

    private enum ReadCategory {
        UNSET,
        NORMAL,
        SPAN,
        SINGLE,
        DIFF_CHROMOSOME,
        UNMAPPED,
        SECONDARY,
        CHIMERIC
    }

    private enum LocationInfo {
        FILTERED,
        PROXIMITY,
        INTERSECT,
        STRADDLE,
        CLIP,
    }

    private enum Region {
        UNSET,
        BP1,
        BP2
    }

    private static class ReadInfo {
        public SAMRecord Read = null;
        public BreakPointInspectorApplication.Region Region = BreakPointInspectorApplication.Region.UNSET;
        public ReadCategory Category = ReadCategory.UNSET;
        public LocationInfo Location = LocationInfo.FILTERED;
    }

    private static class ReadCollection {
        public ArrayList<ReadInfo> Infos = new ArrayList<>();
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(BAM_PATH, true, "input BAM");
        options.addOption(BREAK_POINT, true, "position of break point in chrX:123456 format");
        options.addOption(PROXIMITY, true, "base distance around breakpoint");
        options.addOption(SV_LEN, true, "length of the SV to inspect");
        options.addOption(INCLUDE_PROPER, false, "include proper reads in output");
        options.addOption(INCLUDE_FILTERED, false, "included filtered reads in output");
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
                    return "INWARDS";
                case RF:
                    return "OUTWARDS";
                case TANDEM:
                    return "TANDEM";
            }
        }
        return "UNALIGNED";
    }

    private static class OrientationStats {
        public int Inwards = 0;
        public int Outwards = 0;
        public int Tandem = 0;

        public int Total() {
            return Inwards + Outwards + Tandem;
        }
    }

    private static OrientationStats[] BP1_STATS = { new OrientationStats(), new OrientationStats() };
    private static OrientationStats[] BP2_STATS = { new OrientationStats(), new OrientationStats() };
    private static OrientationStats FILTERED = new OrientationStats();

    private static void applyToStats(final ReadInfo info) {
        OrientationStats stats = FILTERED;

        // pick breakpoint
        OrientationStats[] array = { FILTERED, FILTERED };
        if (info.Region == Region.BP1) {
            array = BP1_STATS;
        } else if (info.Region == Region.BP2) {
            array = BP2_STATS;
        }
        // pick type of location
        switch (info.Location) {
            case PROXIMITY:
                stats = array[0];
                break;
            case CLIP:
                stats = array[1];
                break;
        }
        // increment orientation counter
        if (isOrientable(info.Read)) {
            final SamPairUtil.PairOrientation orientation = SamPairUtil.getPairOrientation(info.Read);
            switch (orientation) {
                case FR:
                    stats.Inwards++;
                    break;
                case RF:
                    stats.Outwards++;
                    break;
                case TANDEM:
                    stats.Tandem++;
                    break;
            }
        }
    }

    private static String toString(final OrientationStats stats) {
        List<Integer> list = Arrays.asList(stats.Inwards, stats.Outwards, stats.Tandem, stats.Total());
        return list.stream().map(Object::toString).collect(Collectors.joining("\t"));
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

    public static void main(final String... args) throws ParseException, IOException {

        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        // grab arguments
        final String bamPath = cmd.getOptionValue(BAM_PATH);
        final String breakPoint = cmd.getOptionValue(BREAK_POINT);
        final int range = Integer.parseInt(cmd.getOptionValue(PROXIMITY, "500"));
        final int svLen = Integer.parseInt(cmd.getOptionValue(SV_LEN, "0"));
        final boolean outputProperPairs = cmd.hasOption(INCLUDE_PROPER);
        final boolean outputFilteredReads = cmd.hasOption(INCLUDE_FILTERED);

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
        int refIndex = reader.getFileHeader().getSequenceIndex(chromosome);
        if (refIndex < 0 && !chromosome.startsWith("chr"))
            refIndex = reader.getFileHeader().getSequenceIndex("chr" + chromosome);
        if (refIndex < 0) {
            System.out.println("Could not find chromosome in file");
            System.exit(1);
        }
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
            final ReadInfo info = new ReadInfo();

            info.Read = read;
            collection.Infos.add(info);

            // if unmapped there's nothing to do
            if (read.getReadUnmappedFlag()) {
                info.Region = Region.UNSET;
                info.Category = ReadCategory.UNMAPPED;
                info.Location = LocationInfo.FILTERED;
                continue;
            }

            // determine the region
            final boolean closerToRegion1 =
                    Math.abs(read.getAlignmentStart() - location1) < Math.abs(read.getAlignmentStart() - location2);
            info.Region = closerToRegion1 ? Region.BP1 : Region.BP2;
            final int location = closerToRegion1 ? location1 : location2;

            if (read.getProperPairFlag()) {
                info.Category = ReadCategory.NORMAL;
                // only care about normals in Location1
                if (readIntersectsLocation(read, location)) {
                    info.Location = LocationInfo.INTERSECT;
                } else if (pairStraddlesLocation(read, location)) {
                    info.Location = LocationInfo.STRADDLE;
                } else {
                    info.Location = LocationInfo.FILTERED;
                }
            } else if (read.getMateUnmappedFlag()) {
                info.Category = ReadCategory.SINGLE;
                info.Location = LocationInfo.FILTERED;
            } else if (read.getReferenceIndex() != read.getMateReferenceIndex()) {
                info.Category = ReadCategory.DIFF_CHROMOSOME;
                info.Location = LocationInfo.FILTERED;
            } else {

                // determine type
                if (read.getSupplementaryAlignmentFlag())
                    info.Category = ReadCategory.CHIMERIC;
                else if (read.getNotPrimaryAlignmentFlag()) {
                    info.Category = ReadCategory.SECONDARY;
                    info.Location = LocationInfo.FILTERED;
                } else if (read.getReadPairedFlag())
                    info.Category = ReadCategory.SPAN;

                // determine classification

                final boolean clipped = read.getUnclippedStart() != read.getAlignmentStart()
                        || read.getUnclippedEnd() != read.getAlignmentEnd();
                if (clipped && read.getAlignmentStart() - 1 == location) {
                    info.Location = LocationInfo.CLIP;
                } else if (clipped && read.getAlignmentEnd() == location) {
                    info.Location = LocationInfo.CLIP;
                } else if (readIntersectsLocation(read, location)) {
                    info.Location = LocationInfo.CLIP; // map this to clip intentionally
                } else {
                    info.Location = LocationInfo.PROXIMITY;
                }
            }
        }

        // print  header
        // @formatter:off
        System.out.println(
                String.join("\t",
                        "READ_NAME",
                        "REGION",
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
                        "MATE_ALIGNMENT_START"
                ));
        // @formatter:on

        // manta style categorisation
        int negativePairs = 0;
        int positivePairs = 0;
        int negativeSplitReads = 0;
        int positiveSplitReads = 0;
        int readsMissingMate = 0;

        final List<Integer> insertSizes = new ArrayList<>();
        for (final ReadCollection reads : readMap.values()) {

            // make sure we at least have the primary pairs
            final List<ReadInfo> pair = reads.Infos.stream().filter(i -> !i.Read.getNotPrimaryAlignmentFlag()).collect(
                    Collectors.toList());
            if (pair.size() != 2) {
                readsMissingMate++;
                continue;
            }

            // update the stats
            final boolean normal = pair.get(0).Category == ReadCategory.NORMAL;
            if (pair.get(0).Location == LocationInfo.FILTERED || pair.get(1).Location == LocationInfo.FILTERED) {
                // we won't count these at all
            } else if (normal) {
                if (pair.get(0).Region == Region.BP1) {
                    if (pair.get(0).Location == LocationInfo.STRADDLE
                            && pair.get(1).Location == LocationInfo.STRADDLE) {
                        negativePairs++;
                    } else if (pair.get(0).Location == LocationInfo.INTERSECT
                            || pair.get(1).Location == LocationInfo.INTERSECT) {
                        negativePairs++;
                        negativeSplitReads++;
                    }
                }
            } else {
                // we should check by orientation too
                if (pair.get(0).Region != pair.get(1).Region)
                    positivePairs++;

                final List<LocationInfo> clipOrIntersect = Arrays.asList(LocationInfo.CLIP, LocationInfo.INTERSECT);
                if (clipOrIntersect.contains(pair.get(0).Location) || clipOrIntersect.contains(pair.get(1).Location))
                    positiveSplitReads++;
            }

            // output the reads
            for (final ReadInfo info : reads.Infos) {

                // determine if we want to output this read
                if (info.Category == ReadCategory.NORMAL) {
                    switch (info.Location) {
                        case STRADDLE:
                        case INTERSECT:
                            if (!outputProperPairs) {
                                continue;
                            }
                            break;
                        default:
                            continue;
                    }
                } else if (info.Location == LocationInfo.FILTERED) {
                    if (!outputFilteredReads)
                        continue;
                } else {
                    applyToStats(info);
                }

                // @formatter:off
                final SAMRecord r = info.Read;
                String output = String.join("\t",
                        r.getReadName(),
                        info.Region.toString(),
                        String.join(",", info.Category.toString(), info.Location.toString()),
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
        System.out.println("PR\tSR");
        System.out.println(
                Integer.toString(negativePairs) + ":" + Integer.toString(positivePairs) + "\t" + Integer.toString(
                        negativeSplitReads) + ":" + Integer.toString(positiveSplitReads));

        // @formatter:off
        System.out.println(String.join("\t", "-READS-", "INWARDS", "OUTWARDS", "TANDEM", "TOTAL"));
        System.out.println(String.join("\t", "BP1_PROXIMITY", toString(BP1_STATS[0])));
        System.out.println(String.join("\t", "BP1_CLIPPED", toString(BP1_STATS[1])));
        System.out.println(String.join("\t", "BP2_PROXIMITY", toString(BP2_STATS[0])));
        System.out.println(String.join("\t", "BP2_CLIPPED", toString(BP2_STATS[1])));
        //System.out.println(String.join("\t", "FILTERED"));
        //System.out.println(String.join("\t", "TOTAL"));

        System.out.println();
        System.out.println("READS_MISSING_MATE\t" + readsMissingMate);
    }
}
