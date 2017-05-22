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
import java.util.stream.Stream;

import htsjdk.samtools.*;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import org.jetbrains.annotations.NotNull;

public class BreakPointInspectorApplication {

    private static final String REF_PATH = "ref";
    private static final String TUMOR_PATH = "tumor";
    private static final String BREAK_POINT1 = "bp1";
    private static final String BREAK_POINT2 = "bp2";
    private static final String PROXIMITY = "proximity";
    private static final String SV_LEN = "svlen";

    private static final int MANTA_REQ_PAIR_MIN = 50;
    private static final int MANTA_REQ_SPLIT_MIN = 15;

    private enum ReadCategory {
        UNSET,
        NORMAL,
        SPAN,
        MATE_UNMAPPED,
        UNMAPPED,
        SECONDARY,
        CHIMERIC
    }

    private enum Overlap {
        FILTERED,
        PROXIMITY,
        INTERSECT,
        STRADDLE,
        CLIP,
    }

    private enum Region {
        OTHER,
        BP1,
        BP2
    }

    private enum ClipSide {
        NONE,
        LEFT_CLIP,
        RIGHT_CLIP
    }

    private static class ClipInfo {
        public ClipSide Side = ClipSide.NONE;
        public int Length = 0;
    }

    private static class ReadInfo {
        public SAMRecord Read = null;
        public BreakPointInspectorApplication.Region Region = BreakPointInspectorApplication.Region.OTHER;
        public ReadCategory Category = ReadCategory.UNSET;
        public Overlap Location = Overlap.FILTERED;
        public ClipInfo Clipping = new ClipInfo();
    }

    private static class ReadCollection {
        public ArrayList<ReadInfo> Infos = new ArrayList<>();
    }

    private static class Location {
        public int ReferenceIndex = -1;
        public int Position = -1;

        public static Location parseLocationString(final String location, final SAMSequenceDictionary dictionary)
                throws RuntimeException {
            final Location result = new Location();

            final String[] split = location.split(":");
            if (split.length != 2)
                throw new RuntimeException(location + " is not a valid location string");

            final String chromosome = split[0];
            try {
                result.Position = Integer.parseInt(split[1]);
            } catch (NumberFormatException e) {
                throw new RuntimeException(location + " is not a valid location string");
            }

            // query the position
            result.ReferenceIndex = dictionary.getSequenceIndex(chromosome);
            if (result.ReferenceIndex < 0) {
                if (!chromosome.startsWith("chr"))
                    result.ReferenceIndex = dictionary.getSequenceIndex("chr" + chromosome);
                else
                    result.ReferenceIndex = dictionary.getSequenceIndex(chromosome.substring(3));
            }
            if (result.ReferenceIndex < 0) {
                throw new RuntimeException(chromosome + " is not in the BAM");
            }

            return result;
        }
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(REF_PATH, true, "the Reference BAM");
        options.addOption(TUMOR_PATH, true, "the Tumor BAM");
        options.addOption(BREAK_POINT1, true, "position of first break point in chrX:123456 format");
        options.addOption(BREAK_POINT2, true, "position of second break point in chrX:123456 format (optional)");
        options.addOption(PROXIMITY, true, "base distance around breakpoint");
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
        ArrayList<String> names = new ArrayList<>();
        for (final SAMFlag flag : flags) {
            names.add(flag.name());
        }
        return String.join("|", names);
    }

    private static boolean isOrientable(final SAMRecord read) {
        return read.getReadPairedFlag() && !read.getReadUnmappedFlag() && !read.getMateUnmappedFlag();
    }

    private static void printHelpAndExit(final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Break-Point-Inspector", "Inspect structural variants", options, "", true);
        System.exit(1);
    }

    private static boolean readIntersectsLocation(final SAMRecord read, final Location location) {
        assert !read.getReadUnmappedFlag();
        return read.getReferenceIndex() == location.ReferenceIndex && read.getAlignmentStart() <= location.Position
                && read.getAlignmentEnd() >= location.Position;
    }

    private static boolean pairStraddlesLocation(final SAMRecord read, final Location location) {
        assert !read.getReadUnmappedFlag() && !read.getMateUnmappedFlag();
        if (read.getInferredInsertSize() > 0)
            return read.getAlignmentStart() <= location.Position
                    && (read.getAlignmentStart() + read.getInferredInsertSize()) >= location.Position;
        else
            return read.getMateAlignmentStart() <= location.Position
                    && (read.getMateAlignmentStart() - read.getInferredInsertSize()) >= location.Position;
    }

    private static ClipInfo getClipInfo(final SAMRecord read) {
        final ClipInfo result = new ClipInfo();
        final Cigar cigar = read.getCigar();
        switch (cigar.getFirstCigarElement().getOperator()) {
            case S:
            case H:
                result.Side = ClipSide.LEFT_CLIP;
                result.Length = cigar.getFirstCigarElement().getLength();
                return result;
        }
        switch (cigar.getLastCigarElement().getOperator()) {
            case S:
            case H:
                result.Side = ClipSide.RIGHT_CLIP;
                result.Length = cigar.getLastCigarElement().getLength();
                return result;
        }
        return result;
    }

    private static Region determineRegion(final SAMRecord read, final Location location1, final Location location2) {
        if (location1.ReferenceIndex != location2.ReferenceIndex) {
            if (read.getReferenceIndex() == location1.ReferenceIndex)
                return Region.BP1;
            else if (read.getReferenceIndex() == location2.ReferenceIndex)
                return Region.BP2;
        } else if (read.getReferenceIndex() == location1.ReferenceIndex) {
            if (Math.abs(read.getAlignmentStart() - location1.Position) < Math.abs(
                    read.getAlignmentStart() - location2.Position))
                return Region.BP1;
            else
                return Region.BP2;
        }
        return Region.OTHER;
    }

    private static class ClassifiedReadResults {
        public Map<String, ReadCollection> ReadMap = new Hashtable<>();
    }

    private static class BreakPointStats {
        public int PR_Only_Normal = 0;
        public int PR_SR_Normal = 0;
        public int PR_Only_Support = 0;
        public int PR_SR_Support = 0;
        public int SR_Only_Support = 0;
        public int Unmapped_Mate = 0;
        public int Diff_Variant = 0;

        public static List<String> GetHeader() {
            return Arrays.asList("PR_ONLY_NORMAL", "PR_SR_NORMAL", "PR_ONLY_SUPPORT", "PR_SR_SUPPORT",
                    "SR_ONLY_SUPPORT", "UNMAPPED_MATE", "DIFF_VARIANT");
        }

        public List<Integer> GetData() {
            return Arrays.asList(PR_Only_Normal, PR_SR_Normal, PR_Only_Support, PR_SR_Support, SR_Only_Support,
                    Unmapped_Mate, Diff_Variant);
        }
    }

    private static class OrientationStats {
        public int InnieCount = 0;
        public int OutieCount = 0;
        public int TandemCount = 0;

        public List<Integer> GetData() {
            return Arrays.asList(InnieCount, OutieCount, TandemCount);
        }
    }

    private static class SampleStats {
        public BreakPointStats BP1_Stats = new BreakPointStats();
        public BreakPointStats BP2_Stats = new BreakPointStats();
        public OrientationStats Orientation = new OrientationStats();

        public BreakPointStats Get(final Region r) {
            switch (r) {
                case BP1:
                    return BP1_Stats;
                case BP2:
                    return BP2_Stats;
                default:
                    throw new RuntimeException("invalid stats");
            }
        }

        public void Print(final String prefix) {
            final String header = Stream.concat(BreakPointStats.GetHeader().stream().map(h -> "BP1_" + h),
                    BreakPointStats.GetHeader().stream().map(h -> "BP2_" + h)).map(h -> prefix + h).collect(
                    Collectors.joining("\t"));
            final String data = Stream.concat(BP1_Stats.GetData().stream(), BP2_Stats.GetData().stream()).map(
                    i -> Integer.toString(i)).collect(Collectors.joining("\t"));
            System.out.println(header);
            System.out.println(data);
        }
    }

    private static SampleStats calculateStats(final ClassifiedReadResults queryResult) {

        final SampleStats result = new SampleStats();

        for (final ReadCollection reads : queryResult.ReadMap.values()) {

            final List<ReadInfo> pair = reads.Infos;
            final ReadInfo p0 = pair.get(0);

            // handle odd sizes
            if (pair.size() < 2) {
                if (p0.Category == ReadCategory.MATE_UNMAPPED) {
                    result.Get(p0.Region).Unmapped_Mate++;
                } else if (p0.Category != ReadCategory.NORMAL || p0.Location != Overlap.FILTERED) {
                    result.Get(p0.Region).Diff_Variant++;
                }
                continue;
            } else if (pair.size() > 2) {
                // TODO: secondary / supplementary reads?
                System.err.println("condition pair.size() > 3 not implemented: " + p0.Read.getReadName());
                continue;
            }

            final ReadInfo p1 = pair.get(1);

            // possible two paired but unmapped reads?
            if (p0.Region == Region.OTHER && p1.Region == Region.OTHER) {
                continue;
            } else if (p1.Region == Region.OTHER) { // must be unmapped mate
                result.Get(p0.Region).Unmapped_Mate++;
                continue;
            } else if (p0.Region == Region.OTHER) { // must be unmapped mate
                result.Get(p1.Region).Unmapped_Mate++;
                continue;
            }

            final boolean pairPossible = pair.stream().allMatch(p -> p.Clipping.Length <= MANTA_REQ_PAIR_MIN);

            if (p0.Region != p1.Region) {
                // supports the break point
                boolean interesting = false;
                for (final ReadInfo r : pair) {
                    final boolean splitEvidence =
                            r.Location == Overlap.CLIP && r.Clipping.Length >= MANTA_REQ_SPLIT_MIN;
                    if (pairPossible && splitEvidence) {
                        result.Get(r.Region).PR_SR_Support++;
                    } else if (pairPossible) {
                        result.Get(r.Region).PR_Only_Support++;
                    } else if (splitEvidence) {
                        result.Get(r.Region).SR_Only_Support++;
                    }

                    interesting |= pairPossible || splitEvidence;
                }

                if (interesting) {
                    switch (SamPairUtil.getPairOrientation(p0.Read)) {
                        case FR:
                            result.Orientation.InnieCount++;
                            break;
                        case RF:
                            result.Orientation.OutieCount++;
                            break;
                        case TANDEM:
                            result.Orientation.TandemCount++;
                            break;
                    }
                }
            } else {
                final BreakPointStats stats = p0.Region == Region.BP1 ? result.BP1_Stats : result.BP2_Stats;
                final boolean splitEvidence = pair.stream().anyMatch(
                        p -> p.Location == Overlap.CLIP && p.Clipping.Length >= MANTA_REQ_SPLIT_MIN);
                if (splitEvidence) {
                    stats.SR_Only_Support++;
                } else if (p0.Location == Overlap.STRADDLE && p1.Location == Overlap.STRADDLE) {
                    stats.PR_Only_Normal++;
                } else if (p0.Location == Overlap.INTERSECT || p1.Location == Overlap.INTERSECT) {
                    stats.PR_SR_Normal++;
                    // TODO: does the intersection have to occur within a certain bound of read?
                }
            }
        }

        return result;
    }

    private static ClassifiedReadResults performQueryAndClassify(final SamReader reader,
            final QueryInterval[] intervals, final Location bp1, final Location bp2) {

        final ClassifiedReadResults result = new ClassifiedReadResults();

        // execute query and parse the results
        final SAMRecordIterator results = reader.query(intervals, false);
        while (results.hasNext()) {

            final SAMRecord read = results.next();
            final ReadCollection collection = result.ReadMap.computeIfAbsent(read.getReadName(),
                    k -> new ReadCollection());
            final ReadInfo info = new ReadInfo();

            info.Read = read;
            collection.Infos.add(info);

            // if unmapped there's nothing to do
            if (read.getReadUnmappedFlag()) {
                info.Region = Region.OTHER;
                info.Category = ReadCategory.UNMAPPED;
                info.Location = Overlap.FILTERED;
                continue;
            }

            info.Clipping = getClipInfo(read);
            info.Region = determineRegion(read, bp1, bp2);
            if (info.Region == Region.OTHER) {
                throw new RuntimeException("read from unexpected region");
            }

            final Location bp = info.Region == Region.BP1 ? bp1 : bp2;
            final boolean clipped = info.Clipping.Side != ClipSide.NONE;
            final boolean normal = read.getProperPairFlag();

            if (read.getProperPairFlag()) {
                info.Category = ReadCategory.NORMAL;

                if (!clipped) {
                    if (readIntersectsLocation(read, bp)) {
                        info.Location = Overlap.INTERSECT;
                    } else if (pairStraddlesLocation(read, bp)) {
                        info.Location = Overlap.STRADDLE;
                    } else {
                        info.Location = Overlap.FILTERED;
                    }
                    continue;
                }

            } else {
                // determine type
                if (read.getMateUnmappedFlag()) {
                    info.Category = ReadCategory.MATE_UNMAPPED;
                } else if (read.getSupplementaryAlignmentFlag()) {
                    info.Category = ReadCategory.CHIMERIC;
                } else if (read.getNotPrimaryAlignmentFlag()) {
                    info.Category = ReadCategory.SECONDARY;
                } else if (read.getReadPairedFlag())
                    info.Category = ReadCategory.SPAN;
            }

            // determine classification
            if (clipped && read.getAlignmentStart() - 1 == bp.Position) {
                info.Location = Overlap.CLIP;
            } else if (clipped && read.getAlignmentEnd() == bp.Position) {
                info.Location = Overlap.CLIP;
            } else if (readIntersectsLocation(read, bp)) {
                info.Location = Overlap.CLIP; // map this to clip intentionally
            } else {
                info.Location = Overlap.PROXIMITY;
            }

            // if normal read and we don't clip, then it's not interesting
            if (normal && info.Location == Overlap.PROXIMITY)
                info.Location = Overlap.FILTERED;
        }

        return result;
    }

    public static void main(final String... args) throws ParseException, IOException {

        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        // grab arguments
        final String refBAM = cmd.getOptionValue(REF_PATH);
        final String tumorBAM = cmd.getOptionValue(TUMOR_PATH);
        final String bp1String = cmd.getOptionValue(BREAK_POINT1);
        final String bp2String = cmd.getOptionValue(BREAK_POINT2);
        final int range = Integer.parseInt(cmd.getOptionValue(PROXIMITY, "500"));
        final int svLen = Integer.parseInt(cmd.getOptionValue(SV_LEN, "0"));

        if (refBAM == null || tumorBAM == null || bp1String == null)
            printHelpAndExit(options);

        // load the files
        final File tumorFile = new File(tumorBAM);
        final SamReader tumorReader = SamReaderFactory.makeDefault().open(tumorFile);
        final File refFile = new File(refBAM);
        final SamReader refReader = SamReaderFactory.makeDefault().open(refFile);

        // parse the location strings
        final Location location1 = Location.parseLocationString(bp1String,
                tumorReader.getFileHeader().getSequenceDictionary());
        final Location location2;
        if (bp2String != null) {
            location2 = Location.parseLocationString(bp2String, tumorReader.getFileHeader().getSequenceDictionary());
        } else if (svLen > 0) {
            location2 = new Location();
            location2.ReferenceIndex = location1.ReferenceIndex;
            location2.Position = location1.Position + svLen;
        } else {
            printHelpAndExit(options);
            return;
        }

        // work out the query intervals
        QueryInterval[] queryIntervals = {
                new QueryInterval(location1.ReferenceIndex, Math.max(0, location1.Position - range),
                        location1.Position + range),
                new QueryInterval(location2.ReferenceIndex, Math.max(0, location2.Position - range),
                        location2.Position + range) };
        queryIntervals = QueryInterval.optimizeIntervals(queryIntervals);

        // begin processing

        final ClassifiedReadResults refResult = performQueryAndClassify(refReader, queryIntervals, location1,
                location2);
        final SampleStats refStats = calculateStats(refResult);
        refStats.Print("REF_");

        final ClassifiedReadResults tumorResult = performQueryAndClassify(tumorReader, queryIntervals, location1,
                location2);
        final SampleStats tumorStats = calculateStats(tumorResult);
        tumorStats.Print("TUMOR_");

        System.out.println("INNIE\tOUTIE\tTANDEM");
        System.out.println(tumorStats.Orientation.GetData().stream().map(i -> Integer.toString(i)).collect(
                Collectors.joining("\t")));
    }
}
