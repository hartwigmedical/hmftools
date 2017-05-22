package com.hartwig.hmftools.breakpointinspector;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
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

    private static final String REF_PATH = "ref";
    private static final String TUMOR_PATH = "tumor";
    private static final String BREAK_POINT1 = "bp1";
    private static final String BREAK_POINT2 = "bp2";
    private static final String PROXIMITY = "proximity";
    private static final String SV_LEN = "svlen";

    private static final String INCLUDE_PROPER = "proper";
    private static final String INCLUDE_FILTERED = "filtered";

    private enum ReadCategory {
        UNSET,
        NORMAL,
        SPAN,
        MATE_UNMAPPED,
        DIFF_CHROMOSOME,
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
        options.addOption(BREAK_POINT2, true, "position of second break point in chrX:123456 format");
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

    private static void printHelpAndExit(final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Break-Point-Inspector", "Retrieve reads from an indexed BAM", options, "", true);
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

    private static class VariantQueryResult {
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
    }

    private static class SampleStats {
        public BreakPointStats BP1_Stats = new BreakPointStats();
        public BreakPointStats BP2_Stats = new BreakPointStats();

        public BreakPointStats Get(final Region r) {
            switch (r) {
                case OTHER:
                    throw new RuntimeException("invalid stats");
                case BP1:
                    return BP1_Stats;
                case BP2:
                    return BP2_Stats;
                default:
                    throw new RuntimeException("invalid stats");
            }
        }

        public void Print() {
            // @formatter:off
            System.out.println(String.join("\t",
                    "BP1_PR_ONLY_NORMAL",
                    "BP1_PR_SR_NORMAL",
                    "BP1_PR_ONLY_SUPPORT",
                    "BP1_PR_SR_SUPPORT",
                    "BP1_SR_ONLY_SUPPORT",
                    // BP2
                    "BP2_PR_ONLY_NORMAL",
                    "BP2_PR_SR_NORMAL",
                    "BP2_PR_ONLY_SUPPORT",
                    "BP2_PR_SR_SUPPORT",
                    "BP2_SR_ONLY_SUPPORT"
            ));
            System.out.println(String.join("\t",
                    Integer.toString(BP1_Stats.PR_Only_Normal),
                    Integer.toString(BP1_Stats.PR_SR_Normal),
                    Integer.toString(BP1_Stats.PR_Only_Support),
                    Integer.toString(BP1_Stats.PR_SR_Support),
                    Integer.toString(BP1_Stats.SR_Only_Support),
                    // BP2
                    Integer.toString(BP2_Stats.PR_Only_Normal),
                    Integer.toString(BP2_Stats.PR_SR_Normal),
                    Integer.toString(BP2_Stats.PR_Only_Support),
                    Integer.toString(BP2_Stats.PR_SR_Support),
                    Integer.toString(BP2_Stats.SR_Only_Support)
            ));
            // @formatter:on
            System.out.println();
        }
    }

    private static SampleStats processResult(final VariantQueryResult queryResult) {

        final SampleStats result = new SampleStats();

        final TreeMultimap<Integer, String> outputMap = TreeMultimap.create();
        for (final ReadCollection reads : queryResult.ReadMap.values()) {

            // output the reads
            for (final ReadInfo info : reads.Infos) {
                final SAMRecord r = info.Read;
                if (info.Category == ReadCategory.NORMAL && info.Location == Overlap.FILTERED)
                    continue;

                // @formatter:off
                String output = String.join("\t",
                        "T",
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

            final List<ReadInfo> pair = reads.Infos;
            final ReadInfo p0 = pair.get(0);

            // handle odd sizes
            if (pair.size() < 2) {
                if (p0.Region == Region.OTHER) {
                    continue;
                }

                switch (p0.Category) {
                    case MATE_UNMAPPED:
                        result.Get(p0.Region).Unmapped_Mate++;
                        continue;
                    case DIFF_CHROMOSOME:
                        result.Get(p0.Region).Diff_Variant++;
                        continue;
                }
                continue;
            } else if (pair.size() > 2) {
                // TODO: secondary / supplementary reads?
                continue;
            }

            final ReadInfo p1 = pair.get(1);

            // TODO: how would this happen?
            if (p0.Region == Region.OTHER && p1.Region == Region.OTHER) {
                continue;
            } else if (p1.Region == Region.OTHER) {
                if (p1.Category == ReadCategory.DIFF_CHROMOSOME)
                    result.Get(p0.Region).Diff_Variant++;
                continue;
            } else if (p0.Region == Region.OTHER) {
                if (p0.Category == ReadCategory.DIFF_CHROMOSOME)
                    result.Get(p1.Region).Diff_Variant++;
                continue;
            }

            final boolean clipped = p0.Clipping.Length > 0 || p1.Clipping.Length > 0;
            final boolean pairPossible = p0.Clipping.Length <= 50 && p1.Clipping.Length <= 50;
            final boolean splitPossible = p0.Clipping.Length >= 15 || p1.Clipping.Length >= 15;

            final BreakPointStats stats0 = p0.Region == Region.BP1 ? result.BP1_Stats : result.BP2_Stats;
            final BreakPointStats stats1 = p1.Region == Region.BP1 ? result.BP1_Stats : result.BP2_Stats;

            if (p0.Region != p1.Region) {
                // supports the break point
                if (pairPossible && !clipped) {
                    stats0.PR_Only_Support++;
                    stats1.PR_Only_Support++;
                } else {
                    if (p0.Location == Overlap.CLIP) {
                        if (pairPossible && splitPossible) {
                            stats0.PR_SR_Support++;
                            stats1.PR_SR_Support++;
                        } else if (splitPossible) {
                            stats0.SR_Only_Support++;
                        }
                    }
                    if (p1.Location == Overlap.CLIP) {
                        if (pairPossible && splitPossible) {
                            stats0.PR_SR_Support++;
                            stats1.PR_SR_Support++;
                        } else if (splitPossible) {
                            stats1.SR_Only_Support++;
                        }
                    }

                    // TODO: should also check orientation
                }
            } else {
                // normal
                if (p0.Location == Overlap.CLIP || p1.Location == Overlap.CLIP) {
                    if (splitPossible)
                        stats0.SR_Only_Support++;
                } else if (p0.Location == Overlap.STRADDLE && p1.Location == Overlap.STRADDLE) {
                    stats0.PR_Only_Normal++;
                } else if (p0.Location == Overlap.INTERSECT || p1.Location == Overlap.INTERSECT) {
                    stats0.PR_SR_Normal++;
                }
            }
        }

        for (String s : outputMap.values())
            System.err.println(s);

        return result;
    }

    private static VariantQueryResult performQuery(final SamReader reader, final QueryInterval[] intervals,
            final Location bp1, final Location bp2) {

        final VariantQueryResult result = new VariantQueryResult();

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
                info.Category = ReadCategory.DIFF_CHROMOSOME;
                info.Location = Overlap.FILTERED;
                continue;
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
                if (read.getSupplementaryAlignmentFlag()) {
                    info.Category = ReadCategory.CHIMERIC;
                } else if (read.getNotPrimaryAlignmentFlag()) {
                    info.Category = ReadCategory.SECONDARY;
                } else if (read.getMateUnmappedFlag()) {
                    info.Category = ReadCategory.MATE_UNMAPPED;
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

        // query the position
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

        QueryInterval[] queryIntervals = {
                new QueryInterval(location1.ReferenceIndex, Math.max(0, location1.Position - range),
                        location1.Position + range),
                new QueryInterval(location2.ReferenceIndex, Math.max(0, location2.Position - range),
                        location2.Position + range) };
        queryIntervals = QueryInterval.optimizeIntervals(queryIntervals);

        // print  header
        // @formatter:off
        System.err.println(
                String.join("\t",
                        "SAMPLE",
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

        final VariantQueryResult tumorResult = performQuery(tumorReader, queryIntervals, location1, location2);
        final SampleStats tumorStats = processResult(tumorResult);
        tumorStats.Print();

        // final VariantQueryResult refResult = performQuery(refReader, queryIntervals, location1, location2);
        // processResult(refResult);
    }
}
