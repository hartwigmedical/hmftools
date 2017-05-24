package com.hartwig.hmftools.breakpointinspector;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import htsjdk.samtools.*;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import static com.hartwig.hmftools.breakpointinspector.Util.*;
import static com.hartwig.hmftools.breakpointinspector.Stats.*;

import org.jetbrains.annotations.NotNull;

public class BreakPointInspectorApplication {

    private static final String REF_PATH = "ref";
    private static final String TUMOR_PATH = "tumor";
    private static final String BREAK_POINT1 = "bp1";
    private static final String BREAK_POINT2 = "bp2";
    private static final String PROXIMITY = "proximity";
    private static final String SV_LEN = "svlen";
    private static final String VCF = "vcf";

    private static final int MANTA_REQ_PAIR_MIN = 50;
    private static final int MANTA_REQ_SPLIT_MIN = 15;

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(REF_PATH, true, "the Reference BAM");
        options.addOption(TUMOR_PATH, true, "the Tumor BAM");
        options.addOption(BREAK_POINT1, true, "position of first break point in chrX:123456 format");
        options.addOption(BREAK_POINT2, true, "position of second break point in chrX:123456 format (optional)");
        options.addOption(PROXIMITY, true, "base distance around breakpoint");
        options.addOption(SV_LEN, true, "length of the SV to inspect");
        options.addOption(VCF, true, "VCF file to batch inspect");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
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
            case H:
                result.HardClipped = true;
            case S:
                result.Side = ClipSide.LEFT_CLIP;
                result.Length = cigar.getFirstCigarElement().getLength();
                return result;
        }
        switch (cigar.getLastCigarElement().getOperator()) {
            case H:
                result.HardClipped = true;
            case S:
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

    private static boolean isMate(final SAMRecord read, final SAMRecord mate) {
        return read.getReadName().equals(mate.getReadName())
                && read.getMateReferenceIndex() == mate.getReferenceIndex()
                && read.getMateAlignmentStart() == mate.getAlignmentStart();
    }

    private static Stats.ClipStats calculateClippingStats(final ClassifiedReadResults queryResult) {
        final Stats.ClipStats result = new Stats.ClipStats();
        for (final NamedReadCollection collection : queryResult.ReadMap.values()) {
            for (final ReadInfo info : collection.Reads) {
                final SAMRecord read = info.Read;
                if (info.Location == Overlap.CLIP) {
                    // TODO: handle multiple clip sides
                    final Location alignment = Location.fromSAMRecord(read, info.Clipping.Side == ClipSide.LEFT_CLIP);
                    final Stats.Clip clip = result.LocationMap.computeIfAbsent(alignment, k -> new Stats.Clip() {{
                        Side = info.Clipping.Side;
                    }});
                    if (info.Clipping.HardClipped) {
                        clip.HardClippedReads.add(read);
                    } else {
                        final String clippedSequence;
                        if (info.Clipping.Side == ClipSide.LEFT_CLIP) {
                            clippedSequence = read.getReadString().substring(0, info.Clipping.Length);
                        } else {
                            clippedSequence = read.getReadString().substring(
                                    read.getReadLength() - info.Clipping.Length);
                        }

                        if (clippedSequence.length() > clip.LongestClipSequence.length() && clippedSequence.contains(
                                clip.LongestClipSequence)) {
                            // the existing sequence supports the new sequence
                            clip.LongestClipSequence = clippedSequence;
                        } else if (!clip.LongestClipSequence.contains(clippedSequence)) {
                            // this read does not support the existing sequence
                            continue;
                        }

                        clip.Reads.add(read);
                    }
                }
            }
        }
        return result;
    }

    private static Stats.Sample calculateEvidenceStats(final ClassifiedReadResults queryResult) {
        final Stats.Sample result = new Stats.Sample();
        for (final NamedReadCollection collection : queryResult.ReadMap.values()) {
            // consider the pairings
            for (final ReadInfo p0 : collection.Reads) {
                // find the mate
                final ReadInfo p1 = collection.Reads.stream().filter(i -> isMate(p0.Read, i.Read)).findFirst().orElse(
                        null);

                // single read
                if (p1 == null) {
                    if (p0.Category == ReadCategory.MATE_UNMAPPED) {
                        result.Get(p0.Breakpoint).Unmapped_Mate++;
                    } else if (p0.Category != ReadCategory.NORMAL || p0.Location != Overlap.FILTERED) {
                        result.Get(p0.Breakpoint).Diff_Variant++;
                    }
                    continue;
                }

                // don't consider pairs twice from the reverse pairing
                if (p0.Read.getInferredInsertSize() < 0) {
                    continue;
                }

                // possible two paired but unmapped reads?
                if (p0.Breakpoint == Region.OTHER && p1.Breakpoint == Region.OTHER) {
                    continue;
                } else if (p1.Breakpoint == Region.OTHER) { // must be unmapped mate
                    result.Get(p0.Breakpoint).Unmapped_Mate++;
                    continue;
                } else if (p0.Breakpoint == Region.OTHER) { // must be unmapped mate
                    result.Get(p1.Breakpoint).Unmapped_Mate++;
                    continue;
                }

                final List<ReadInfo> pair = Arrays.asList(p0, p1);
                final boolean pairPossible = pair.stream().allMatch(p -> p.Clipping.Length <= MANTA_REQ_PAIR_MIN);

                if (p0.Breakpoint != p1.Breakpoint) {
                    // supports the break point
                    boolean interesting = false;
                    for (final ReadInfo r : pair) {
                        final boolean splitEvidence =
                                r.Location == Overlap.CLIP && r.Clipping.Length >= MANTA_REQ_SPLIT_MIN;
                        if (pairPossible && splitEvidence) {
                            result.Get(r.Breakpoint).PR_SR_Support++;
                        } else if (pairPossible) {
                            result.Get(r.Breakpoint).PR_Only_Support++;
                        } else if (splitEvidence) {
                            result.Get(r.Breakpoint).SR_Only_Support++;
                        }

                        interesting |= pairPossible || splitEvidence;
                    }

                    if (interesting) {
                        switch (SamPairUtil.getPairOrientation(p0.Read)) {
                            case FR:
                                result.Orientation_Stats.InnieCount++;
                                break;
                            case RF:
                                result.Orientation_Stats.OutieCount++;
                                break;
                            case TANDEM:
                                result.Orientation_Stats.TandemCount++;
                                break;
                        }
                    }
                } else {
                    final Stats.BreakPoint stats = p0.Breakpoint == Region.BP1 ? result.BP1_Stats : result.BP2_Stats;
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
        }
        return result;
    }

    private static Stats.Sample calculateStats(final ClassifiedReadResults queryResult) {
        final Stats.Sample result = calculateEvidenceStats(queryResult);
        result.Clipping_Stats = calculateClippingStats(queryResult);
        return result;
    }

    private static ClassifiedReadResults performQueryAndClassify(final SamReader reader,
            final QueryInterval[] intervals, final Location bp1, final Location bp2) {

        final ClassifiedReadResults result = new ClassifiedReadResults();

        // execute query and parse the results
        final SAMRecordIterator results = reader.query(intervals, false);
        while (results.hasNext()) {

            final SAMRecord read = results.next();
            final NamedReadCollection collection = result.ReadMap.computeIfAbsent(read.getReadName(),
                    k -> new NamedReadCollection());
            final ReadInfo info = new ReadInfo();

            info.Read = read;
            collection.Reads.add(info);

            // if unmapped there's nothing to do
            if (read.getReadUnmappedFlag()) {
                info.Breakpoint = Region.OTHER;
                info.Category = ReadCategory.UNMAPPED;
                info.Location = Overlap.FILTERED;
                continue;
            }

            info.Clipping = getClipInfo(read);
            info.Breakpoint = determineRegion(read, bp1, bp2);
            if (info.Breakpoint == Region.OTHER) {
                throw new RuntimeException("read from unexpected region");
            }

            final Location bp = info.Breakpoint == Region.BP1 ? bp1 : bp2;
            final boolean normal = read.getProperPairFlag();

            if (normal) {
                info.Category = ReadCategory.NORMAL;

                if (info.Clipping.Side == ClipSide.NONE) {
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
            // TODO: double check clipping conditions
            if (info.Clipping.Side == ClipSide.LEFT_CLIP && read.getAlignmentStart() - 1 == bp.Position) {
                info.Location = Overlap.CLIP;
            } else if (info.Clipping.Side == ClipSide.RIGHT_CLIP && read.getAlignmentEnd() == bp.Position) {
                info.Location = Overlap.CLIP;
            } else {
                info.Location = Overlap.PROXIMITY;
            }

            // if normal read and we don't clip, then it's not interesting
            if (normal && info.Location == Overlap.PROXIMITY)
                info.Location = Overlap.FILTERED;
        }
        results.close();

        return result;
    }

    private static void processStructuralVariant(final List<String> headers, final SamReader refReader,
            final SamReader tumorReader, final Location location1, final Location location2, final int range) {
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
        final Sample refStats = calculateStats(refResult);

        final ClassifiedReadResults tumorResult = performQueryAndClassify(tumorReader, queryIntervals, location1,
                location2);
        final Sample tumorStats = calculateStats(tumorResult);

        final ArrayList<String> data = new ArrayList<>(headers);
        data.addAll(refStats.GetData());
        data.addAll(tumorStats.GetData());
        data.add(tumorStats.Clipping_Stats.toString());
        System.out.println(String.join("\t", data));
    }

    public static void main(final String... args) throws ParseException, IOException {

        final Options options = createOptions();
        try {
            final CommandLine cmd = createCommandLine(options, args);

            // grab arguments
            final String refBAM = cmd.getOptionValue(REF_PATH);
            final String tumorBAM = cmd.getOptionValue(TUMOR_PATH);
            final String bp1String = cmd.getOptionValue(BREAK_POINT1);
            final String bp2String = cmd.getOptionValue(BREAK_POINT2);
            final String vcfPath = cmd.getOptionValue(VCF);
            final int range = Integer.parseInt(cmd.getOptionValue(PROXIMITY, "500"));

            if (refBAM == null || tumorBAM == null || (vcfPath != null && bp1String != null) || (bp2String != null
                    && cmd.getOptionValue(SV_LEN) != null))
                printHelpAndExit(options);

            // load the files
            final File tumorFile = new File(tumorBAM);
            final SamReader tumorReader = SamReaderFactory.makeDefault().open(tumorFile);
            final File refFile = new File(refBAM);
            final SamReader refReader = SamReaderFactory.makeDefault().open(refFile);

            // output the header
            final ArrayList<String> header = new ArrayList<>(
                    Arrays.asList("ID", "MANTA_BP1", "MANTA_BP2", "MANTA_SVLEN", "MANTA_HOMSEQ", "MANTA_INSSEQ"));
            header.addAll(prefixList(Sample.GetHeader(), "REF_"));
            header.addAll(prefixList(Sample.GetHeader(), "TUMOR_"));
            header.add("TUMOR_CLIP_INFO");
            System.out.println(String.join("\t", header));

            if (vcfPath != null) {
                final File vcfFile = new File(vcfPath);
                final VCFFileReader vcfReader = new VCFFileReader(vcfFile, false);

                for (final VariantContext variant : vcfReader) {
                    if (variant.isFiltered())
                        continue;

                    final String location = variant.getContig() + ":" + Integer.toString(variant.getStart());
                    final Location location1 = Location.parseLocationString(location,
                            tumorReader.getFileHeader().getSequenceDictionary());
                    final Location location2;

                    final String variantType = variant.getAttributeAsString("SVTYPE", "");
                    switch (variantType) {
                        case "DEL":
                        case "DUP":
                        case "INV":
                            final int svLen = Math.abs(variant.getAttributeAsInt("SVLEN", 0));
                            location2 = location1.add(svLen);
                            break;
                        case "BND":
                            final String call = variant.getAlternateAllele(0).getDisplayString();
                            final String[] split = call.split("[\\]\\[]");
                            location2 = Location.parseLocationString(split[1],
                                    tumorReader.getFileHeader().getSequenceDictionary());
                            break;
                        default:
                            System.err.println(variant.getID() + " : UNEXPECTED SVTYPE=" + variantType);
                            continue;
                    }

                    final List<String> extraHeaders = Arrays.asList(variant.getID(), location1.toString(),
                            location2.toString(), variant.getAttributeAsString("SVLEN", "."),
                            variant.getAttributeAsString("HOMSEQ", "."),
                            variant.getAttributeAsString("SVINSSEQ", "."));
                    processStructuralVariant(extraHeaders, refReader, tumorReader, location1, location2, range);
                }
            } else {
                final int svLen = Integer.parseInt(cmd.getOptionValue(SV_LEN, "0"));

                // parse the location strings
                final Location location1 = Location.parseLocationString(bp1String,
                        tumorReader.getFileHeader().getSequenceDictionary());
                final Location location2;
                if (bp2String != null) {
                    location2 = Location.parseLocationString(bp2String,
                            tumorReader.getFileHeader().getSequenceDictionary());
                } else if (svLen > 0) {
                    location2 = location1.add(svLen);
                } else {
                    printHelpAndExit(options);
                    System.exit(1);
                    return;
                }

                final List<String> extraHeaders = Arrays.asList("manual", location1.toString(), location2.toString(),
                        svLen > 0 ? Integer.toString(svLen) : ".", ".", ".");
                processStructuralVariant(extraHeaders, refReader, tumorReader, location1, location2, range);
            }

        } catch (ParseException e) {
            printHelpAndExit(options);
            System.exit(1);
        }
    }
}
