package com.hartwig.hmftools.breakpointinspector;

import static com.hartwig.hmftools.breakpointinspector.ReadHelpers.determineRegion;
import static com.hartwig.hmftools.breakpointinspector.ReadHelpers.getClips;
import static com.hartwig.hmftools.breakpointinspector.ReadHelpers.isClipped;
import static com.hartwig.hmftools.breakpointinspector.ReadHelpers.isMate;
import static com.hartwig.hmftools.breakpointinspector.ReadHelpers.pairStraddlesLocation;
import static com.hartwig.hmftools.breakpointinspector.ReadHelpers.readIntersectsLocation;
import static com.hartwig.hmftools.breakpointinspector.Stats.BreakpointStats;
import static com.hartwig.hmftools.breakpointinspector.Stats.ClipStats;
import static com.hartwig.hmftools.breakpointinspector.Stats.SampleStats;
import static com.hartwig.hmftools.breakpointinspector.Util.ClassifiedReadResults;
import static com.hartwig.hmftools.breakpointinspector.Util.HMFVariantContext;
import static com.hartwig.hmftools.breakpointinspector.Util.Location;
import static com.hartwig.hmftools.breakpointinspector.Util.NamedReadCollection;
import static com.hartwig.hmftools.breakpointinspector.Util.Overlap;
import static com.hartwig.hmftools.breakpointinspector.Util.ReadCategory;
import static com.hartwig.hmftools.breakpointinspector.Util.ReadInfo;
import static com.hartwig.hmftools.breakpointinspector.Util.Region;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

import com.google.common.collect.Lists;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

class Analysis {

    private static Set<Integer> tumorWrittenReads = new HashSet<>();
    private static Set<Integer> refWrittenReads = new HashSet<>();

    private static boolean assessDEL(final HMFVariantContext ctx, Pair<ReadInfo, ReadInfo> pair) {

        final Location adjustedBP1 = ctx.Uncertainty1 == null ? ctx.MantaBP1 : ctx.MantaBP1.add(ctx.Uncertainty1.End);
        final Location adjustedBP2 = ctx.Uncertainty2 == null ? ctx.MantaBP2 : ctx.MantaBP2.add(ctx.Uncertainty2.Start);

        boolean pairEvidence;
        // side
        pairEvidence = pair.getLeft().Read.getAlignmentStart() <= adjustedBP1.Position;
        pairEvidence &= adjustedBP2.Position <= pair.getRight().Read.getAlignmentEnd();
        // orientation
        pairEvidence &= !pair.getLeft().Read.getReadNegativeStrandFlag();
        pairEvidence &= pair.getRight().Read.getReadNegativeStrandFlag();

        return pairEvidence;
    }

    private static boolean assessDUP(final HMFVariantContext ctx, final Pair<ReadInfo, ReadInfo> pair) {

        final Location adjustedBP1 = ctx.Uncertainty1 == null ? ctx.MantaBP1 : ctx.MantaBP1.add(ctx.Uncertainty1.Start);
        final Location adjustedBP2 = ctx.Uncertainty2 == null ? ctx.MantaBP2 : ctx.MantaBP2.add(ctx.Uncertainty2.End);

        boolean pairEvidence;
        // side
        pairEvidence = adjustedBP1.Position <= pair.getLeft().Read.getAlignmentEnd();
        pairEvidence &= pair.getRight().Read.getAlignmentStart() <= adjustedBP2.Position;
        // orientation
        pairEvidence &= pair.getLeft().Read.getReadNegativeStrandFlag();
        pairEvidence &= !pair.getRight().Read.getReadNegativeStrandFlag();

        return pairEvidence;
    }

    private static boolean assessINV(final HMFVariantContext ctx, final Pair<ReadInfo, ReadInfo> pair, boolean inv3) {

        boolean pairEvidence;
        if (inv3) {
            final Location adjustedBP1 = ctx.Uncertainty1 == null ? ctx.MantaBP1 : ctx.MantaBP1.add(ctx.Uncertainty1.End);
            final Location adjustedBP2 = ctx.Uncertainty2 == null ? ctx.MantaBP2 : ctx.MantaBP2.add(ctx.Uncertainty2.End);

            pairEvidence = pair.getLeft().Read.getAlignmentStart() <= adjustedBP1.Position;
            pairEvidence &= pair.getRight().Read.getAlignmentStart() <= adjustedBP2.Position;
            pairEvidence &= !pair.getLeft().Read.getReadNegativeStrandFlag(); // forward strand x ---->
            pairEvidence &= !pair.getRight().Read.getReadNegativeStrandFlag(); // forward strand x ---->
        } else {
            final Location adjustedBP1 = ctx.Uncertainty1 == null ? ctx.MantaBP1 : ctx.MantaBP1.add(ctx.Uncertainty1.Start);
            final Location adjustedBP2 = ctx.Uncertainty2 == null ? ctx.MantaBP2 : ctx.MantaBP2.add(ctx.Uncertainty2.Start);

            pairEvidence = adjustedBP1.Position <= pair.getLeft().Read.getAlignmentEnd();
            pairEvidence &= adjustedBP2.Position <= pair.getRight().Read.getAlignmentEnd();
            pairEvidence &= pair.getLeft().Read.getReadNegativeStrandFlag(); // reverse strand <---- x
            pairEvidence &= pair.getRight().Read.getReadNegativeStrandFlag(); // reverse strand <---- x
        }

        return pairEvidence;
    }

    private static boolean assessPR(final HMFVariantContext ctx, final Pair<ReadInfo, ReadInfo> pair) {
        switch (ctx.Type) {
            case DEL:
                return assessDEL(ctx, pair);
            case DUP:
                return assessDUP(ctx, pair);
            case INV3:
                return assessINV(ctx, pair, true);
            case INV5:
                return assessINV(ctx, pair, false);
        }
        return false;
    }

    private static void determineBreakpoints(final HMFVariantContext ctx, final SampleStats result, final ClipStats bpClips1,
            final ClipStats bpClips2) {

        // TODO: needs finesse

        // 1: look at clipped pairs first
        if (result.BP1 == null) {
            result.BP1 = bpClips1.LocationMap.entrySet()
                    .stream()
                    .filter(e -> e.getKey().closeTo(ctx.MantaBP1, ctx.Uncertainty1))
                    .max(Comparator.comparingInt(a -> a.getValue().Reads.size()))
                    .map(Map.Entry::getKey)
                    .orElse(null);
        }
        if (result.BP2 == null) {
            result.BP2 = bpClips2.LocationMap.entrySet()
                    .stream()
                    .filter(e -> e.getKey().closeTo(ctx.MantaBP2, ctx.Uncertainty2))
                    .max(Comparator.comparingInt(a -> a.getValue().Reads.size()))
                    .map(Map.Entry::getKey)
                    .orElse(null);
        }

        // 2: then we should look at unpaired clips to determine breakpoint
        if (result.BP1 == null) {
            result.BP1 = result.Clipping_Stats.LocationMap.entrySet()
                    .stream()
                    .filter(e -> e.getKey().closeTo(ctx.MantaBP1, ctx.Uncertainty1))
                    .max(Comparator.comparingInt(a -> a.getValue().Reads.size()))
                    .map(Map.Entry::getKey)
                    .orElse(null);
        }
        if (result.BP2 == null) {
            result.BP2 = result.Clipping_Stats.LocationMap.entrySet()
                    .stream()
                    .filter(e -> e.getKey().closeTo(ctx.MantaBP2, ctx.Uncertainty2))
                    .max(Comparator.comparingInt(a -> a.getValue().Reads.size()))
                    .map(Map.Entry::getKey)
                    .orElse(null);
        }

        // TODO: should we display Manta BP?
    }

    private static SampleStats calculateStats(final HMFVariantContext ctx, final ClassifiedReadResults queryResult) {

        final SampleStats result = new SampleStats();
        result.BP1 = ctx.BP1;
        result.BP2 = ctx.BP2;

        final Stats.ClipStats clipPRSR1 = new Stats.ClipStats();
        final Stats.ClipStats clipPRSR2 = new Stats.ClipStats();

        final List<Pair<ReadInfo, ReadInfo>> spanningPairs = Lists.newArrayList();
        final List<Pair<ReadInfo, ReadInfo>> localPairs = Lists.newArrayList();

        for (final NamedReadCollection collection : queryResult.values()) {
            for (final ReadInfo r0 : collection) {

                // update clipping stats
                result.Clipping_Stats.addToClippingStats(r0.Read);

                // find the mate
                final ReadInfo r1 = collection.stream().filter(r -> isMate(r0.Read, r.Read)).findFirst().orElse(null);

                // single read
                if (r1 == null) {
                    if (r0.Category == ReadCategory.MATE_UNMAPPED) {
                        result.Get(r0.Breakpoint).Unmapped_Mate++;
                    } else if (r0.Category != ReadCategory.NORMAL && r0.Location != Overlap.FILTERED) { // TODO: check condition
                        result.Get(r0.Breakpoint).Diff_Variant++;
                    }
                    continue;
                }

                // special case for BND pairings, we want them in same order as manta breakpoint
                if (r0.Read.getInferredInsertSize() == 0 && !r0.Read.getReferenceIndex().equals(ctx.MantaBP1.ReferenceIndex)) {
                    continue;
                }

                // don't consider pairs twice from the reverse pairing
                if (r0.Read.getInferredInsertSize() < 0) {
                    continue;
                }

                // possible two paired but unmapped reads?
                if (r0.Breakpoint == Region.OTHER && r1.Breakpoint == Region.OTHER) {
                    continue;
                } else if (r1.Breakpoint == Region.OTHER) { // must be unmapped mate
                    result.Get(r0.Breakpoint).Unmapped_Mate++;
                    continue;
                } else if (r0.Breakpoint == Region.OTHER) { // must be unmapped mate
                    result.Get(r1.Breakpoint).Unmapped_Mate++;
                    continue;
                }

                // at this stage, we have a legitimate pair
                final Pair<ReadInfo, ReadInfo> pair = Pair.of(r0, r1);

                // supports the break point
                if (r0.Breakpoint != r1.Breakpoint) {
                    if (assessPR(ctx, pair)) {
                        spanningPairs.add(pair);
                        clipPRSR1.addToClippingStats(r0.Read);
                        clipPRSR2.addToClippingStats(r1.Read);
                    } else {
                        result.BP1_Stats.Diff_Variant++;
                        result.BP2_Stats.Diff_Variant++;
                    }
                } else {
                    // p0 == p1
                    localPairs.add(pair);
                }
            }
        }

        determineBreakpoints(ctx, result, clipPRSR1, clipPRSR2);

        // look at PR evidence
        for (final Pair<ReadInfo, ReadInfo> pair : spanningPairs) {
            for (final ReadInfo r : Arrays.asList(pair.getLeft(), pair.getRight())) {
                final boolean sr =
                        getClips(r.Read).anyMatch(c -> c.Alignment.closeTo(r.Breakpoint == Region.BP1 ? result.BP1 : result.BP2));
                // TODO: also check side of clip
                if (sr) {
                    result.Get(r.Breakpoint).PR_SR_Support++;
                } else {
                    result.Get(r.Breakpoint).PR_Only_Support++;
                }
            }
        }

        // look at normal or SR evidence
        for (final Pair<ReadInfo, ReadInfo> pair : localPairs) {
            final BreakpointStats stats = result.Get(pair.getLeft().Breakpoint);
            final boolean sr = Stream.of(pair.getLeft(), pair.getRight())
                    .anyMatch(
                            r -> getClips(r.Read).anyMatch(c -> c.Alignment.closeTo(r.Breakpoint == Region.BP1 ? result.BP1 : result.BP2)));
            // TODO: also check side of clip
            if (sr) {
                stats.SR_Only_Support++;
            } else if (pair.getLeft().Location == Overlap.STRADDLE && pair.getRight().Location == Overlap.STRADDLE) {
                stats.PR_Only_Normal++;
            } else if (pair.getLeft().Location == Overlap.INTERSECT || pair.getRight().Location == Overlap.INTERSECT) {
                stats.PR_SR_Normal++;
            }
        }

        return result;
    }

    private static ClassifiedReadResults classifyReads(final HMFVariantContext ctx, final List<SAMRecord> reads) {

        final ClassifiedReadResults result = new ClassifiedReadResults();

        for (final SAMRecord read : reads) {

            final ReadInfo info = new ReadInfo(read);
            final NamedReadCollection collection = result.computeIfAbsent(read.getReadName(), k -> new NamedReadCollection());
            collection.add(info);

            // if unmapped there's nothing to do
            if (read.getReadUnmappedFlag()) {
                info.Breakpoint = Region.OTHER;
                info.Category = ReadCategory.UNMAPPED;
                info.Location = Overlap.FILTERED;
                continue;
            }

            info.Breakpoint = determineRegion(read, ctx.MantaBP1, ctx.MantaBP2);
            if (info.Breakpoint == Region.OTHER) {
                throw new RuntimeException("read from unexpected region");
            }

            final Location bp = info.Breakpoint == Region.BP1 ? ctx.MantaBP1 : ctx.MantaBP2;
            final boolean normal = read.getProperPairFlag();
            final boolean clipped = isClipped(info.Read);

            if (normal) {
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
                } else if (read.getReadPairedFlag()) {
                    info.Category = ReadCategory.SPAN;
                }
            }

            // determine classification
            if (clipped) {
                info.Location = Overlap.CLIP;
            } else {
                info.Location = Overlap.PROXIMITY;
            }
        }

        return result;
    }

    private static List<SAMRecord> performQuery(final SamReader reader, final QueryInterval[] intervals) {
        final List<SAMRecord> output = Lists.newArrayList();

        final SAMRecordIterator results = reader.query(intervals, false);
        while (results.hasNext()) {
            final SAMRecord read = results.next();
            output.add(read);
        }
        results.close();

        return output;
    }

    static StructuralVariantResult processStructuralVariant(final SamReader refReader, @Nullable final SAMFileWriter refWriter,
            final SamReader tumorReader, @Nullable final SAMFileWriter tumorWriter, final HMFVariantContext ctx, final int range) {

        // perform query for reads

        QueryInterval[] queryIntervals =
                { new QueryInterval(ctx.MantaBP1.ReferenceIndex, Math.max(0, ctx.MantaBP1.Position - range), ctx.MantaBP1.Position + range),
                        new QueryInterval(ctx.MantaBP2.ReferenceIndex, Math.max(0, ctx.MantaBP2.Position - range),
                                ctx.MantaBP2.Position + range) };
        queryIntervals = QueryInterval.optimizeIntervals(queryIntervals);

        final StructuralVariantResult result = new StructuralVariantResult();
        final List<SAMRecord> tumorReads = performQuery(tumorReader, queryIntervals);
        final List<SAMRecord> refReads = performQuery(refReader, queryIntervals);

        // write evidence

        if (tumorWriter != null) {
            for (final SAMRecord read : tumorReads) {
                final int hash = read.getSAMString().hashCode();
                if (!tumorWrittenReads.contains(hash)) {
                    tumorWriter.addAlignment(read);
                    tumorWrittenReads.add(hash);
                }
            }
        }
        if (refWriter != null) {
            for (final SAMRecord read : refReads) {
                final int hash = read.getSAMString().hashCode();
                if (!refWrittenReads.contains(hash)) {
                    refWriter.addAlignment(read);
                    refWrittenReads.add(hash);
                }
            }
        }

        // processing

        result.TumorStats = calculateStats(ctx, classifyReads(ctx, tumorReads));

        // TODO: better way to propogate this?
        ctx.BP1 = result.TumorStats.BP1;
        ctx.BP2 = result.TumorStats.BP2;

        result.RefStats = calculateStats(ctx, classifyReads(ctx, refReads));

        // filtering

        result.Filter = Filter.getFilterString(ctx, result.TumorStats, result.RefStats);

        // output

        return result;
    }

    static class StructuralVariantResult {
        SampleStats TumorStats;
        SampleStats RefStats;
        String Filter;
    }
}
