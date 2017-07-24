package com.hartwig.hmftools.breakpointinspector;

import static com.hartwig.hmftools.breakpointinspector.Stats.SampleStats;
import static com.hartwig.hmftools.breakpointinspector.Util.HMFVariantContext;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Stream;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

class Analysis {

    private static class AlignmentList extends ArrayList<SAMRecord> {
    }

    private static class ReadMap extends HashMap<String, AlignmentList> {
    }

    private static class PairedReads extends ArrayList<Pair<SAMRecord, SAMRecord>> {
    }

    static class StructuralVariantResult {
        SampleStats TumorStats = new SampleStats();
        SampleStats RefStats = new SampleStats();
        String Filter = "";
    }

    private static Pair<Integer, Integer> orientation(final Pair<SAMRecord, SAMRecord> pair) {
        return Pair.of(pair.getLeft().getReadNegativeStrandFlag() ? -1 : 1,
                !pair.getRight().getReadNegativeStrandFlag() ? -1 : 1);
    }

    private static void determineBreakpoints(final HMFVariantContext ctx, final PairedReads pairs) {

        final Pair<Integer, Integer> ctxOrientation = Pair.of(ctx.OrientationBP1, ctx.OrientationBP2);

        // extract all interesting pairs

        final PairedReads interesting = new PairedReads();

        for (final Pair<SAMRecord, SAMRecord> pair : pairs) {

            final boolean correctOrientation = orientation(pair).equals(ctxOrientation);
            final boolean proper = Stream.of(pair.getLeft(), pair.getRight()).anyMatch(SAMRecord::getProperPairFlag);
            final boolean clipped = Stream.of(pair.getLeft(), pair.getRight()).anyMatch(r -> r.getCigar().isClipped());

            // TODO: more precise definition of interesting?
            if ((correctOrientation && !proper) || (proper && clipped)) {
                interesting.add(pair);
            }
        }

        // TODO: filter for sequences within Manta uncertainty

        // examine clipping

        Clipping bp1_clipping = new Clipping();
        Clipping bp2_clipping = new Clipping();
        for (final Pair<SAMRecord, SAMRecord> pair : interesting) {
            if (ctx.OrientationBP1 > 0) {
                bp1_clipping.add(Clipping.getRightClip(pair.getLeft()));
                bp1_clipping.add(Clipping.getRightClip(pair.getRight()));
            } else {
                bp1_clipping.add(Clipping.getLeftClip(pair.getLeft()));
                bp1_clipping.add(Clipping.getLeftClip(pair.getLeft()));
            }

            if (ctx.OrientationBP2 > 0) {
                bp2_clipping.add(Clipping.getLeftClip(pair.getRight()));
                bp2_clipping.add(Clipping.getLeftClip(pair.getLeft()));
            } else {
                bp2_clipping.add(Clipping.getRightClip(pair.getRight()));
                bp2_clipping.add(Clipping.getRightClip(pair.getLeft()));
            }
        }

        final ClipStats bp1_candidate = bp1_clipping.getTopSequence();
        final ClipStats bp2_candidate = bp2_clipping.getTopSequence();

        if (bp1_candidate == null) {
            // TODO: fail-over logic
        } else {
            ctx.BP1 = ctx.OrientationBP1 > 0 ? bp1_candidate.Alignment.add(1) : bp1_candidate.Alignment.add(-1);
        }

        if (bp2_candidate == null) {
            // TODO: fail-over logic
        } else {
            ctx.BP2 = ctx.OrientationBP2 > 0 ? bp2_candidate.Alignment.add(1) : bp2_candidate.Alignment.add(-1);
        }
    }

    private static PairedReads pairs(final ReadMap reads) {

        final PairedReads pairs = new PairedReads();

        for (final AlignmentList alignments : reads.values()) {
            for (int i = 0; i < alignments.size(); ++i) {
                final SAMRecord r0 = alignments.get(i);
                for (int j = i + 1; j < alignments.size(); ++j) {
                    final SAMRecord r1 = alignments.get(j);
                    if (ReadHelpers.isMate(r0, r1)) {
                        pairs.add(Pair.of(r0, r1));
                    }
                }
            }
        }

        return pairs;
    }

    private static ReadMap readsByName(final List<SAMRecord> alignments) {
        final ReadMap result = new ReadMap();
        alignments.forEach(a -> result.computeIfAbsent(a.getReadName(), k -> new AlignmentList()).add(a));
        return result;
    }

    private static List<SAMRecord> performQuery(final SamReader reader, final QueryInterval[] intervals) {
        return reader.queryOverlapping(intervals).toList();
    }

    private final static Set<Integer> tumorWrittenReads = new HashSet<>();
    private final static Set<Integer> refWrittenReads = new HashSet<>();
    static StructuralVariantResult processStructuralVariant(final SamReader refReader, @Nullable final SAMFileWriter refWriter,
            final SamReader tumorReader, @Nullable final SAMFileWriter tumorWriter, final HMFVariantContext ctx, final int range) {

        // perform query for reads

        QueryInterval[] queryIntervals =
                { new QueryInterval(ctx.MantaBP1.ReferenceIndex, Math.max(0, ctx.MantaBP1.Position - range), ctx.MantaBP1.Position + range),
                        new QueryInterval(ctx.MantaBP2.ReferenceIndex, Math.max(0, ctx.MantaBP2.Position - range),
                                ctx.MantaBP2.Position + range) };
        queryIntervals = QueryInterval.optimizeIntervals(queryIntervals);

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

        final ReadMap tumorReadMap = readsByName(tumorReads);
        final PairedReads tumorPairedReads = pairs(tumorReadMap);
        determineBreakpoints(ctx, tumorPairedReads);

        final StructuralVariantResult result = new StructuralVariantResult();
        result.Filter = Filter.getFilterString(ctx, result.TumorStats, result.RefStats);
        return result;
    }
}
