package com.hartwig.hmftools.breakpointinspector;

import static com.hartwig.hmftools.breakpointinspector.ReadHelpers.getClips;
import static com.hartwig.hmftools.breakpointinspector.Util.prefixList;
import static com.hartwig.hmftools.breakpointinspector.Util.toStrings;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Ordering;
import com.google.common.collect.TreeMultimap;

import htsjdk.samtools.SAMRecord;

class Stats {

    static class BreakPoint {
        int PR_Only_Normal = 0;
        int PR_SR_Normal = 0;
        int PR_Only_Support = 0;
        int PR_SR_Support = 0;
        int SR_Only_Support = 0;
        int Unmapped_Mate = 0;
        int Diff_Variant = 0;

        static List<String> GetHeader() {
            return Arrays.asList("PR_ONLY_NORMAL", "PR_SR_NORMAL", "PR_ONLY_SUPPORT", "PR_SR_SUPPORT",
                    "SR_ONLY_SUPPORT", "UNMAPPED_MATE", "DIFF_VARIANT");
        }

        List<Integer> GetData() {
            return Arrays.asList(PR_Only_Normal, PR_SR_Normal, PR_Only_Support, PR_SR_Support, SR_Only_Support,
                    Unmapped_Mate, Diff_Variant);
        }
    }

    static class Clip {
        String LongestClipSequence = "";
        List<SAMRecord> Reads = new ArrayList<>();
        List<SAMRecord> HardClippedReads = new ArrayList<>();

        @Override
        public String toString() {
            return LongestClipSequence + "," + Reads.size();
        }
    }

    static class ClipStats {
        Map<Util.Location, Clip> LocationMap = new HashMap<>();

        @Override
        public String toString() {

            final TreeMultimap<Integer, String> sortedClips = TreeMultimap.create(Collections.reverseOrder(),
                    Ordering.natural());
            for (final Map.Entry<Util.Location, Clip> kv : LocationMap.entrySet()) {
                final Util.Location alignment = kv.getKey();
                final Clip stats = kv.getValue();
                if (stats.Reads.size() < 2)
                    continue; // skip if we only have hard clips
                sortedClips.put(stats.Reads.size(), alignment + "," + stats);
            }

            return String.join(";", sortedClips.values());
        }

        void addToClippingStats(final SAMRecord read) {
            for (final Util.ClipInfo clip : getClips(read)) {
                final Stats.Clip stats = LocationMap.computeIfAbsent(clip.Alignment, k -> new Stats.Clip());
                if (clip.HardClipped) {
                    stats.HardClippedReads.add(read);
                } else {
                    if (clip.Sequence.length() > stats.LongestClipSequence.length() && (
                            clip.Sequence.startsWith(stats.LongestClipSequence) || clip.Sequence.endsWith(
                                    stats.LongestClipSequence))) {
                        // the existing sequence supports the new sequence
                        stats.LongestClipSequence = clip.Sequence;
                    } else if (!(stats.LongestClipSequence.startsWith(clip.Sequence)
                            || stats.LongestClipSequence.endsWith(clip.Sequence))) {
                        // this read does not support the existing sequence
                        continue;
                    }
                    stats.Reads.add(read);
                }
            }
        }
    }

    static class Sample {
        Util.Location BP1;
        BreakPoint BP1_Stats = new BreakPoint();
        Util.Location BP2;
        BreakPoint BP2_Stats = new BreakPoint();
        ClipStats Clipping_Stats = new ClipStats();

        BreakPoint Get(final Util.Region r) {
            switch (r) {
                case BP1:
                    return BP1_Stats;
                case BP2:
                    return BP2_Stats;
                default:
                    throw new RuntimeException("invalid stats");
            }
        }

        static List<String> GetHeader() {
            final ArrayList<String> header = new ArrayList<>();
            header.addAll(prefixList(BreakPoint.GetHeader(), "BP1_"));
            header.addAll(prefixList(BreakPoint.GetHeader(), "BP2_"));
            return header;
        }

        List<String> GetData() {
            final ArrayList<String> data = new ArrayList<>();
            data.addAll(toStrings(BP1_Stats.GetData()));
            data.addAll(toStrings(BP2_Stats.GetData()));
            return data;
        }
    }
}
