package com.hartwig.hmftools.breakpointinspector;

import static com.hartwig.hmftools.breakpointinspector.Util.prefixList;
import static com.hartwig.hmftools.breakpointinspector.Util.toStrings;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;

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

    static class Orientation {
        int InnieCount = 0;
        int OutieCount = 0;
        int TandemCount = 0;

        static List<String> GetHeader() {
            return Arrays.asList("INNIE", "OUTIE", "TANDEM");
        }

        List<Integer> GetData() {
            return Arrays.asList(InnieCount, OutieCount, TandemCount);
        }
    }

    static class Clip {
        Util.ClipSide Side = Util.ClipSide.NONE;
        String LongestClipSequence = "";
        List<SAMRecord> Reads = new ArrayList<>();
        List<SAMRecord> HardClippedReads = new ArrayList<>();

        @Override
        public String toString() {
            return (Side == Util.ClipSide.RIGHT_CLIP ? "*" : "") + LongestClipSequence + (
                    Side == Util.ClipSide.LEFT_CLIP ? "*" : "") + "," + Reads.size();
        }
    }

    static class ClipStats {
        Map<Util.Location, Clip> LocationMap = new Hashtable<>();

        @Override
        public String toString() {
            final TreeMultimap<Util.Location, String> sortedClips = TreeMultimap.create();
            for (final Map.Entry<Util.Location, Clip> kv : LocationMap.entrySet()) {
                final Util.Location alignment = kv.getKey();
                final Clip stats = kv.getValue();
                sortedClips.put(alignment, alignment + "," + stats);
            }

            return String.join(";", sortedClips.values());
        }
    }

    static class Sample {
        BreakPoint BP1_Stats = new BreakPoint();
        BreakPoint BP2_Stats = new BreakPoint();
        Orientation Orientation_Stats = new Orientation();
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
            header.addAll(Orientation.GetHeader());
            return header;
        }

        List<String> GetData() {
            final ArrayList<String> data = new ArrayList<>();
            data.addAll(toStrings(BP1_Stats.GetData()));
            data.addAll(toStrings(BP2_Stats.GetData()));
            data.addAll(toStrings(Orientation_Stats.GetData()));
            return data;
        }
    }
}
