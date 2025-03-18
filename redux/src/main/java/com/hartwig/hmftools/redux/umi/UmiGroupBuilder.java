package com.hartwig.hmftools.redux.umi;

import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.concurrent.atomic.AtomicReference;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.redux.common.DuplicateGroup;
import com.hartwig.hmftools.redux.common.ReadInfo;

import htsjdk.samtools.SAMRecord;

public class UmiGroupBuilder
{
    private final UmiGroupBuilder_0 mBuilder0;
    private final UmiGroupBuilder_1 mBuilder1;

    public static final AtomicReference<Set<String>> mSmallestBad = new AtomicReference<>(null);

    public static final Set<String> DEBUG_READS = Sets.newHashSet(
            "A01524:289:HFJLTDRX3:2:2266:19153:33880:TTGGCTT_TGTCGTT\t99\t22\t16345850\t46\t143M\t=\t16345892\t185\tAAAAGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGT\tFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFF:FFFFFFF:FFFFF,,FFFFFFFFFFF\tXA:Z:14,+20077969,143M,1;14,-19495348,143M,1;17,+29541405,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:TTGGCTT_TGTCGTT\tNM:i:0\tAS:i:143\tXS:i:138",
            "A01524:289:HFJLTDRX3:1:2107:6470:28197:ATCTCCT_TGTCGTC\t99\t22\t16345854\t46\t143M\t=\t16345893\t182\tGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFF:FFFFFF:FFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-19495344,143M,1;14,+20077973,143M,1;17,+29541409,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:ATCTCCT_TGTCGTC\tNM:i:0\tAS:i:143\tXS:i:138",
            "A01524:289:HFJLTDRX3:2:2223:16957:6621:GTGTCAT_AGCATCT\t163\t22\t16345854\t46\t143M\t=\t16345892\t181\tGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFF:FFFFFFFFFFFFFF:FF:FFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\tXA:Z:14,-19495344,143M,1;14,+20077973,143M,1;17,+29541409,143M,4;\tMC:Z:143M\tMD:Z:143\tNM:i:0\tAS:i:143\tXS:i:138",
            "A01524:289:HFJLTDRX3:2:2264:27127:21699:TGTCGTC_ATCTCCT\t163\t22\t16345854\t46\t143M\t=\t16345893\t182\tGGTTGGATAGCTATTATCCTGAGTCTTATGTCTGATACCATGTTTTTGTTTTGTTTTTAGAGTCTTATATTTAAAGAAAAAGTAACAAGCCTTAAATTTAAAGAAAAACCTACAGGCTTGGAGACAAGAAGCTGTAAGTATCT\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF:,FFFFFFFFF:F::FFFFFFFF:FFFFF:FFFFFFFFFFFF\tXA:Z:14,+20077973,143M,1;14,-19495344,143M,1;17,+29541409,143M,4;\tMC:Z:143M\tMD:Z:143\tUI:Z:ATCTCCT_TGTCGTC\tNM:i:0\tAS:i:143\tXS:i:138"
            );

    public UmiGroupBuilder(final SequencingType sequencing, final UmiConfig config, final UmiStatistics stats)
    {
        mBuilder0 = new UmiGroupBuilder_0(sequencing, config, stats);
        mBuilder1 = new UmiGroupBuilder_1(sequencing, config);
    }

    public List<DuplicateGroup> processUmiGroups(
            final List<DuplicateGroup> duplicateGroups, final List<ReadInfo> singleFragments, boolean captureStats)
    {
        final Set<String> inputReads = Stream.concat(duplicateGroups.stream().flatMap(x -> x.reads().stream()), singleFragments.stream().map(x -> x.read()))
                .map(x -> x.getSAMString())
                .collect(Collectors.toCollection(Sets::newHashSet));

        List<DuplicateGroup> duplicateGroups0 = duplicateGroups.stream().map(DuplicateGroup::deepCopy).collect(Collectors.toList());
        List<DuplicateGroup> duplicateGroups1 = duplicateGroups.stream().map(DuplicateGroup::deepCopy).collect(Collectors.toList());

        List<ReadInfo> singleFragments0 = singleFragments.stream().map(ReadInfo::deepCopy).collect(Collectors.toList());
        List<ReadInfo> singleFragments1 = singleFragments.stream().map(ReadInfo::deepCopy).collect(Collectors.toList());

        List<DuplicateGroup> umiGroups0 = mBuilder0.processUmiGroups(duplicateGroups0, singleFragments0, captureStats);
        List<DuplicateGroup> umiGroups1 = mBuilder1.processUmiGroups(duplicateGroups1, singleFragments1, captureStats);

        Map<Set<String>, Integer> diff = Maps.newHashMap();
        for(ReadInfo singleRead : singleFragments1)
        {
            Set<String> key = Sets.newHashSet(singleRead.read().getSAMString());
            diff.merge(key, 1, Integer::sum);
        }

        for(DuplicateGroup umiGroup : umiGroups1)
        {
            Set<String> key = umiGroup.reads().stream().map(SAMRecord::getSAMString).collect(Collectors.toCollection(Sets::newHashSet));
            diff.merge(key, 1, Integer::sum);
        }

        for(ReadInfo singleRead : singleFragments0)
        {
            Set<String> key = Sets.newHashSet(singleRead.read().getSAMString());
            diff.merge(key, -1, Integer::sum);
        }

        for(DuplicateGroup umiGroup : umiGroups0)
        {
            Set<String> key = umiGroup.reads().stream().map(SAMRecord::getSAMString).collect(Collectors.toCollection(Sets::newHashSet));
            diff.merge(key, -1, Integer::sum);
        }

        Map<Set<String>, Integer> tmpDiff = Maps.newHashMap();
        diff.entrySet().stream().filter(x -> x.getValue() != 0).forEach(x -> tmpDiff.put(x.getKey(), x.getValue()));
        diff = tmpDiff;

        if(!diff.isEmpty() && Objects.equals(System.getProperty("env"), "test"))
        {
            throw new RuntimeException("UMI Group Builder fails");
        }

        if(!diff.isEmpty())
        {
            // Skip supp reads
            final long suppReadCount = inputReads.stream()
                    .map(x -> Integer.parseInt(x.split("\t")[1]))
                    .filter(flag -> ((flag & 2048) != 0))
                    .count();

            // TODO: test supps
            if(suppReadCount == 0)
            {
                mSmallestBad.updateAndGet(currentRef -> {
                    if(currentRef == null)
                        return inputReads;

                    if(currentRef.size() <= inputReads.size())
                        return currentRef;

                    return inputReads;
                });
            }
        }

        singleFragments.clear();
        singleFragments.addAll(singleFragments0);
        return umiGroups0;
    }
}
