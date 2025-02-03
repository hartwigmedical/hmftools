package com.hartwig.hmftools.esvee.assembly.alignment;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.MAX_OBSERVED_CONCORDANT_FRAG_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.LOCAL_INDEL;
import static com.hartwig.hmftools.esvee.common.CommonUtils.isIndel;
import static com.hartwig.hmftools.esvee.common.CommonUtils.isShortLocalDelDupIns;
import static com.hartwig.hmftools.esvee.common.SvConstants.DEFAULT_MAX_CONCORDANT_FRAG_LENGTH;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.PhaseSet;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;

public class AlignmentFragments
{
    private final AssemblyAlignment mAssemblyAlignment;

    public AlignmentFragments(final AssemblyAlignment assemblyAlignment, final List<String> combinedSampleIds)
    {
        mAssemblyAlignment = assemblyAlignment;

        for(Breakend breakend : assemblyAlignment.breakends())
        {
            // rather than use the genome position of a read vs the aligned breakend position, use its position in the assembly
            List<BreakendSupport> sampleSupport = breakend.sampleSupport();

            combinedSampleIds.forEach(x -> sampleSupport.add(new BreakendSupport()));
        }
    }

    private class FragmentReads
    {
        public final List<SupportRead> PrimaryReads;
        public final List<SupportRead> DuplicateReads; // identical reads from another assembly
        public int FragmentLength;
        public boolean Processed;

        public FragmentReads(final SupportRead read)
        {
            PrimaryReads = Lists.newArrayList(read);
            DuplicateReads = Lists.newArrayList();
            FragmentLength = -1;
            Processed = false;
        }

        public void addRead(final SupportRead read)
        {
            if(PrimaryReads.stream().anyMatch(x -> x.flags() == read.flags()))
                DuplicateReads.add(read);
            else
                PrimaryReads.add(read);
        }

        public void setReadDetails(final SupportRead read)
        {
            read.setInferredFragmentLength(FragmentLength);

            SupportRead matchingRead = PrimaryReads.stream().filter(x -> x.flags() == read.flags()).findFirst().orElse(null);

            if(matchingRead != null)
                read.setBreakendSupportType(matchingRead.breakendSupportType());
        }

        public void setFragmentLength()
        {
            FragmentLength = PrimaryReads.stream().mapToInt(x -> x.inferredFragmentLength()).max().orElse(-1);
            PrimaryReads.forEach(x -> x.setInferredFragmentLength(FragmentLength));
            DuplicateReads.forEach(x -> x.setInferredFragmentLength(FragmentLength));
        }

        public SupportRead findSpecificRead(boolean getFirst, boolean useSupplementaries)
        {
            SupportRead read = PrimaryReads.stream().filter(x -> x.firstInPair() == getFirst && !x.isSupplementary()).findFirst().orElse(null);

            if(read != null)
                return read;

            if(!useSupplementaries)
                return null;

            return PrimaryReads.stream().filter(x -> x.firstInPair() == getFirst).findFirst().orElse(null);
        }

        public String toString()
        {
            return format("%s: reads(%d) duplicates(%d) fragmentLength(%d) processed(%s)",
                    PrimaryReads.get(0).id(), PrimaryReads.size(), DuplicateReads.size(), FragmentLength, Processed);
        }
    }

    public void allocateBreakendSupport()
    {
        Map<String,FragmentReads> fragmentMap = Maps.newHashMap();

        List<JunctionAssembly> assemblies;

        if(mAssemblyAlignment.phaseSet() != null)
        {
            assemblies = Lists.newArrayList(mAssemblyAlignment.phaseSet().allAssemblies());

            for(PhaseSet mergedPhaseSet : mAssemblyAlignment.phaseSet().mergedPhaseSets())
            {
                assemblies.addAll(mergedPhaseSet.allAssemblies());
            }
        }
        else
        {
            assemblies = mAssemblyAlignment.assemblies();
        }

        for(JunctionAssembly assembly : assemblies)
        {
            for(SupportRead read : assembly.support())
            {
                // only check split status and fragment length from the primary pair, but take a supplementary where the primary wasn't captured
                // the expectation is that the primary also forms a junction assembly and that the two of them have formed a link,
                // so the primary's support for the link will be captured
                FragmentReads fragmentReads = fragmentMap.get(read.id());

                if(fragmentReads == null)
                {
                    fragmentReads = new FragmentReads(read);
                    fragmentMap.put(read.id(), fragmentReads);
                    continue;
                }

                // if the group has already been processed, then just fill in the read details for this read
                if(fragmentReads.Processed)
                {
                    fragmentReads.setReadDetails(read);
                    continue;
                }

                fragmentReads.addRead(read);

                SupportRead firstRead = fragmentReads.findSpecificRead(true, false);
                SupportRead secondRead = fragmentReads.findSpecificRead(false, false);

                if(firstRead == null || secondRead == null)
                    continue;

                fragmentReads.Processed = true;

                // find associated breakends
                processSupportReads(firstRead, secondRead);

                // apply to all
                fragmentReads.setFragmentLength();
            }
        }

        // handle incomplete fragment reads
        for(FragmentReads fragmentReads : fragmentMap.values())
        {
            if(fragmentReads.Processed)
                continue;

            SupportRead firstRead = fragmentReads.findSpecificRead(true, true);
            SupportRead secondRead = fragmentReads.findSpecificRead(false, true);

            processSupportReads(firstRead, secondRead);

            fragmentReads.setFragmentLength();
        }
    }

    private void processSupportReads(final SupportRead firstRead, final SupportRead secondRead)
    {
        List<ReadBreakendMatch> firstBreakendMatches = firstRead != null ? findReadBreakendMatch(firstRead) : Collections.emptyList();
        List<ReadBreakendMatch> secondBreakendMatches = secondRead != null ? findReadBreakendMatch(secondRead) : Collections.emptyList();

        if(!firstBreakendMatches.isEmpty() && !secondBreakendMatches.isEmpty())
        {
            processCompleteFragment(firstRead, secondRead, firstBreakendMatches, secondBreakendMatches);
        }
        else if(!firstBreakendMatches.isEmpty())
        {
            processSoloRead(firstRead, firstBreakendMatches);
        }
        else if(!secondBreakendMatches.isEmpty())
        {
            processSoloRead(secondRead, secondBreakendMatches);
        }
    }

    private void processCompleteFragment(
            final SupportRead firstRead, final SupportRead secondRead,
            final List<ReadBreakendMatch> firstBreakendMatches, final List<ReadBreakendMatch> secondBreakendMatches)
    {
        int lowerIndex = min(firstRead.fullAssemblyIndexStart(), secondRead.fullAssemblyIndexStart());
        int upperIndex = max(firstRead.fullAssemblyIndexEnd(), secondRead.fullAssemblyIndexEnd());

        int fragmentLength = upperIndex - lowerIndex + 1;
        firstRead.setInferredFragmentLength(fragmentLength);
        secondRead.setInferredFragmentLength(fragmentLength);

        int forwardReads = 0;
        int reverseReads = 0;

        if(firstRead.type() == SupportType.JUNCTION)
        {
            if(firstRead.orientation().isForward())
                ++forwardReads;
            else
                ++reverseReads;
        }

        if(secondRead.type() == SupportType.JUNCTION)
        {
            if(secondRead.orientation().isForward())
                ++forwardReads;
            else
                ++reverseReads;
        }

        Set<Breakend> breakends = Sets.newHashSet();

        // add each breakend and its pair only once to ensure they are both updated with the same split/discordant status
        addUniqueBreakends(breakends, firstBreakendMatches);
        addUniqueBreakends(breakends, secondBreakendMatches);

        boolean isShortIndel = breakends.stream().allMatch(x -> x.isShortLocalDelDupIns());

        boolean isSplitSupport = firstBreakendMatches.stream().anyMatch(x -> x.IsSplit)
                || secondBreakendMatches.stream().anyMatch(x -> x.IsSplit);

        if(!isSplitSupport && isShortIndel)
            return;

        firstRead.setBreakendSupportType(isSplitSupport ? SupportType.JUNCTION : SupportType.DISCORDANT);
        secondRead.setBreakendSupportType(isSplitSupport ? SupportType.JUNCTION : SupportType.DISCORDANT);

        for(Breakend breakend : breakends)
        {
            breakend.updateBreakendSupport(firstRead.sampleIndex(), isSplitSupport, forwardReads, reverseReads);
            breakend.addInferredFragmentLength(fragmentLength, true);

            if(!breakend.isSingle())
            {
                Breakend otherBreakend = breakend.otherBreakend();
                otherBreakend.updateBreakendSupport(firstRead.sampleIndex(), isSplitSupport, forwardReads, reverseReads);
                otherBreakend.addInferredFragmentLength(fragmentLength, true);
            }
        }
    }

    private static void addUniqueBreakends(final Set<Breakend> breakends, final List<ReadBreakendMatch> readBreakendMatches)
    {
        for(ReadBreakendMatch breakendMatch : readBreakendMatches)
        {
            if(breakends.contains(breakendMatch.Breakend))
                continue;

            if(!breakendMatch.Breakend.isSingle() && breakends.contains(breakendMatch.Breakend.otherBreakend()))
                continue;

            breakends.add(breakendMatch.Breakend);
        }
    }

    private void processSoloRead(final SupportRead read, final List<ReadBreakendMatch> readBreakendMatches)
    {
        int forwardReads = 0;
        int reverseReads = 0;

        if(read.type() == SupportType.JUNCTION)
        {
            if(read.orientation().isForward())
                ++forwardReads;
            else
                ++reverseReads;
        }

        Set<Breakend> breakends = Sets.newHashSet();
        addUniqueBreakends(breakends, readBreakendMatches);

        int inferredFragmentLength = -1;

        // since these reads are missing a mate, manually calculate their fragment length factoring in the simple SV type if they are local
        boolean isShortIndel = breakends.stream().allMatch(x -> x.isShortLocalDelDupIns());
        boolean setValidFragmentLength = false;
        StructuralVariantType svType = null;

        int indelLength = 0;

        if(isLocalIndel())
        {
            indelLength = mAssemblyAlignment.phaseSet().assemblyLinks().get(0).length();
            svType = mAssemblyAlignment.phaseSet().assemblyLinks().get(0).svType();
        }
        else if(!read.isDiscordant() && breakends.size() <= 2)
        {
            indelLength = breakends.iterator().next().svLength();
            svType = breakends.iterator().next().svType();
        }

        if(svType != null && isIndel(svType) && indelLength != 0 && read.insertSize() > 0)
        {
            inferredFragmentLength = calcIndelSoloReadFragmentLength(read, svType, indelLength);

            setValidFragmentLength = inferredFragmentLength <= MAX_OBSERVED_CONCORDANT_FRAG_LENGTH;

            if(setValidFragmentLength)
                read.setInferredFragmentLength(inferredFragmentLength);
        }

        boolean isSplitSupport = readBreakendMatches.stream().anyMatch(x -> x.IsSplit);

        if(!isSplitSupport && isShortIndel)
            return;

        read.setBreakendSupportType(isSplitSupport ? SupportType.JUNCTION : SupportType.DISCORDANT);

        for(Breakend breakend : breakends)
        {
            breakend.updateBreakendSupport(read.sampleIndex(), isSplitSupport, forwardReads, reverseReads);
            breakend.addInferredFragmentLength(inferredFragmentLength, setValidFragmentLength);

            if(!breakend.isSingle())
            {
                breakend.otherBreakend().updateBreakendSupport(read.sampleIndex(), isSplitSupport, forwardReads, reverseReads);
                breakend.otherBreakend().addInferredFragmentLength(inferredFragmentLength, setValidFragmentLength);
            }
        }
    }

    private int calcIndelSoloReadFragmentLength(final SupportRead read, final StructuralVariantType svType, int indelLength)
    {
        int fragmentLength = abs(read.insertSize());

        if(mAssemblyAlignment.assemblies().stream().allMatch(x -> x.outcome() == LOCAL_INDEL))
        {
            if(read.isPairedRead() && svType == DEL)
            {
                JunctionAssembly refAssembly = mAssemblyAlignment.assemblies().stream().filter(x -> x.supportCount() == 0).findFirst().orElse(null);

                if(refAssembly == null)
                    return 0;

                if((refAssembly.isForwardJunction() && read.orientation().isReverse())
                || (refAssembly.isReverseJunction() && read.orientation().isForward()))
                {
                    // the delete falls within the two reads, so subtract its length as usual
                    fragmentLength -= indelLength;
                    return fragmentLength;
                }
            }

            // soft-clip bases which were realigned locally extend the fragment length beyond the other junction if the mate is not across
            // from the junction
            JunctionAssembly juncAssembly = mAssemblyAlignment.assemblies().stream().filter(x -> x.supportCount() > 0).findFirst().orElse(null);

            if(juncAssembly.isForwardJunction())
                fragmentLength += read.rightClipLength();
            else
                fragmentLength += read.leftClipLength();
        }
        else
        {
            if(svType == DEL)
            {
                // remove the deleted length from the aligned fragment length
                fragmentLength -= indelLength;
            }
            else
            {
                // add any inserted bases
                fragmentLength += indelLength;
            }
        }

        return fragmentLength;
    }

    private boolean isLocalIndel()
    {
        if(mAssemblyAlignment.assemblies().size() != 2)
            return false;

        if(mAssemblyAlignment.phaseSet() == null || mAssemblyAlignment.phaseSet().assemblyLinks().size() != 1)
            return false;

        if(mAssemblyAlignment.assemblies().stream().allMatch(x -> x.indel()))
            return true;

        if(mAssemblyAlignment.assemblies().stream().allMatch(x -> x.outcome() == LOCAL_INDEL))
            return true;

        // otherwise check the characteristics of the link
        StructuralVariantType svType = mAssemblyAlignment.phaseSet().assemblyLinks().get(0).svType();
        int svLength = mAssemblyAlignment.phaseSet().assemblyLinks().get(0).length();

        return isShortLocalDelDupIns(svType, svLength);
    }

    private class ReadBreakendMatch
    {
        public final boolean IsSplit;
        public final Breakend Breakend;

        public ReadBreakendMatch(final Breakend breakend, final boolean isSplit)
        {
            IsSplit = isSplit;
            Breakend = breakend;
        }

        public String toString() { return format("%s breakend(%s)", IsSplit ? "split" : "disc", Breakend); }
    }

    private List<ReadBreakendMatch> findReadBreakendMatch(final SupportRead read)
    {
        Breakend closestDiscordantBreakend = null;
        int closestDiscJunctionDistance = INVALID_DISCORANT_DISTANCE;

        List<ReadBreakendMatch> breakendMatches = Lists.newArrayListWithCapacity(2);

        for(Breakend breakend : mAssemblyAlignment.breakends())
        {
            if(readSpansJunction(breakend, read))
            {
                breakendMatches.add(new ReadBreakendMatch(breakend, true));
                continue;
            }

            int discordantDistance = readDiscordantBreakendDistance(breakend, read);

            if(discordantDistance != INVALID_DISCORANT_DISTANCE)
            {
                if(closestDiscJunctionDistance == INVALID_DISCORANT_DISTANCE || discordantDistance < closestDiscJunctionDistance)
                {
                    closestDiscJunctionDistance = discordantDistance;
                    closestDiscordantBreakend = breakend;
                }
            }
        }

        if(closestDiscordantBreakend != null)
            breakendMatches.add(new ReadBreakendMatch(closestDiscordantBreakend, false));

        return breakendMatches;
    }

    private boolean readSpansJunction(final Breakend breakend, final SupportRead read)
    {
        // first an aligned junction read
        if(read.unclippedStart() < breakend.Position && read.unclippedEnd() > breakend.Position)
            return true;

        // next a misaligned junction read - crossing the segment boundary
        int readSeqIndexStart = read.fullAssemblyIndexStart();
        int readSeqIndexEnd = read.fullAssemblyIndexEnd();

        int fullSequenceEndIndex = mAssemblyAlignment.fullSequenceLength() - 1;

        // look for a read crossing any segment boundary
        for(BreakendSegment segment : breakend.segments())
        {
            int segmentSeqStart = segment.Alignment.sequenceStart();

            if(segmentSeqStart > 0 && readSeqIndexStart < segmentSeqStart && readSeqIndexEnd > segmentSeqStart)
                return true;

            int segmentSeqEnd = segment.Alignment.sequenceEnd();

            if(segmentSeqEnd < fullSequenceEndIndex && readSeqIndexStart < segmentSeqEnd && readSeqIndexEnd > segmentSeqEnd)
                return true;

            if(segment.indelSeqenceIndices() != null)
            {
                int[] indelSequenceIndices = segment.indelSeqenceIndices();

                if(readSeqIndexStart < indelSequenceIndices[0] && readSeqIndexEnd > indelSequenceIndices[0])
                    return true;

                if(readSeqIndexStart < indelSequenceIndices[1] && readSeqIndexEnd > indelSequenceIndices[1])
                    return true;
            }
        }

        /*
        for(BreakendSegment segment : breakend.segments())
        {
            int segmentBreakendIndex = segment.SegmentOrient.isForward() ? segment.Alignment.sequenceEnd() : segment.Alignment.sequenceStart();

            if(readSeqIndexStart < segmentBreakendIndex && readSeqIndexEnd > segmentBreakendIndex)
                return true;
        }
        */

        return false;
    }

    private static final int INVALID_DISCORANT_DISTANCE = -1;

    private static int readDiscordantBreakendDistance(final Breakend breakend, final SupportRead read)
    {
        // first check for an aligned discordant read
        if(breakend.Chromosome.equals(read.chromosome()) && read.orientation() == breakend.Orient)
        {
            if(breakend.Orient.isForward())
            {
                int maxPosition = max(breakend.maxPosition(), breakend.Position);

                if(read.alignmentEnd() <= maxPosition && read.alignmentStart() >= breakend.Position - DEFAULT_MAX_CONCORDANT_FRAG_LENGTH)
                    return maxPosition - read.alignmentEnd();
            }
            else
            {
                int minPosition = min(breakend.minPosition(), breakend.Position);

                if(read.alignmentStart() >= minPosition && read.alignmentEnd() <= breakend.Position + DEFAULT_MAX_CONCORDANT_FRAG_LENGTH)
                    return read.alignmentStart() - minPosition;
            }
        }

        int readSeqIndexStart = read.fullAssemblyIndexStart();
        int readSeqIndexEnd = read.fullAssemblyIndexEnd();

        for(BreakendSegment segment : breakend.segments())
        {
            // read must lie within the segment for it be applicable
            if(!positionsWithin(readSeqIndexStart, readSeqIndexEnd, segment.Alignment.sequenceStart(), segment.Alignment.sequenceEnd()))
                continue;

            int discordantDistance;

            if(read.fullAssemblyOrientation().isForward())
            {
                // calculate the distance between the read's end and the segment end
                discordantDistance = segment.Alignment.sequenceEnd() - readSeqIndexEnd;
            }
            else
            {
                discordantDistance = readSeqIndexStart - segment.Alignment.sequenceStart();
            }

            if(discordantDistance > 0 && discordantDistance < DEFAULT_MAX_CONCORDANT_FRAG_LENGTH)
                return discordantDistance;
        }

        /*
        for(BreakendSegment segment : breakend.segments())
        {
            if(read.fullAssemblyOrientation() != segment.SegmentOrient)
                continue;

            if(segment.SegmentOrient.isForward())
            {
                int segmentBreakendIndex = segment.Alignment.sequenceEnd();

                if(readSeqIndexEnd <= segmentBreakendIndex && readSeqIndexEnd >= segmentBreakendIndex - DEFAULT_DISCORDANT_FRAGMENT_LENGTH)
                    return segmentBreakendIndex - readSeqIndexEnd;
            }
            else
            {
                int segmentBreakendIndex = segment.Alignment.sequenceStart();

                if(readSeqIndexStart >= segmentBreakendIndex && readSeqIndexStart <= segmentBreakendIndex + DEFAULT_DISCORDANT_FRAGMENT_LENGTH)
                    return readSeqIndexStart - segmentBreakendIndex;
            }
        }
        */

        return INVALID_DISCORANT_DISTANCE;
    }
}

