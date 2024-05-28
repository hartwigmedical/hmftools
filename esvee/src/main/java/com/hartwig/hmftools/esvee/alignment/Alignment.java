package com.hartwig.hmftools.esvee.alignment;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.types.ThreadTask.mergePerfCounters;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;
import com.hartwig.hmftools.esvee.assembly.types.ThreadTask;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;

public class Alignment
{
    private final AssemblyConfig mConfig;

    private final Aligner mAligner;
    private final AlignmentWriter mWriter;
    private final AlignmentCache mAlignmentCache;

    public Alignment(final AssemblyConfig config, final Aligner aligner)
    {
        mConfig = config;
        mAligner = aligner;
        mWriter = new AlignmentWriter(mConfig);
        mAlignmentCache = new AlignmentCache(config.AlignmentFile);
    }

    public void close() { mWriter.close(); }

    public static boolean skipUnlinkedJunctionAssembly(final JunctionAssembly assembly)
    {
        // apply filters on what to run alignment on
        if(assembly.outcome() == AssemblyOutcome.DUP_BRANCHED
        || assembly.outcome() == AssemblyOutcome.SECONDARY
        || assembly.outcome() == AssemblyOutcome.REMOTE_REGION)
        {
            // since identical to or associated with other links
            return true;
        }

        return false;
    }

    public void run(final List<AssemblyAlignment> assemblyAlignments, final List<PerformanceCounter> perfCounters)
    {
        if(mAligner == null && !mAlignmentCache.enabled())
            return;

        int singleAssemblies = (int)assemblyAlignments.stream().filter(x -> x.svType() == SGL).count();
        int linkedAssemblies = assemblyAlignments.size() - singleAssemblies;

        SV_LOGGER.info("running alignment for {} assemblies, linked({}) single({})",
                assemblyAlignments.size(), linkedAssemblies, singleAssemblies);

        Queue<AssemblyAlignment> assemblyAlignmentQueue = new ConcurrentLinkedQueue<>();

        assemblyAlignments.forEach(x -> assemblyAlignmentQueue.add(x));

        List<Thread> threadTasks = new ArrayList<>();
        List<AssemblerAlignerTask> alignerTasks = Lists.newArrayList();

        int taskCount = min(mConfig.Threads, assemblyAlignments.size());

        for(int i = 0; i < taskCount; ++i)
        {
            AssemblerAlignerTask assemblerAlignerTask = new AssemblerAlignerTask(assemblyAlignmentQueue);
            alignerTasks.add(assemblerAlignerTask);
            threadTasks.add(assemblerAlignerTask);
        }

        if(!runThreadTasks(threadTasks))
            System.exit(1);

        SV_LOGGER.debug("requeried supp alignments({})", alignerTasks.stream().mapToInt(x ->x.requeriedSuppCount()).sum());

        SV_LOGGER.info("alignment complete");

        mergePerfCounters(perfCounters, alignerTasks.stream().collect(Collectors.toList()));
    }

    private class AssemblerAlignerTask extends ThreadTask
    {
        private final Queue<AssemblyAlignment> mAssemblyAlignments;
        private final int mAssemblyAlignmentCount;
        private int mRequeriedSuppCount;

        public AssemblerAlignerTask(final Queue<AssemblyAlignment> assemblyAlignments)
        {
            super("AssemblerAlignment");
            mAssemblyAlignments = assemblyAlignments;
            mAssemblyAlignmentCount = assemblyAlignments.size();
            mRequeriedSuppCount = 0;
        }

        private static final int LOG_COUNT = 10000;

        public int requeriedSuppCount() { return mRequeriedSuppCount; }

        @Override
        public void run()
        {
            while(true)
            {
                try
                {
                    int remainingCount = mAssemblyAlignments.size();
                    int processedCount = mAssemblyAlignmentCount - remainingCount;

                    mPerfCounter.start();

                    ++processedCount;

                    AssemblyAlignment assemblyAlignment = mAssemblyAlignments.remove();

                    processAssembly(assemblyAlignment);

                    if((processedCount % LOG_COUNT) == 0)
                    {
                        SV_LOGGER.debug("processed {} assembly alignments", processedCount);
                    }

                    stopCheckLog(assemblyAlignment.info(), mConfig.PerfLogTime);
                }
                catch(NoSuchElementException e)
                {
                    SV_LOGGER.trace("all alignment tasks complete");
                    break;
                }
                catch(Exception e)
                {
                    e.printStackTrace();
                    System.exit(1);
                }
            }
        }

        private void processAssembly(final AssemblyAlignment assemblyAlignment)
        {
            List<AlignData> alignments;

            if(mAlignmentCache.enabled())
            {
                alignments = mAlignmentCache.findAssemblyAlignments(assemblyAlignment.info());
            }
            else
            {
                List<BwaMemAlignment> bwaAlignments = mAligner.alignSequence(assemblyAlignment.fullSequence().getBytes());

                alignments = bwaAlignments.stream()
                        .map(x -> AlignData.from(x, mConfig.RefGenVersion))
                        .filter(x -> x != null).collect(Collectors.toList());

                alignments = requerySupplementaryAlignments(assemblyAlignment, alignments);
            }

            processAlignmentResults(assemblyAlignment, alignments);

            AlignmentWriter.writeAssemblyAlignment(mWriter.alignmentWriter(), assemblyAlignment, alignments);
            AlignmentWriter.writeAlignmentDetails(mWriter.alignmentDetailsWriter(), assemblyAlignment, alignments);
        }

        private List<AlignData> requerySupplementaryAlignments(final AssemblyAlignment assemblyAlignment, final List<AlignData> alignments)
        {
            // re-alignment supplementaries to get a more reliable map quality
            if(alignments.stream().noneMatch(x -> x.isSupplementary()))
                return alignments;

            List<AlignData> newAlignments = Lists.newArrayList();

            for(AlignData alignData : alignments)
            {
                if(!alignData.isSupplementary())
                {
                    newAlignments.add(alignData);
                    continue;
                }

                newAlignments.addAll(requeryAlignment(assemblyAlignment, alignData));
            }

            return newAlignments;
        }

        private List<AlignData> requeryAlignment(final AssemblyAlignment assemblyAlignment, final AlignData alignData)
        {
            ++mRequeriedSuppCount;

            String fullSequence = assemblyAlignment.fullSequence();

            alignData.setFullSequenceData(fullSequence, assemblyAlignment.fullSequenceLength());

            String alignmentSequence = fullSequence.substring(alignData.sequenceStart(), alignData.sequenceEnd() + 1);

            List<BwaMemAlignment> requeryBwaAlignments = mAligner.alignSequence(alignmentSequence.getBytes());

            List<AlignData> requeryAlignments = requeryBwaAlignments.stream()
                    .map(x -> AlignData.from(x, mConfig.RefGenVersion))
                    .filter(x -> x != null).collect(Collectors.toList());

            List<AlignData> convertedAlignments = Lists.newArrayList();

            for(AlignData rqAlignment : requeryAlignments)
            {
                rqAlignment.setFullSequenceData(alignmentSequence, alignmentSequence.length());

                // eg:
                // alignData = {AlignData@3240} "10:2543491-2543563 72S73M fwd seq(72-145 adj=72-144) score(58) flags(2048) mapQual(55 align=73 adj=73)"
                // rqAlignment = {AlignData@3246} "10:2543809-2543878 3S70M fwd seq(3-73 adj=3-72) score(65) flags(0) mapQual(17 align=70 adj=70)"

                AlignData convertedAlignment = new AlignData(
                        rqAlignment.RefLocation,
                        rqAlignment.rawSequenceStart(),
                        rqAlignment.rawSequenceEnd(),
                        rqAlignment.MapQual, rqAlignment.Score, rqAlignment.Flags, rqAlignment.Cigar, rqAlignment.NMatches,
                        rqAlignment.XaTag, rqAlignment.MdTag);

                // restore values to be in terms of the original sequence
                int rqSeqOffsetStart = rqAlignment.sequenceStart();
                int adjSequenceStart = alignData.sequenceStart() + rqSeqOffsetStart;
                int rqSeqOffsetEnd = alignmentSequence.length() - 1 - rqAlignment.sequenceEnd();
                int adjSequenceEnd = alignData.sequenceEnd() - rqSeqOffsetEnd;
                convertedAlignment.setRequeriedSequenceCoords(adjSequenceStart, adjSequenceEnd);

                convertedAlignments.add(convertedAlignment);
            }

            return convertedAlignments;
        }

        private void processAlignmentResults(
                final AssemblyAlignment assemblyAlignment, final List<AlignData> alignments)
        {
            BreakendBuilder breakendBuilder = new BreakendBuilder(mConfig.RefGenome, assemblyAlignment);
            breakendBuilder.formBreakends(alignments);

            for(JunctionAssembly assembly : assemblyAlignment.assemblies())
            {
                boolean matched = false;

                for(Breakend breakend : assemblyAlignment.breakends())
                {
                    if(breakend.matches(assembly.junction().Chromosome, assembly.junction().Position, assembly.junction().Orient))
                    {
                        assembly.setAlignmentOutcome(AlignmentOutcome.MATCH);
                        matched = true;
                        break;
                    }
                }

                if(!matched)
                    assembly.setAlignmentOutcome(AlignmentOutcome.NO_MATCH);
            }

            if(assemblyAlignment.breakends().isEmpty())
                assemblyAlignment.assemblies().forEach(x -> x.setAlignmentOutcome(AlignmentOutcome.NO_RESULT));

            allocateSupport(assemblyAlignment);
        }

        private void allocateSupport(final AssemblyAlignment assemblyAlignment)
        {
            List<String> combinedSampleIds = mConfig.combinedSampleIds();

            // build up a map of read ID to the set of breakends it supports, and the top type of support (split then discordant)
            Map<String,BreakendFragmentSupport> fragmentSupportMap = Maps.newHashMap();

            for(Breakend breakend : assemblyAlignment.breakends())
            {
                // rather than use the genome position of a read vs the aligned breakend position, use its position in the assembly
                List<BreakendSupport> sampleSupport = breakend.sampleSupport();

                combinedSampleIds.forEach(x -> sampleSupport.add(new BreakendSupport()));

                for(JunctionAssembly assembly : assemblyAlignment.assemblies())
                {
                    for(SupportRead read : assembly.support())
                    {
                        if(read.type() == SupportType.JUNCTION_MATE) // any reason to count towards strand bias?
                            continue;

                        if(!breakend.Chromosome.equals(read.chromosome()))
                            continue;

                        boolean isSplitFragment = false;
                        boolean isDiscFragment = false;

                        BreakendSupport support = sampleSupport.get(read.sampleIndex());

                        if(breakend.readSpansJunction(read, false))
                        {
                            isSplitFragment = true;
                        }
                        else if(breakend.isRelatedDiscordantRead(read.alignmentStart(), read.alignmentEnd(), read.orientation()))
                        {
                            isDiscFragment = true;
                        }

                        if(!isSplitFragment && !isDiscFragment)
                            continue;

                        if(read.orientation().isForward())
                            ++support.ForwardReads;
                        else
                            ++support.ReverseReads;

                        BreakendFragmentSupport fragmentSupport = fragmentSupportMap.get(read.id());

                        if(fragmentSupport == null)
                        {
                            fragmentSupport = new BreakendFragmentSupport(read.sampleIndex(), isSplitFragment, breakend);
                            fragmentSupportMap.put(read.id(), fragmentSupport);
                        }
                        else
                        {
                            fragmentSupport.Breakends.add(breakend);
                            fragmentSupport.IsSplit |= isSplitFragment;
                        }
                    }
                }
            }

            // count fragments to both breakends if it is in either
            for(BreakendFragmentSupport fragmentSupport : fragmentSupportMap.values())
            {
                Set<Breakend> processed = Sets.newHashSet();

                for(Breakend breakend : fragmentSupport.Breakends)
                {
                    if(processed.contains(breakend) || (!breakend.isSingle() && processed.contains(breakend.otherBreakend())))
                        continue;

                    processed.add(breakend);

                    boolean allowDiscordantSupport = !breakend.isShortLocalDelDupIns();

                    BreakendSupport support = breakend.sampleSupport().get(fragmentSupport.SampleIndex);

                    if(fragmentSupport.IsSplit)
                        ++support.SplitFragments;
                    else if(allowDiscordantSupport)
                        ++support.DiscordantFragments;

                    if(breakend.isSingle())
                        continue;

                    Breakend otherBreakend = breakend.otherBreakend();

                    BreakendSupport otherSupport = otherBreakend.sampleSupport().get(fragmentSupport.SampleIndex);

                    processed.add(otherBreakend);

                    // assign each read preferably as split over discordant
                    if(fragmentSupport.IsSplit)
                        ++otherSupport.SplitFragments;
                    else if(allowDiscordantSupport)
                        ++otherSupport.DiscordantFragments;
                }
            }
        }
    }

    private class BreakendFragmentSupport
    {
        public final int SampleIndex;
        public boolean IsSplit;
        public final Set<Breakend> Breakends;

        public BreakendFragmentSupport(final int sampleIndex, final boolean isSplit, final Breakend breakend)
        {
            SampleIndex = sampleIndex;
            IsSplit = isSplit;
            Breakends = Sets.newHashSet(breakend);
        }

        public String toString() { return format("%d: %s breakends(%d)", SampleIndex, IsSplit ? "split" : "disc", Breakends.size()); }
    }

    public void calcAssemblyFragmentLengths(final List<List<AssemblyAlignment>> assemblyAlignmentGroups)
    {
        assemblyAlignmentGroups.forEach(x -> calcInferredFragmentLengths(x));
    }

    private void calcInferredFragmentLengths(final List<AssemblyAlignment> assemblyAlignments)
    {
        Map<String, BreakendFragmentData> fragmentSupportMap = Maps.newHashMap();

        for(AssemblyAlignment assemblyAlignment : assemblyAlignments)
        {
            for(Breakend breakend : assemblyAlignment.breakends())
            {
                for(JunctionAssembly assembly : assemblyAlignment.assemblies())
                {
                    for(SupportRead read : assembly.support())
                    {
                        if(read.isSupplementary())
                            continue;

                        if(!breakend.Chromosome.equals(read.chromosome()))
                            continue;

                        boolean isSplitFragment = false;
                        boolean isDiscFragment = false;

                        if(breakend.readSpansJunction(read, false))
                        {
                            isSplitFragment = true;
                        }
                        else if(breakend.isRelatedDiscordantRead(read.alignmentStart(), read.alignmentEnd(), read.orientation()))
                        {
                            isDiscFragment = true;
                        }

                        if(!isSplitFragment && !isDiscFragment)
                            continue;

                        BreakendFragmentData fragmentSupport = fragmentSupportMap.get(read.id());

                        if(fragmentSupport == null)
                        {
                            fragmentSupport = new BreakendFragmentData(read, breakend);
                            fragmentSupportMap.put(read.id(), fragmentSupport);
                        }
                        else
                        {
                            fragmentSupport.add(read, breakend);
                        }
                    }
                }
            }
        }

        Map<Breakend,LengthData> breakendFragmentLengths = Maps.newHashMap();

        for(Map.Entry<String,BreakendFragmentData> entry : fragmentSupportMap.entrySet())
        {
            BreakendFragmentData breakendFragmentData = entry.getValue();

            int inferredFragmentLength = calcInferredFragmentLength(breakendFragmentData);
            if(inferredFragmentLength != INVALID_FRAGMENT_LENGTH)
            {
                breakendFragmentData.Reads.forEach(x -> x.setInferredFragmentLength(inferredFragmentLength));

                for(Breakend breakend : breakendFragmentData.Breakends)
                {
                    LengthData lengthData = breakendFragmentLengths.get(breakend);

                    if(lengthData == null)
                    {
                        lengthData = new LengthData();
                        breakendFragmentLengths.put(breakend, lengthData);
                    }

                    ++lengthData.Count;
                    lengthData.TotalLength += inferredFragmentLength;
                }
            }
        }

        for(Map.Entry<Breakend,LengthData> entry : breakendFragmentLengths.entrySet())
        {
            Breakend breakend = entry.getKey();
            LengthData lengthData = entry.getValue();
            breakend.setAverageInferredFragmentLength(lengthData.averageLength());
        }
    }

    private class LengthData
    {
        public int Count;
        public int TotalLength;

        public LengthData()
        {
            Count = 0;
            TotalLength = 0;
        }

        public int averageLength() { return (int)round(TotalLength / (double)Count); }
    }

    public static final int INVALID_FRAGMENT_LENGTH = -1;

    private static Breakend findRelatedBreakend(final SupportRead read, final List<Breakend> breakends)
    {
        for(Breakend breakend : breakends)
        {
            if(breakend.readSpansJunction(read, false))
                return breakend;

            if(breakend.isRelatedDiscordantRead(read.alignmentStart(), read.alignmentEnd(), read.orientation()))
                return breakend;
        }

        return null;
    }

    public static int calcInferredFragmentLength(final BreakendFragmentData breakendFragmentData)
    {
        if(breakendFragmentData.Reads.size() != 2 || breakendFragmentData.Breakends.isEmpty())
            return INVALID_FRAGMENT_LENGTH;

        SupportRead read1 = breakendFragmentData.Reads.get(0);
        SupportRead read2 = breakendFragmentData.Reads.get(1);

        Breakend breakend1 = findRelatedBreakend(read1, breakendFragmentData.Breakends);
        Breakend breakend2 = findRelatedBreakend(read2, breakendFragmentData.Breakends);

        if(breakend1 == null || breakend2 == null)
            return INVALID_FRAGMENT_LENGTH;

        if(breakend1 == breakend2)
        {
            return max(read1.unclippedEnd(), read2.unclippedEnd()) - min(read1.unclippedStart(), read2.unclippedStart());
        }

        int junctionDistance = breakend1.Orient.isForward() ?
                breakend1.Position - read1.unclippedStart() : read1.unclippedEnd() - breakend1.Position;

        if(breakend1.otherBreakend() == breakend2)
        {
            junctionDistance += breakend2.Orient.isForward() ?
                    breakend2.Position - read2.unclippedStart() : read2.unclippedEnd() - breakend2.Position;

            // immediately across a junction, so take the distance from each
            return junctionDistance;
        }

        // factor in chained breakend links
        Breakend prevBreakend = breakend1;
        Breakend nextBreakend;
        boolean nextIsFacing = false;

        if(read1.orientation() == breakend1.Orient)
        {
            nextBreakend = breakend1.otherBreakend();
            nextIsFacing = false;
        }
        else
        {
            nextBreakend = !breakend1.facingBreakends().isEmpty() ? breakend1.facingBreakends().get(0) : null;
            nextIsFacing = true;
        }

        while(nextBreakend != null)
        {
            if(nextBreakend == breakend2)
            {
                if(breakend2.Orient.isForward())
                {
                    if(nextIsFacing)
                        junctionDistance += max(breakend2.Position, read2.unclippedEnd()) - prevBreakend.Position;
                    else
                        junctionDistance += breakend2.Position - read2.unclippedStart();
                }
                else
                {
                    if(nextIsFacing)
                        junctionDistance += prevBreakend.Position - min(breakend2.Position, read2.unclippedStart());
                    else
                        junctionDistance += read2.unclippedEnd() - breakend2.Position;
                }

                break;
            }

            if(nextIsFacing)
            {
                junctionDistance += abs(prevBreakend.Position - nextBreakend.Position);
                prevBreakend = nextBreakend;
                nextBreakend = nextBreakend.otherBreakend();
                nextIsFacing = false;
            }
            else
            {
                prevBreakend = nextBreakend;
                nextBreakend = !nextBreakend.facingBreakends().isEmpty() ? nextBreakend.facingBreakends().get(0) : null;
                nextIsFacing = true;
            }
        }

        return junctionDistance;
    }

    private class BreakendFragmentData
    {
        public final List<Breakend> Breakends;
        public final List<SupportRead> Reads;

        public BreakendFragmentData(final SupportRead read, final Breakend breakend)
        {
            Breakends = Lists.newArrayList(breakend);
            Reads = Lists.newArrayList(read);
        }

        public void add(final SupportRead read, final Breakend breakend)
        {
            if(Reads.stream().noneMatch(x -> x.flags() == read.flags()))
                Reads.add(read);

            if(!Breakends.contains(breakend))
                Breakends.add(breakend);
        }

        public String toString() { return format("breakends(%d) reads(%d)", Breakends.size(), Reads.size()); }
    }
}
