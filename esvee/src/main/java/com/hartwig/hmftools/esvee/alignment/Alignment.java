package com.hartwig.hmftools.esvee.alignment;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
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
        // apply filters on what to bother aligning
        if(assembly.outcome() == AssemblyOutcome.DUP_BRANCHED
        || assembly.outcome() == AssemblyOutcome.DUP_SPLIT
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

        SV_LOGGER.info("alignment complete");

        mergePerfCounters(perfCounters, alignerTasks.stream().collect(Collectors.toList()));
    }

    private class AssemblerAlignerTask extends ThreadTask
    {
        private final Queue<AssemblyAlignment> mAssemblyAlignments;
        private final int mAssemblyAlignmentCount;

        public AssemblerAlignerTask(final Queue<AssemblyAlignment> assemblyAlignments)
        {
            super("AssemblerAlignment");
            mAssemblyAlignments = assemblyAlignments;
            mAssemblyAlignmentCount = assemblyAlignments.size();
        }

        private static final int LOG_COUNT = 10000;

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

                    mPerfCounter.stop();
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
            }

            processAlignmentResults(assemblyAlignment, alignments);

            AlignmentWriter.writeAssemblyAlignment(mWriter.alignmentWriter(), assemblyAlignment, alignments);
            AlignmentWriter.writeAlignmentDetails(mWriter.alignmentDetailsWriter(), assemblyAlignment, alignments);
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
                    if(breakend.matches(assembly.junction().Chromosome, assembly.junction().Position, assembly.junction().Orientation))
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
                            isSplitFragment = true;
                        else
                            isDiscFragment = true;

                        if(!isSplitFragment && !isDiscFragment)
                            continue;

                        if(read.orientation() == POS_ORIENT)
                            ++support.ForwardReads;
                        else
                            ++support.ReverseReads;

                        BreakendFragmentSupport fragmentSupport = fragmentSupportMap.get(read.id());

                        if(fragmentSupport == null)
                        {
                            int inferredFragLength = assemblyAlignment.calcInferredFragmentLength(assembly, read);

                            fragmentSupport = new BreakendFragmentSupport(read.sampleIndex(), isSplitFragment, breakend, inferredFragLength);
                            fragmentSupportMap.put(read.id(), fragmentSupport);
                        }
                        else
                        {
                            fragmentSupport.addBreakend(breakend, isSplitFragment);
                        }
                    }
                }
            }

            int inferredFragmentCount = 0;
            int inferredFragmentTotal = 0;

            for(Map.Entry<String,BreakendFragmentSupport> entry : fragmentSupportMap.entrySet())
            {
                BreakendFragmentSupport fragmentSupport = entry.getValue();

                if(fragmentSupport.InferredFragmentLength > 0)
                {
                    ++inferredFragmentCount;
                    inferredFragmentTotal += fragmentSupport.InferredFragmentLength;
                }
            }

            int averageInferredFragLength = inferredFragmentCount > 0 ? (int)round(inferredFragmentTotal / (double)inferredFragmentCount) : 0;

            for(Breakend breakend : assemblyAlignment.breakends())
            {
                breakend.setAverageInferredFragmentLength(averageInferredFragLength);
            }

            // count fragments to both breakends if it is in either
            for(BreakendFragmentSupport fragmentSupport : fragmentSupportMap.values())
            {
                for(Breakend breakend : fragmentSupport.Breakends)
                {
                    boolean allowDiscordantSupport = !breakend.isShortLocalDelDupIns();

                    BreakendSupport support = breakend.sampleSupport().get(fragmentSupport.SampleIndex);

                    if(fragmentSupport.IsSplit)
                        ++support.SplitFragments;
                    else if(allowDiscordantSupport)
                        ++support.DiscordantFragments;

                    Breakend otherBreakend = breakend.otherBreakend();

                    if(otherBreakend != null && !fragmentSupport.Breakends.contains(otherBreakend))
                    {
                        BreakendSupport otherSupport = otherBreakend.sampleSupport().get(fragmentSupport.SampleIndex);

                        // assign each read preferably as split over discordant
                        if(fragmentSupport.IsSplit)
                            ++otherSupport.SplitFragments;
                        else if(allowDiscordantSupport)
                            ++otherSupport.DiscordantFragments;
                    }
                }
            }
        }
    }

    private class BreakendFragmentSupport
    {
        public final int SampleIndex;
        public boolean IsSplit;
        public final Set<Breakend> Breakends;
        public final int InferredFragmentLength;

        public BreakendFragmentSupport(final int sampleIndex, final boolean isSplit, final Breakend breakend, final int inferredFragLength)
        {
            SampleIndex = sampleIndex;
            IsSplit = isSplit;
            Breakends = Sets.newHashSet(breakend);
            InferredFragmentLength = inferredFragLength;
        }

        public void addBreakend(final Breakend breakend, boolean isSplit)
        {
            Breakends.add(breakend);
            IsSplit |= isSplit;
        }

        public String toString() { return format("%d: %s breakends(%d)", SampleIndex, IsSplit ? "split" : "disc", Breakends.size()); }
    }
}
