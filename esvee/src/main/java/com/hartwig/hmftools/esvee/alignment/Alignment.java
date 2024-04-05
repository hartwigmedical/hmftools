package com.hartwig.hmftools.esvee.alignment;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PROXIMATE_DEL_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.types.ThreadTask.mergePerfCounters;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;

import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome;
import com.hartwig.hmftools.esvee.assembly.types.AssemblySupport;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;
import com.hartwig.hmftools.esvee.assembly.types.ThreadTask;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;

public class Alignment
{
    private final AssemblyConfig mConfig;

    private final AlignmentWriter mWriter;
    private final Aligner mAligner;
    private final List<Thread> mThreadTasks;

    public Alignment(final AssemblyConfig config, final Aligner aligner)
    {
        mConfig = config;
        mAligner = aligner;
        mWriter = new AlignmentWriter(mConfig);
        mThreadTasks = Lists.newArrayList();
    }

    private static final int MIN_ALIGN_LENGTH = MIN_VARIANT_LENGTH * 2;
    private static final int MIN_SUPPORT_COUNT = 4;

    public static boolean alignJunctionAssembly(final JunctionAssembly assembly)
    {
        // apply filters on what to bother aligning

        if(assembly.refBaseTrimLength() < MIN_ALIGN_LENGTH)
            return false;

        if(assembly.extensionLength() < MIN_ALIGN_LENGTH)
            return false;

        if(assembly.outcome() == AssemblyOutcome.DUP_BRANCHED
        || assembly.outcome() == AssemblyOutcome.DUP_SPLIT
        || assembly.outcome() == AssemblyOutcome.SECONDARY)
        {
            // since identical to or associated with other links
            return false;
        }

        if(assembly.supportCount() < MIN_SUPPORT_COUNT)
            return false;

        return true;
    }

    public static boolean alignAssemblyLink(final AssemblyLink assemblyLink)
    {
        if(assemblyLink.type() != LinkType.SPLIT)
            return false;

        if(assemblyLink.svType() == DEL || assemblyLink.svType() == DUP)
        {
            if(assemblyLink.length() < PROXIMATE_DEL_LENGTH)
                return false;
        }

        int combinedSequenceLength = assemblyLink.first().refBaseTrimLength() + assemblyLink.second().refBaseTrimLength()
                + assemblyLink.insertedBases().length() - assemblyLink.overlapBases().length();

        if(combinedSequenceLength < MIN_ALIGN_LENGTH)
            return false;

        Set<String> uniqueFrags = Sets.newHashSet();

        for(int i = 0; i <= 1; ++i)
        {
            JunctionAssembly assembly = (i == 0) ? assemblyLink.first() : assemblyLink.second();

            for(AssemblySupport support : assembly.support())
            {
                if(support.type() == SupportType.JUNCTION_MATE)
                    continue;

                if(uniqueFrags.contains(support.read().id()))
                    continue;

                uniqueFrags.add(support.read().id());

                if(uniqueFrags.size() >= MIN_SUPPORT_COUNT)
                    return true;
            }
        }

        return false;
    }

    public void run(final List<AssemblyAlignment> assemblyAlignments, final List<PerformanceCounter> perfCounters)
    {
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
            String fullSequence = assemblyAlignment.fullSequence();

            List<BwaMemAlignment> alignmentResults = mAligner.alignSequence(fullSequence.getBytes());

            if(alignmentResults.isEmpty())
                return;

            AlignmentWriter.writeAssemblyAlignment(mWriter.alignmentWriter(), assemblyAlignment, fullSequence, alignmentResults);

            AlignmentWriter.writeAlignmentDetails(mWriter.alignmentDetailsWriter(), assemblyAlignment, alignmentResults);
        }
    }
}
