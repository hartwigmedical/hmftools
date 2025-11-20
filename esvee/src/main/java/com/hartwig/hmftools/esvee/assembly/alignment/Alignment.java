package com.hartwig.hmftools.esvee.assembly.alignment;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.perf.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.isWeakIndelBasedUnlinkedAssembly;
import static com.hartwig.hmftools.esvee.assembly.types.ThreadTask.mergePerfCounters;

import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.perf.PerformanceCounter;
import com.hartwig.hmftools.esvee.assembly.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.output.AlignmentWriter;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.common.perf.TaskQueue;
import com.hartwig.hmftools.esvee.common.WriteType;

public class Alignment
{
    private final AssemblyConfig mConfig;

    private final Aligner mAligner;
    private final AlignmentWriter mWriter;

    public Alignment(final AssemblyConfig config, final Aligner aligner)
    {
        mConfig = config;
        mAligner = aligner;
        mWriter = new AlignmentWriter(mConfig);
    }

    public void close() { mWriter.close(); }

    public static boolean skipUnlinkedJunctionAssembly(final JunctionAssembly assembly)
    {
        // apply filters on what to run alignment on
        if(assembly.outcome() == AssemblyOutcome.DUP_BRANCHED
        || assembly.outcome() == AssemblyOutcome.SECONDARY
        || assembly.outcome() == AssemblyOutcome.SUPP_ONLY)
        {
            // since identical to or associated with other links
            return true;
        }

        if(isWeakIndelBasedUnlinkedAssembly(assembly))
            return true;

        return false;
    }

    public void run(final List<AssemblyAlignment> assemblyAlignments, final List<PerformanceCounter> perfCounters)
    {
        if(mAligner == null)
            return;

        int singleAssemblies = (int) assemblyAlignments.stream().filter(x -> x.assemblies().size() == 1).count();
        int linkedAssemblies = assemblyAlignments.size() - singleAssemblies;

        SV_LOGGER.info("running alignment for {} assemblies, linked({}) single({})",
                assemblyAlignments.size(), linkedAssemblies, singleAssemblies);

        Queue<AssemblyAlignment> assemblyAlignmentQueue = new ConcurrentLinkedQueue<>();
        assemblyAlignments.forEach(x -> assemblyAlignmentQueue.add(x));

        TaskQueue taskQueue = new TaskQueue(assemblyAlignmentQueue, "assembly alignments", 10000);

        List<Thread> threadTasks = new ArrayList<>();
        List<AssemblyAligner> alignerTasks = Lists.newArrayList();

        int taskCount = min(mConfig.Threads, assemblyAlignments.size());

        for(int i = 0; i < taskCount; ++i)
        {
            AssemblyAligner alignerTask = new AssemblyAligner(mConfig, mAligner, mWriter, taskQueue);
            alignerTasks.add(alignerTask);
            threadTasks.add(alignerTask);
        }

        if(!runThreadTasks(threadTasks))
            System.exit(1);

        int requeriedSuppCount = alignerTasks.stream().mapToInt(x -> x.requeriedSuppCount()).sum();
        int requeriedSoftClipCount = alignerTasks.stream().mapToInt(x -> x.requeriedSoftClipCount()).sum();

        if(requeriedSuppCount > 0 || requeriedSoftClipCount > 0)
        {
            SV_LOGGER.debug("requeried supp alignments({}) soft-clips({})", requeriedSuppCount, requeriedSoftClipCount);
        }

        SV_LOGGER.info("alignment complete");

        mergePerfCounters(perfCounters, alignerTasks.stream().collect(Collectors.toList()));
    }

    protected static void writeAssemblyData(
            final AlignmentWriter writer, final AssemblyConfig config,
            final AssemblyAlignment assemblyAlignment, final List<AlignData> alignments, final List<AlignData> requeriedAlignments)
    {
        if(!assemblyAlignment.isValid())
            return;

        if(config.WriteTypes.contains(WriteType.PHASED_ASSEMBLY))
            AlignmentWriter.writePhasedAssembly(writer.alignmentWriter(), assemblyAlignment);

        if(config.WriteTypes.contains(WriteType.ALIGNMENT))
        {
            List<AlignData> alignmentsToWrite;

            if(!requeriedAlignments.isEmpty())
            {
                alignmentsToWrite = Lists.newArrayList(alignments);
                alignmentsToWrite.addAll(requeriedAlignments);
            }
            else
            {
                alignmentsToWrite = alignments;
            }

            AlignmentWriter.writeAlignmentDetails(writer.alignmentDetailsWriter(), assemblyAlignment, alignmentsToWrite);
        }
    }
}
