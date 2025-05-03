package com.hartwig.hmftools.esvee.assembly.alignment;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.perf.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ALIGNMENT_REQUERY_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.hasDominantIndelReadAssembly;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.NO_LINK;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.UNSET;
import static com.hartwig.hmftools.esvee.assembly.types.ThreadTask.mergePerfCounters;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.perf.PerformanceCounter;
import com.hartwig.hmftools.esvee.assembly.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.output.AlignmentWriter;
import com.hartwig.hmftools.esvee.assembly.output.WriteType;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.ThreadTask;
import com.hartwig.hmftools.common.perf.TaskQueue;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;

import htsjdk.samtools.CigarOperator;

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

        if(!assembly.junction().indelBased() && (assembly.outcome() == UNSET || assembly.outcome() == NO_LINK))
        {
            if(hasDominantIndelReadAssembly(assembly))
                return true;
        }

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
        List<AssemblerAlignerTask> alignerTasks = Lists.newArrayList();

        int taskCount = min(mConfig.Threads, assemblyAlignments.size());

        for(int i = 0; i < taskCount; ++i)
        {
            AssemblerAlignerTask assemblerAlignerTask = new AssemblerAlignerTask(taskQueue);
            alignerTasks.add(assemblerAlignerTask);
            threadTasks.add(assemblerAlignerTask);
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

    private class AssemblerAlignerTask extends ThreadTask
    {
        private final TaskQueue mAssemblyAlignments;
        private int mRequeriedSuppCount;
        private int mRequeriedSoftClipCount;

        public AssemblerAlignerTask(final TaskQueue assemblyAlignments)
        {
            super("AssemblerAlignment");
            mAssemblyAlignments = assemblyAlignments;
            mRequeriedSuppCount = 0;
            mRequeriedSoftClipCount = 0;
        }

        public int requeriedSuppCount() { return mRequeriedSuppCount; }
        public int requeriedSoftClipCount() { return mRequeriedSoftClipCount; }

        @Override
        public void run()
        {
            while(true)
            {
                try
                {
                    mPerfCounter.start();

                    AssemblyAlignment assemblyAlignment = (AssemblyAlignment)mAssemblyAlignments.removeItem();

                    processAssembly(assemblyAlignment);

                    stopCheckLog(
                            format("alignment count(%d) assemblies(%s)", assemblyAlignment.assemblies().size(), assemblyAlignment.info()),
                            mConfig.PerfLogTime);
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
            if(!assemblyAlignment.isValid())
            {
                SV_LOGGER.warn("assembly alignment({}) invalid, skipping", assemblyAlignment);
                return;
            }

            if(assemblyAlignment.isMerged())
            {
                writeAssemblyData(assemblyAlignment, Collections.emptyList(), Collections.emptyList());
                return;
            }

            List<BwaMemAlignment> bwaAlignments = mAligner.alignSequence(assemblyAlignment.fullSequence().getBytes());

            List<AlignData> alignments = bwaAlignments.stream()
                    .map(x -> AlignData.from(x, mConfig.RefGenVersion))
                    .filter(x -> x != null).collect(Collectors.toList());

            List<AlignData> requeriedAlignments = Lists.newArrayList();

            alignments = requerySupplementaryAlignments(assemblyAlignment, alignments, requeriedAlignments);

            alignments = requerySoftClipAlignments(assemblyAlignment, alignments, requeriedAlignments);

            processAlignmentResults(assemblyAlignment, alignments);

            AlignmentFragments alignmentFragments = new AlignmentFragments(assemblyAlignment, mConfig.combinedSampleIds());
            alignmentFragments.allocateBreakendSupport();

            writeAssemblyData(assemblyAlignment, alignments, requeriedAlignments);
        }

        private List<AlignData> requerySoftClipAlignments(
                final AssemblyAlignment assemblyAlignment, final List<AlignData> alignments, final List<AlignData> requeriedAlignments)
        {
            // re-align supplementaries to get a more reliable map quality
            if(alignments.size() != 1)
                return alignments;

            AlignData alignment = alignments.get(0);

            int softClipLength = max(alignment.leftSoftClipLength(), alignment.rightSoftClipLength());

            if(softClipLength <= ALIGNMENT_REQUERY_SOFT_CLIP_LENGTH)
                return alignments;

            String softClipBases = "";
            String newCigar = "";
            boolean isLeftClip = false;

            if(alignment.leftSoftClipLength() > alignment.rightSoftClipLength())
            {
                softClipBases = assemblyAlignment.fullSequence().substring(0, alignment.leftSoftClipLength());
                newCigar = alignment.cigar().substring(alignment.cigar().indexOf(CigarOperator.S.toString()) + 1);
                isLeftClip = true;
            }
            else
            {
                int fullLength = assemblyAlignment.fullSequenceLength();
                softClipBases = assemblyAlignment.fullSequence().substring(fullLength - alignment.rightSoftClipLength());
                newCigar = alignment.cigar().substring(0, alignment.cigar().lastIndexOf(CigarOperator.M.toString()) + 1);
            }

            ++mRequeriedSoftClipCount;

            List<BwaMemAlignment> requeryBwaAlignments = mAligner.alignSequence(softClipBases.getBytes());

            List<AlignData> newAlignments = requeryBwaAlignments.stream()
                    .map(x -> AlignData.from(x, mConfig.RefGenVersion))
                    .filter(x -> x != null).collect(Collectors.toList());

            if(newAlignments.isEmpty())
                return alignments;

            requeriedAlignments.add(alignment);

            SV_LOGGER.trace("assembly({}) requeried single alignment({}) with long soft-clip", assemblyAlignment, alignment);

            AlignData refAlignment = new AlignData(
                    alignment.refLocation(), alignment.rawSequenceStart(), alignment.rawSequenceEnd(),
                    alignment.mapQual(), alignment.score(), alignment.flags(), newCigar, alignment.nMatches(), alignment.xaTag(), alignment.mdTag());

            if(isLeftClip)
                newAlignments.add(refAlignment);
            else
                newAlignments.add(0, refAlignment);

            return newAlignments;
        }

        private List<AlignData> requerySupplementaryAlignments(
                final AssemblyAlignment assemblyAlignment, final List<AlignData> alignments, final List<AlignData> requeriedAlignments)
        {
            // re-align supplementaries to get a more reliable map quality
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

                requeriedAlignments.add(alignData);
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
                        rqAlignment.refLocation(),
                        rqAlignment.rawSequenceStart(),
                        rqAlignment.rawSequenceEnd(),
                        rqAlignment.mapQual(), rqAlignment.score(), rqAlignment.flags(), rqAlignment.cigar(), rqAlignment.nMatches(),
                        rqAlignment.xaTag(), rqAlignment.mdTag());

                // restore values to be in terms of the original sequence
                int rqSeqOffsetStart = rqAlignment.sequenceStart();
                int adjSequenceStart = alignData.sequenceStart() + rqSeqOffsetStart;
                int rqSeqOffsetEnd = alignmentSequence.length() - 1 - rqAlignment.sequenceEnd();
                int adjSequenceEnd = alignData.sequenceEnd() - rqSeqOffsetEnd;
                convertedAlignment.setRequeriedSequenceCoords(adjSequenceStart, adjSequenceEnd);

                convertedAlignment.setSoftClipLengths(
                        convertedAlignment.leftSoftClipLength() + alignData.leftSoftClipLength(),
                        convertedAlignment.rightSoftClipLength() + alignData.rightSoftClipLength());

                convertedAlignments.add(convertedAlignment);

                alignData.markDroppedOnRequery();
            }

            return convertedAlignments;
        }

        private void processAlignmentResults(final AssemblyAlignment assemblyAlignment, final List<AlignData> alignments)
        {
            if(alignments.isEmpty())
                return;

            BreakendBuilder breakendBuilder = new BreakendBuilder(mConfig.RefGenome, assemblyAlignment);
            breakendBuilder.formBreakends(alignments);
        }
    }

    private void writeAssemblyData(
            final AssemblyAlignment assemblyAlignment, final List<AlignData> alignments, final List<AlignData> requeriedAlignments)
    {
        if(!assemblyAlignment.isValid())
            return;

        if(mConfig.WriteTypes.contains(WriteType.PHASED_ASSEMBLY))
            AlignmentWriter.writeAssemblyAlignment(mWriter.alignmentWriter(), assemblyAlignment);

        if(mConfig.WriteTypes.contains(WriteType.ALIGNMENT))
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

            AlignmentWriter.writeAlignmentDetails(mWriter.alignmentDetailsWriter(), assemblyAlignment, alignmentsToWrite);
        }
    }
}
