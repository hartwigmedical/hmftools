package com.hartwig.hmftools.esvee.old;

import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.TASK_LOG_COUNT;
import static com.hartwig.hmftools.esvee.read.ReadUtils.flipRead;

import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.common.RegionOfInterest;
import com.hartwig.hmftools.esvee.common.ThreadTask;
import com.hartwig.hmftools.esvee.read.BamReader;
import com.hartwig.hmftools.esvee.read.Read;

public class SupportScanner extends ThreadTask
{
    private final Queue<AlignedAssembly> mAlignedAssemblyQueue;
    private final int mAlignedAssemblyCount;

    private final SvConfig mConfig;
    private final BamReader mBamReader;

    // private final Counter mExtraSupportCounter; // previously took and incremented Counters.ExtraScannedSupport
    private int mExtraSupportCount;
    private final SupportChecker mSupportChecker;

    private final List<AlignedAssembly> mResults;

    public static List<SupportScanner> createThreadTasks(
            final List<AlignedAssembly> alignedAssemblies, final List<BamReader> bamReaders, final SvConfig config,
            final int taskCount, final List<Thread> threadTasks)
    {
        List<SupportScanner> supportScanners = Lists.newArrayList();

        Queue<AlignedAssembly> alignedAssemblyQueue = new ConcurrentLinkedQueue<>();
        alignedAssemblyQueue.addAll(alignedAssemblies);

        for(int i = 0; i < taskCount; ++i)
        {
            BamReader bamReader = bamReaders.get(i);

            SupportScanner supportScanner = new SupportScanner(config, bamReader, alignedAssemblyQueue);
            supportScanners.add(supportScanner);
            threadTasks.add(supportScanner);
        }

        SV_LOGGER.debug("splitting {} aligned assemblies for rescanning sets across {} threads", alignedAssemblies.size(), taskCount);

        return supportScanners;
    }

    public SupportScanner(final SvConfig config, final BamReader bamReader, final Queue<AlignedAssembly> alignedAssemblyQueue)
    {
        super("SupportScanner");

        mConfig = config;
        mAlignedAssemblyQueue = alignedAssemblyQueue;
        mAlignedAssemblyCount = alignedAssemblyQueue.size();
        mBamReader = bamReader;
        mExtraSupportCount = 0;
        mSupportChecker = new SupportChecker();

        mResults = Lists.newArrayList();
    }

    @Override
    public void run()
    {
        while(true)
        {
            try
            {
                int remainingCount = mAlignedAssemblyQueue.size();
                int processedCount = mAlignedAssemblyCount - remainingCount;

                AlignedAssembly alignedAssembly = mAlignedAssemblyQueue.remove();

                mPerfCounter.start();

                AlignedAssembly rescannedAssembly = rescanSupport(alignedAssembly);
                mResults.add(rescannedAssembly);

                // mPerfCounter.stop();
                stopCheckLog(alignedAssembly.toString(), mConfig.PerfLogTime);

                if(processedCount > 0 && (processedCount % TASK_LOG_COUNT) == 0)
                {
                    SV_LOGGER.info("processed support for {} aligned assemblies, remaining({})", processedCount, remainingCount);
                }
            }
            catch(NoSuchElementException e)
            {
                SV_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    public int addedReadCount() { return mExtraSupportCount; }

    /*
    public static List<AlignedAssembly> mergeRescannedAssemblies(final List<SupportScanner> supportScanners)
    {
        List<AlignedAssembly> combined = Lists.newArrayList();
        supportScanners.forEach(x -> combined.addAll(x.mResults));
        return combined;
    }
    */

    public AlignedAssembly rescanSupport(final AlignedAssembly assembly)
    {
        if(assembly.getAlignmentBlocks().size() < 2)
            return assembly; // No variant, don't care about support

        final List<RegionOfInterest> alignmentRegions = assembly.getAlignmentBlocks().stream()
                .filter(Alignment::isMapped)
                .map(block -> new RegionOfInterest(block.Chromosome, block.ReferenceStartPosition,
                        block.ReferenceStartPosition + block.Length))
                .collect(Collectors.toList());
        final List<RegionOfInterest> deduplicated = RegionOfInterest.tryMerge(alignmentRegions);

        /* FIXME:
        mBamReader.sliceBam();

        final List<Read> potentialSupport = deduplicated.stream()
                .flatMap(region -> mContext.SAMSource.findReadsContaining(region).stream())
                .distinct()
                .collect(Collectors.toList());

        for(Read read : potentialSupport)
            tryAddSupport(assembly, read);
        */

        // CHECK is there any need to create a new assembly object if just new support is being added?

        return assembly;
    }

    private void tryAddSupport(final AlignedAssembly assembly, final Read read)
    {
        if(assembly.containsSupport(read))
            return;

        if(assembly.tryAddSupport(mSupportChecker, read))
        {
            ++mExtraSupportCount;
        }
        else
        {
            final Read inverted = flipRead(read);
            if(assembly.tryAddSupport(mSupportChecker, inverted))
            {
                ++mExtraSupportCount;
            }
        }
    }
}
