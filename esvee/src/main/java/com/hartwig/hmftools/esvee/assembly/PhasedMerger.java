package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.TASK_LOG_COUNT;

import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.common.ThreadTask;
import com.hartwig.hmftools.esvee.sequence.ExtendedAssembly;

public class PhasedMerger extends ThreadTask
{
    private final Queue<Set<ExtendedAssembly>> mPrimaryPhaseSetQueue;
    private final int mPrimaryPhaseSetCount;

    private final AssemblyMerger mAssemblyMerger;
    private final SvConfig mConfig;

    private final List<Set<ExtendedAssembly>> mPhasedResults;

    public static List<PhasedMerger> createThreadTasks(
            final List<Set<ExtendedAssembly>> primaryPhaseSets, final SvConfig config, final int taskCount, final List<Thread> threadTasks)
    {
        List<PhasedMerger> phasedMergerTasks = Lists.newArrayList();

        Queue<Set<ExtendedAssembly>> primaryPhaseQueue = new ConcurrentLinkedQueue<>();
        primaryPhaseQueue.addAll(primaryPhaseSets);

        int primaryPhaseSetCount = primaryPhaseSets.size();

        for(int i = 0; i < taskCount; ++i)
        {
            PhasedMerger phasedMerger = new PhasedMerger(config, primaryPhaseQueue);
            phasedMergerTasks.add(phasedMerger);
            threadTasks.add(phasedMerger);
        }

        SV_LOGGER.debug("splitting {} primary phase sets across {} threads", primaryPhaseSetCount, taskCount);

        return phasedMergerTasks;
    }

    public PhasedMerger(final SvConfig config, final Queue<Set<ExtendedAssembly>> primaryPhaseSetQueue)
    {
        super("PhasedMerger");
        mPrimaryPhaseSetQueue = primaryPhaseSetQueue;
        mPrimaryPhaseSetCount = primaryPhaseSetQueue.size();

        mConfig = config;
        mAssemblyMerger = new AssemblyMerger();

        mPhasedResults = Lists.newArrayList();
    }

    @Override
    public void run()
    {
        while(true)
        {
            try
            {
                int remainingCount = mPrimaryPhaseSetQueue.size();
                int processedCount = mPrimaryPhaseSetCount - remainingCount;

                Set<ExtendedAssembly> primaryPhaseSet = mPrimaryPhaseSetQueue.remove();

                mPerfCounter.start();
                processPhaseSet(primaryPhaseSet);

                mPerfCounter.stop();
                // stopCheckLog(primaryAssembly.toString(), mConfig.PerfLogTime);

                if(processedCount > 0 && (processedCount % TASK_LOG_COUNT) == 0)
                {
                    SV_LOGGER.info("processed {} primary phase sets, remaining({})", processedCount, remainingCount);
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

    public void processPhaseSet(final Set<ExtendedAssembly> primaryPhaseSet)
    {
        Set<ExtendedAssembly> mergedPhaseSets = mAssemblyMerger.primaryPhasedMerging(primaryPhaseSet);

        mPhasedResults.add(mergedPhaseSets);
    }

    public List<Set<ExtendedAssembly>> phasedResults() { return mPhasedResults; }

    public static List<Set<ExtendedAssembly>> mergePhasedResults(final List<PhasedMerger> phasedMergers)
    {
        List<Set<ExtendedAssembly>> combinedResults = Lists.newArrayList();
        phasedMergers.forEach(x -> combinedResults.addAll(x.phasedResults()));
        return combinedResults;
    }
}
