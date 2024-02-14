package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.common.PhaseGroup;
import com.hartwig.hmftools.esvee.common.PhaseSet;
import com.hartwig.hmftools.esvee.common.ThreadTask;

public class PhaseSetTask extends ThreadTask
{
    private final SvConfig mConfig;
    private final Queue<PhaseGroup> mPhaseGroups;
    private final int mPhaseGroupCount;

    public PhaseSetTask(final SvConfig config, final Queue<PhaseGroup> phaseGroups)
    {
        super("PhaseSets");

        mConfig = config;
        mPhaseGroups = phaseGroups;
        mPhaseGroupCount = mPhaseGroups.size();
    }

    public static List<PhaseSetTask> createThreadTasks(
            final SvConfig config, final List<PhaseGroup> phaseGroups, final int taskCount, final List<Thread> threadTasks)
    {
        List<PhaseSetTask> phaseSetTasks = Lists.newArrayList();

        Queue<PhaseGroup> phaseGroupQueue = new ConcurrentLinkedQueue<>();
        phaseGroupQueue.addAll(phaseGroups);

        for(int i = 0; i < taskCount; ++i)
        {
            PhaseSetTask phaseSetTask = new PhaseSetTask(config, phaseGroupQueue);
            phaseSetTasks.add(phaseSetTask);
            threadTasks.add(phaseSetTask);
        }

        SV_LOGGER.debug("splitting {} phase groups across {} threads", phaseGroups.size(), taskCount);

        return phaseSetTasks;
    }

    private static final int PHASE_GROUP_LOG_COUNT = 1000;

    @Override
    public void run()
    {
        while(true)
        {
            try
            {
                int remainingCount = mPhaseGroups.size();
                int processedCount = mPhaseGroupCount - remainingCount;

                PhaseGroup phaseGroup = mPhaseGroups.remove();

                mPerfCounter.start();

                // where there are more than 2 assemblies, start with the ones with the most support and overlapping junction reads
                PhaseSetBuilder phaseSetBuilder = new PhaseSetBuilder(phaseGroup);

                phaseSetBuilder.buildPhaseSets();

                // also set phase set IDs
                int phaseSetId = 0;
                for(PhaseSet phaseSet : phaseGroup.phaseSets())
                {
                    phaseSet.setId(phaseSetId++);
                }

                stopCheckLog(phaseGroup.toString(), mConfig.PerfLogTime);

                if(processedCount > 0 && (processedCount % PHASE_GROUP_LOG_COUNT) == 0)
                {
                    SV_LOGGER.debug("processed {} phase groups, remaining({})", processedCount, remainingCount);
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

}
