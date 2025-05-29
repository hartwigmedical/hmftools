package com.hartwig.hmftools.esvee.assembly.phase;

import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;

import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.assembly.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.read.BamReader;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.PhaseGroup;
import com.hartwig.hmftools.esvee.assembly.types.ThreadTask;
import com.hartwig.hmftools.common.perf.TaskQueue;

public class PhaseSetTask extends ThreadTask
{
    private final AssemblyConfig mConfig;
    private final TaskQueue mPhaseGroups;

    private final RemoteReadExtractor mRemoteReadExtractor;

    public PhaseSetTask(final AssemblyConfig config, final BamReader bamReader, TaskQueue phaseGroups)
    {
        super("PhaseSets");

        mConfig = config;
        mPhaseGroups = phaseGroups;

        mRemoteReadExtractor = new RemoteReadExtractor(bamReader);
    }

    public static List<PhaseSetTask> createThreadTasks(
            final AssemblyConfig config, final List<PhaseGroup> phaseGroups, final List<BamReader> bamReaders,
            final int taskCount, final List<Thread> threadTasks)
    {
        List<PhaseSetTask> phaseSetTasks = Lists.newArrayList();

        Queue<PhaseGroup> phaseGroupQueue = new ConcurrentLinkedQueue<>();
        phaseGroupQueue.addAll(phaseGroups);

        TaskQueue taskQueue = new TaskQueue(phaseGroupQueue, "phase groups", 10000);

        for(int i = 0; i < taskCount; ++i)
        {
            PhaseSetTask phaseSetTask = new PhaseSetTask(config, bamReaders.get(i), taskQueue);
            phaseSetTasks.add(phaseSetTask);
            threadTasks.add(phaseSetTask);
        }

        if(taskCount > 1)
        {
            SV_LOGGER.debug("splitting {} phase groups across {} threads", phaseGroups.size(), taskCount);
        }

        return phaseSetTasks;
    }

    @Override
    public void run()
    {
        while(true)
        {
            try
            {
                PhaseGroup phaseGroup = (PhaseGroup)mPhaseGroups.removeItem();

                if(mConfig.PhaseProcessingLimit > 0 && phaseGroup.assemblyCount() > mConfig.PhaseProcessingLimit)
                    continue;

                mPerfCounter.start();

                // where there are more than 2 assemblies, start with the ones with the most support and overlapping junction reads
                PhaseSetBuilder phaseSetBuilder = new PhaseSetBuilder(mConfig.RefGenome, mRemoteReadExtractor, phaseGroup);

                try
                {
                    phaseSetBuilder.buildPhaseSets();
                    phaseGroup.finalisePhaseSetAlignments();
                }
                catch(Exception e)
                {
                    SV_LOGGER.error("failed building phase group({}) sets", phaseGroup);

                    for(JunctionAssembly assembly : phaseGroup.assemblies())
                    {
                        SV_LOGGER.info("assembly: {}", assembly);
                    }

                    e.printStackTrace();
                    System.exit(1);
                }

                stopCheckLog(format("phaseGroupId(%d) assemblies(%d)", phaseGroup.id(), phaseGroup.assemblyCount()), mConfig.PerfLogTime);
            }
            catch(NoSuchElementException e)
            {
                SV_LOGGER.trace("all phase tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    public int totalRemoteReadsMatched() { return mRemoteReadExtractor.remoteReadsMatched(); }
}
