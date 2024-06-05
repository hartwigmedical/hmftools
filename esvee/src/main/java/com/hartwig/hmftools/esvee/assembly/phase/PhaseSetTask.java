package com.hartwig.hmftools.esvee.assembly.phase;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.read.BamReader;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.PhaseGroup;
import com.hartwig.hmftools.esvee.assembly.types.PhaseSet;
import com.hartwig.hmftools.esvee.assembly.types.ThreadTask;

public class PhaseSetTask extends ThreadTask
{
    private final AssemblyConfig mConfig;
    private final Queue<PhaseGroup> mPhaseGroups;
    private final int mPhaseGroupCount;

    private final RemoteRegionAssembler mRemoteRegionAssembler;

    public PhaseSetTask(final AssemblyConfig config, final BamReader bamReader, final Queue<PhaseGroup> phaseGroups)
    {
        super("PhaseSets");

        mConfig = config;
        mPhaseGroups = phaseGroups;
        mPhaseGroupCount = mPhaseGroups.size();

        mRemoteRegionAssembler = new RemoteRegionAssembler(config.RefGenome, bamReader);
    }

    public static List<PhaseSetTask> createThreadTasks(
            final AssemblyConfig config, final List<PhaseGroup> phaseGroups, final List<BamReader> bamReaders,
            final int taskCount, final List<Thread> threadTasks)
    {
        List<PhaseSetTask> phaseSetTasks = Lists.newArrayList();

        Queue<PhaseGroup> phaseGroupQueue = new ConcurrentLinkedQueue<>();
        phaseGroupQueue.addAll(phaseGroups);

        for(int i = 0; i < taskCount; ++i)
        {
            PhaseSetTask phaseSetTask = new PhaseSetTask(config, bamReaders.get(i), phaseGroupQueue);
            phaseSetTasks.add(phaseSetTask);
            threadTasks.add(phaseSetTask);
        }

        if(taskCount > 1)
        {
            SV_LOGGER.debug("splitting {} phase groups across {} threads", phaseGroups.size(), taskCount);
        }

        return phaseSetTasks;
    }

    private static final int PHASE_GROUP_LOG_COUNT = 10000;

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

                if(mConfig.PhaseProcessingLimit > 0 && phaseGroup.assemblyCount() > mConfig.PhaseProcessingLimit)
                    continue;

                mPerfCounter.start();

                // where there are more than 2 assemblies, start with the ones with the most support and overlapping junction reads
                PhaseSetBuilder phaseSetBuilder = new PhaseSetBuilder(mConfig.RefGenome, mRemoteRegionAssembler, phaseGroup);

                try
                {
                    phaseSetBuilder.buildPhaseSets();
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

                // now that all support has been added, set the supporting read indices for each assembly
                phaseGroup.assemblies().forEach(x -> x.setReadIndices());

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

    public int totalRemoteReadsMatched() { return mRemoteRegionAssembler.totalRemoteReadsMatched(); }
}
