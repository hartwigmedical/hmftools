package com.hartwig.hmftools.esvee.assembly.phase;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.types.ThreadTask.mergePerfCounters;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.JunctionGroup;
import com.hartwig.hmftools.esvee.assembly.types.PhaseGroup;
import com.hartwig.hmftools.esvee.assembly.types.ThreadTask;
import com.hartwig.hmftools.esvee.assembly.output.PhaseGroupBuildWriter;

public class PhaseGroupBuilder
{
    private final AssemblyConfig mConfig;
    private final Map<String, List<JunctionGroup>> mJunctionGroupMap;
    private final List<PhaseGroup> mPhaseGroups;

    private final PhaseGroupBuildWriter mWriter;

    public PhaseGroupBuilder(
            final AssemblyConfig config, final Map<String,List<JunctionGroup>> junctionGroupMap, final PhaseGroupBuildWriter writer)
    {
        mConfig = config;
        mWriter = writer;
        mJunctionGroupMap = junctionGroupMap;
        mPhaseGroups = Lists.newArrayList();

        // setting an index for each junction group (within its chromosome) allows easy access to an assembly's own group during look-ups
        for(List<JunctionGroup> junctionGroups : mJunctionGroupMap.values())
        {
            for(int index = 0; index < junctionGroups.size(); ++index)
            {
                junctionGroups.get(index).setIndex(index);
            }
        }
    }

    public List<PhaseGroup> phaseGroups() { return mPhaseGroups; }

    public void buildGroups(final List<PerformanceCounter> perfCounters)
    {
        List<JunctionGroup> allJunctionGroups = Lists.newArrayList();

        for(List<JunctionGroup> junctionGroups : mJunctionGroupMap.values())
        {
            junctionGroups.stream().filter(x -> !x.junctionAssemblies().isEmpty()).forEach(x -> allJunctionGroups.add(x));
        }

        int junctionGroupCount = allJunctionGroups.size();

        int taskCount = min(mConfig.Threads, junctionGroupCount);

        // stage one is local phase group building
        Queue<JunctionGroup> junctionGroupQueue = new ConcurrentLinkedQueue<>();

        junctionGroupQueue.addAll(allJunctionGroups);

        List<Thread> threadTasks = new ArrayList<>();

        List<LocalGroupBuilder> localBuilderTasks = Lists.newArrayList();

        for(int i = 0; i < taskCount; ++i)
        {
            LocalGroupBuilder groupBuilderTask = new LocalGroupBuilder(mConfig, junctionGroupQueue, mWriter);
            localBuilderTasks.add(groupBuilderTask);
            threadTasks.add(groupBuilderTask);
        }

        if(taskCount > 1)
        {
            SV_LOGGER.debug("phase-grouping {} junction groups across {} threads", junctionGroupCount, taskCount);
        }

        if(!runThreadTasks(threadTasks))
            System.exit(1);

        mergePerfCounters(perfCounters, localBuilderTasks.stream().collect(Collectors.toList()));

        localBuilderTasks.forEach(x -> mPhaseGroups.addAll(x.phaseGroups()));

        junctionGroupQueue.addAll(allJunctionGroups);

        SV_LOGGER.info("building remote phase groups, current group count({})", mPhaseGroups.size());

        threadTasks = new ArrayList<>();

        List<RemoteGroupBuilder> remoteBuilderTasks = Lists.newArrayList();

        for(int i = 0; i < taskCount; ++i)
        {
            RemoteGroupBuilder groupBuilderTask = new RemoteGroupBuilder(mConfig, junctionGroupQueue, mJunctionGroupMap, mWriter);
            remoteBuilderTasks.add(groupBuilderTask);
            threadTasks.add(groupBuilderTask);
        }

        if(!runThreadTasks(threadTasks))
            System.exit(1);

        remoteBuilderTasks.forEach(x -> mPhaseGroups.addAll(x.phaseGroups()));

        // clean-up phase groups which were transferred into another group
        remoteBuilderTasks.forEach(x -> x.removedPhaseGroups().forEach(y -> mPhaseGroups.remove(y)));

        mergePerfCounters(perfCounters, remoteBuilderTasks.stream().collect(Collectors.toList()));

        remoteBuilderTasks.forEach(x -> x.logStats());

        // add unphased assemblies to their own phase group
        for(JunctionGroup junctionGroup : allJunctionGroups)
        {
            for(JunctionAssembly assembly : junctionGroup.junctionAssemblies())
            {
                if(assembly.phaseGroup() == null)
                {
                    mPhaseGroups.add(new PhaseGroup(assembly, null));
                }
            }
        }

        for(int i = 0; i < mPhaseGroups.size(); ++i)
        {
            mPhaseGroups.get(i).setId(i);
        }

        // run validation
        if(mConfig.PerfDebug)
        {
            // check if an assembly is in 2 phase groups
            List<JunctionAssembly> assemblies = Lists.newArrayList();

            for(List<JunctionGroup> junctionGroups : mJunctionGroupMap.values())
            {
                junctionGroups.forEach(x -> assemblies.addAll(x.junctionAssemblies()));
            }

            for(JunctionAssembly assembly : assemblies)
            {
                List<PhaseGroup> phaseGroups = mPhaseGroups.stream().filter(x -> x.assemblies().contains(assembly)).collect(Collectors.toList());

                if(phaseGroups.size() > 1)
                {
                    SV_LOGGER.error("asm({}) in {} phase groups", assembly, phaseGroups.size());
                }
            }
        }
    }

    protected static void linkToPhaseGroups(
            final PhaseGroup phaseGroup, final JunctionAssembly assembly, final JunctionAssembly otherAssembly,
            final Set<PhaseGroup> phaseGroups, @Nullable final List<PhaseGroup> removedPhaseGroups,
            final PhaseGroupBuildWriter writer, @Nullable final String linkType)
    {
        // not thread-safe but called with synchronisation when remote assemblies are linked

        // for phased assemblies, the phase group scenarios are:
        // - assemblies are already in the same phase group
        // - no phase group exists for either so make a new one
        // - the other assembly already has an existing phase group and this one doesn't so just add it
        // - the other assembly already has an existing phase group, so does this one and so transfer this to the other
        if(phaseGroup == null)
        {
            if(otherAssembly.phaseGroup() != null)
            {
                if(writer.enabled())
                    writer.writePhaseGroupBuild(assembly, otherAssembly, linkType, otherAssembly.phaseGroup().assemblyCount());

                otherAssembly.phaseGroup().addAssembly(assembly);
            }
            else
            {
                phaseGroups.add(new PhaseGroup(assembly, otherAssembly));

                if(writer.enabled())
                {
                    writer.writePhaseGroupBuild(assembly, otherAssembly, linkType, 0);
                    writer.writePhaseGroupBuild(otherAssembly, assembly, linkType, 1);
                }
            }
        }
        else
        {
            if(otherAssembly.phaseGroup() != null)
            {
                // transfer from the smaller group to the larger to the other one
                PhaseGroup destPhaseGroup = phaseGroup.assemblyCount() >= otherAssembly.phaseGroup().assemblyCount() ?
                        phaseGroup : otherAssembly.phaseGroup();

                PhaseGroup srcPhaseGroup = destPhaseGroup == phaseGroup ? otherAssembly.phaseGroup() : phaseGroup;

                if(writer.enabled())
                    writer.writePhaseGroupBuild(assembly, otherAssembly, linkType, destPhaseGroup.assemblyCount());

                destPhaseGroup.transferAssemblies(srcPhaseGroup);

                if(removedPhaseGroups != null)
                    removedPhaseGroups.add(srcPhaseGroup);
                else
                    phaseGroups.remove(srcPhaseGroup);
            }
            else
            {
                if(writer.enabled())
                    writer.writePhaseGroupBuild(otherAssembly, assembly, linkType, phaseGroup.assemblyCount());

                phaseGroup.addAssembly(otherAssembly);
            }
        }
    }
}
