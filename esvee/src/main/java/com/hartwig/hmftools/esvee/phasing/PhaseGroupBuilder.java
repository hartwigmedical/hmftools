package com.hartwig.hmftools.esvee.phasing;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.JunctionGroup;
import com.hartwig.hmftools.esvee.common.PhaseGroup;
import com.hartwig.hmftools.esvee.common.ThreadTask;
import com.hartwig.hmftools.esvee.output.PhaseGroupBuildWriter;

public class PhaseGroupBuilder
{
    private final SvConfig mConfig;
    private final Map<String, List<JunctionGroup>> mJunctionGroupMap;
    private final List<PhaseGroup> mPhaseGroups;

    private final List<ThreadTask> mLocalBuilderTasks;
    private final List<ThreadTask> mRemoteBuilderTasks;

    private final PhaseGroupBuildWriter mWriter;

    public PhaseGroupBuilder(
            final SvConfig config, final Map<String,List<JunctionGroup>> junctionGroupMap, final PhaseGroupBuildWriter writer)
    {
        mConfig = config;
        mWriter = writer;
        mJunctionGroupMap = junctionGroupMap;
        mPhaseGroups = Lists.newArrayList();
        mLocalBuilderTasks = Lists.newArrayList();
        mRemoteBuilderTasks = Lists.newArrayList();

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
    public List<ThreadTask> localBuilderTasks() { return mLocalBuilderTasks; }
    public List<ThreadTask> remoteBuilderTasks() { return mRemoteBuilderTasks; }

    public void buildGroups()
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


        localBuilderTasks.forEach(x -> mPhaseGroups.addAll(x.phaseGroups()));
        localBuilderTasks.forEach(x -> mLocalBuilderTasks.add(x));

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

        remoteBuilderTasks.forEach(x -> mRemoteBuilderTasks.add(x));

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

    /* FIXME: needs to be made threadsafe

        if(!mConfig.SkipDiscordant)
        {
            for(JunctionGroup otherJunctionGroup : linkedJunctionGroups)
            {
                // the linked junction group may have discordant read groups which are required by this assembly
                for(DiscordantGroup discordantGroup : otherJunctionGroup.discordantGroups())
                {
                    // ignore any discordant group which actually shares the same reads
                    if(discordantGroup.reads().stream().anyMatch(x -> assembly.hasReadSupport(x)))
                        continue;

                    // as above, check that remote regions overlap before checking reads
                    boolean matched = false;

                    for(RemoteRegion remoteRegion : assembly.remoteRegions())
                    {
                        if(remoteRegion.isSuppOnlyRegion())
                            continue;

                        if(!remoteRegion.overlaps(
                                discordantGroup.chromosome(), discordantGroup.minAlignedPosition(), discordantGroup.maxAlignedPosition(),
                                discordantGroup.orientation()))
                        {
                            continue;
                        }

                        if(remoteRegion.readIds().stream().anyMatch(x -> discordantGroup.hasFragment(x)))
                        {
                            matched = true;
                            break;
                        }
                    }

                    if(!matched)
                        continue;

                    if(phaseGroup == null)
                        phaseGroup = new PhaseGroup(assembly, null);

                    phaseGroup.addDiscordantGroup(discordantGroup);
                }
            }
        }
    */
}
