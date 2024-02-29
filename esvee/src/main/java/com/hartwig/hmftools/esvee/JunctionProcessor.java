package com.hartwig.hmftools.esvee;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.BAM_READ_JUNCTION_BUFFER;
import static com.hartwig.hmftools.esvee.common.JunctionGroup.buildJunctionGroups;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.esvee.alignment.Alignment;
import com.hartwig.hmftools.esvee.alignment.BwaAligner;
import com.hartwig.hmftools.esvee.assembly.PhaseGroupBuilder;
import com.hartwig.hmftools.esvee.assembly.JunctionGroupAssembler;
import com.hartwig.hmftools.esvee.assembly.PhaseSetTask;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.JunctionGroup;
import com.hartwig.hmftools.esvee.common.ThreadTask;
import com.hartwig.hmftools.esvee.output.ResultsWriter;
import com.hartwig.hmftools.esvee.output.WriteType;
import com.hartwig.hmftools.esvee.read.BamReader;
import com.hartwig.hmftools.esvee.output.VcfWriter;
import com.hartwig.hmftools.esvee.read.ReadStats;

public class JunctionProcessor
{
    private final SvConfig mConfig;
    private final ResultsWriter mResultsWriter;
    private final VcfWriter mVcfWriter;

    private final Map<String,List<Junction>> mChrJunctionsMap;
    private final Map<String,List<JunctionGroup>> mJunctionGroupMap;

    private final Alignment mAlignment;

    private final List<BamReader> mBamReaders;

    private final List<PerformanceCounter> mPerfCounters;

    public JunctionProcessor(final SvConfig config)
    {
        mConfig = config;

        mChrJunctionsMap = Maps.newHashMap();
        mJunctionGroupMap = Maps.newHashMap();
        mBamReaders = Lists.newArrayList();

        mAlignment = new Alignment(mConfig, new BwaAligner(mConfig));

        mResultsWriter = new ResultsWriter(mConfig);
        mVcfWriter = new VcfWriter(mConfig);

        mPerfCounters = Lists.newArrayList();
    }

    public boolean loadJunctionFiles()
    {
        for(String junctionFile : mConfig.JunctionFiles)
        {
            Map<String,List<Junction>> newJunctionsMap = Junction.loadJunctions(
                    junctionFile, mConfig.SpecificChrRegions, mConfig.SkipDiscordant);

            if(newJunctionsMap == null)
                return false;

            Junction.mergeJunctions(mChrJunctionsMap, newJunctionsMap);
        }

        // if(mConfig.PerfDebug && !validateJunctionMap(mChrJunctionsMap))
        //    System.exit(1);

        if(mConfig.JunctionFiles.size() > 1)
        {
            SV_LOGGER.debug("merged into {} junctions", mChrJunctionsMap.values().stream().mapToInt(x -> x.size()).sum());
        }

        return !mChrJunctionsMap.isEmpty();
    }

    public void run()
    {
        try
        {
            // create junction groups from existing junctions
            for(Map.Entry<String,List<Junction>> entry : mChrJunctionsMap.entrySet())
            {
                mJunctionGroupMap.put(entry.getKey(), buildJunctionGroups(entry.getValue(), BAM_READ_JUNCTION_BUFFER));
            }

            int junctionCount = mChrJunctionsMap.values().stream().mapToInt(x -> x.size()).sum();

            int taskCount = min(junctionCount, mConfig.Threads);

            loadBamFiles(taskCount);

            runPrimaryAssembly();

            formPhaseGroups();

            alignPhaseSets();

            // write here to show PPG data - may move back later on
            if(mConfig.WriteTypes.contains(WriteType.ASSEMBLIES))
            {
                SV_LOGGER.debug("writing assembly data");

                int assemblyId = 0;
                for(List<JunctionGroup> junctionGroups : mJunctionGroupMap.values())
                {
                    for(JunctionGroup junctionGroup : junctionGroups)
                    {
                        // set ID first so any references have them set - may need or could to be done earlier once
                        // all references between then are established

                        // NOTE: is there a cleaner place to do this? eg prior to alignment after phasing is done?
                        junctionGroup.addBranchedAssemblies();

                        for(JunctionAssembly assembly : junctionGroup.junctionAssemblies())
                        {
                            assembly.setId(assemblyId++);
                        }

                        junctionGroup.junctionAssemblies().forEach(x -> mResultsWriter.writeAssembly(x));
                    }
                }
            }

            if(mConfig.PerfDebug || !mConfig.SpecificChrRegions.hasFilters())
            {
                mPerfCounters.forEach(x -> x.logStats());
            }
        }
        catch(Exception e)
        {
            SV_LOGGER.error("process run error: {}", e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }

    private void loadBamFiles(int taskCount)
    {
        for(int i = 0; i < taskCount; ++i)
        {
            BamReader bamReader = new BamReader(mConfig);
            mBamReaders.add(bamReader);
        }
    }

    private void runPrimaryAssembly()
    {
        int taskCount = mBamReaders.size();

        List<JunctionGroup> junctionGroups = Lists.newArrayList();
        mJunctionGroupMap.values().forEach(x -> junctionGroups.addAll(x));

        // sorted by groups with the latest range, assuming that these will have the most reads and junctions
        Collections.sort(junctionGroups, Comparator.comparingInt(x -> -x.range()));

        List<Thread> threadTasks = new ArrayList<>();

        // Primary Junction Assembly
        List<JunctionGroupAssembler> primaryAssemblyTasks = JunctionGroupAssembler.createThreadTasks(
                junctionGroups, mBamReaders, mConfig, taskCount, threadTasks);

        if(!runThreadTasks(threadTasks))
            System.exit(1);

        int totalJunctionAssemblies = junctionGroups.stream().mapToInt(x -> x.junctionAssemblies().size()).sum();
        int totalDiscordantGroups = junctionGroups.stream().mapToInt(x -> x.discordantGroups().size()).sum();

        ReadStats combinedReadStats = new ReadStats();
        primaryAssemblyTasks.forEach(x -> combinedReadStats.merge(x.readStats()));

        SV_LOGGER.info("created {} junction assemblies, {} discordant groups", totalJunctionAssemblies, totalDiscordantGroups);

        mPerfCounters.add(ThreadTask.mergePerfCounters(primaryAssemblyTasks.stream().collect(Collectors.toList())));

        int totalCachedReads = junctionGroups.stream().mapToInt(x -> x.candidateReadCount()).sum();

        SV_LOGGER.info("cached read count({}) from {} junction groups, stats: {}",
                totalCachedReads, junctionGroups.size(), combinedReadStats);
    }

    private void formPhaseGroups()
    {
        PhaseGroupBuilder phaseGroupBuilder = new PhaseGroupBuilder(mConfig, mJunctionGroupMap);

        phaseGroupBuilder.buildGroups();

        SV_LOGGER.info("building phase sets from {} phase groups", phaseGroupBuilder.phaseGroups().size());

        List<Thread> threadTasks = new ArrayList<>();

        List<PhaseSetTask> phaseSetTasks = PhaseSetTask.createThreadTasks(mConfig, phaseGroupBuilder.phaseGroups(), mConfig.Threads, threadTasks);

        if(!runThreadTasks(threadTasks))
            System.exit(1);

        mPerfCounters.add(ThreadTask.mergePerfCounters(phaseSetTasks.stream().collect(Collectors.toList())));

        SV_LOGGER.info("created {} phase sets", phaseGroupBuilder.phaseGroups().stream().mapToInt(x -> x.phaseSets().size()).sum());
    }

    private void alignPhaseSets()
    {
        /*
        for(PhaseGroup phaseGroup : phaseGroupBuilder.phaseGroups())
        {
            for(PhaseSet phaseSet : phaseGroup.phaseSets())
            {
                mAlignment.processPhaseSet(phaseSet);
            }
        }
        */

    }

    public void close()
    {
        mVcfWriter.close();
        mResultsWriter.close();
    }
}
