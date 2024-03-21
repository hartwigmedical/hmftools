package com.hartwig.hmftools.esvee;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.BAM_READ_JUNCTION_BUFFER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.setAssemblyOutcome;
import static com.hartwig.hmftools.esvee.common.JunctionGroup.buildJunctionGroups;
import static com.hartwig.hmftools.esvee.output.WriteType.ASSEMBLIES;
import static com.hartwig.hmftools.esvee.output.WriteType.ASSEMBLY_BAM;
import static com.hartwig.hmftools.esvee.output.WriteType.READS;

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
import com.hartwig.hmftools.esvee.output.AssemblyReadWriter;
import com.hartwig.hmftools.esvee.output.AssemblyWriter;
import com.hartwig.hmftools.esvee.output.BamWriter;
import com.hartwig.hmftools.esvee.phasing.PhaseGroupBuilder;
import com.hartwig.hmftools.esvee.assembly.JunctionGroupAssembler;
import com.hartwig.hmftools.esvee.phasing.PhaseSetTask;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.JunctionGroup;
import com.hartwig.hmftools.esvee.common.PhaseGroup;
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

    private final Map<String,List<Junction>> mChrJunctionsMap;
    private final Map<String,List<JunctionGroup>> mJunctionGroupMap;

    private final List<BamReader> mBamReaders;

    private final List<PerformanceCounter> mPerfCounters;

    public JunctionProcessor(final SvConfig config)
    {
        mConfig = config;

        mChrJunctionsMap = Maps.newHashMap();
        mJunctionGroupMap = Maps.newHashMap();
        mBamReaders = Lists.newArrayList();

        mResultsWriter = new ResultsWriter(mConfig);

        mPerfCounters = Lists.newArrayList();
    }

    public boolean loadJunctionFiles()
    {
        if(!mConfig.SpecificJunctions.isEmpty())
        {
            for(Junction junction : mConfig.SpecificJunctions)
            {
                List<Junction> chrJunctions = mChrJunctionsMap.get(junction.Chromosome);
                if(chrJunctions == null)
                {
                    chrJunctions = Lists.newArrayList();
                    mChrJunctionsMap.put(junction.Chromosome, chrJunctions);
                }

                chrJunctions.add(junction);
            }

            return true;
        }

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

            List<JunctionAssembly> finalAssemblies = prepareFinalAssemblies();

            writeAssemblyOutput(finalAssemblies);

            writeVariantVcf(finalAssemblies);

            // note this is written after the VCF since writing reassigns the reads to the output BAM, there-by removing their association
            // with the BAM they were read from (ie as used in tumor vs ref counts)
            writeAssemblyBam(finalAssemblies);

            if(mConfig.PerfDebug || (!mConfig.SpecificChrRegions.hasFilters() && mConfig.SpecificJunctions.isEmpty()))
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

        List<JunctionGroupAssembler> primaryAssemblyTasks = JunctionGroupAssembler.createThreadTasks(
                junctionGroups, mBamReaders, mConfig, mResultsWriter, taskCount, threadTasks);

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
        PhaseGroupBuilder phaseGroupBuilder = new PhaseGroupBuilder(mConfig, mJunctionGroupMap, mResultsWriter.phaseGroupBuildWriter());

        phaseGroupBuilder.buildGroups();

        mPerfCounters.add(ThreadTask.mergePerfCounters(phaseGroupBuilder.localBuilderTasks()));
        mPerfCounters.add(ThreadTask.mergePerfCounters(phaseGroupBuilder.remoteBuilderTasks()));

        List<PhaseGroup> phaseGroups = phaseGroupBuilder.phaseGroups();

        if(phaseGroups.isEmpty())
            return;

        SV_LOGGER.info("building phase sets from {} phase groups", phaseGroups.size());

        List<Thread> threadTasks = new ArrayList<>();

        List<PhaseSetTask> phaseSetTasks = PhaseSetTask.createThreadTasks(mConfig, phaseGroupBuilder.phaseGroups(), mConfig.Threads, threadTasks);

        if(!runThreadTasks(threadTasks))
            System.exit(1);

        mPerfCounters.add(ThreadTask.mergePerfCounters(phaseSetTasks.stream().collect(Collectors.toList())));

        SV_LOGGER.info("created {} phase sets", phaseGroups.stream().mapToInt(x -> x.phaseSets().size()).sum());
    }

    private void alignPhaseSets()
    {
        /*
        Alignment alignment = new Alignment(mConfig, new BwaAligner(mConfig));

        for(PhaseGroup phaseGroup : phaseGroupBuilder.phaseGroups())
        {
            for(PhaseSet phaseSet : phaseGroup.phaseSets())
            {
                alignment.processPhaseSet(phaseSet);
            }
        }
        */

    }

    private List<JunctionAssembly> prepareFinalAssemblies()
    {
        List<JunctionAssembly> allAssemblies = Lists.newArrayList();

        int assemblyId = 0;
        for(List<JunctionGroup> junctionGroups : mJunctionGroupMap.values())
        {
            for(JunctionGroup junctionGroup : junctionGroups)
            {
                // set ID first so any references have them set - may need or could to be done earlier once
                // all references between then are established

                // NOTE: is there a cleaner place to do this? eg prior to alignment after phasing is done?
                // junctionGroup.addBranchedAssemblies();

                for(JunctionAssembly assembly : junctionGroup.junctionAssemblies())
                {
                    assembly.setId(assemblyId++);
                    allAssemblies.add(assembly);
                    setAssemblyOutcome(assembly);

                    for(JunctionAssembly branchedAssembly : assembly.branchedAssemblies())
                    {
                        branchedAssembly.setId(assembly.id());
                        allAssemblies.add(branchedAssembly);
                        setAssemblyOutcome(branchedAssembly);
                    }
                }

            }
        }

        // TODO: assembly ID could be set once assemblies have been ordered, and a unique assebly string ID also set in a
        // fairly deterministic way based on final coords, with any duplicates given an extra incrementor

        Collections.sort(allAssemblies, Comparator.comparing(x -> x.junction()));

        return allAssemblies;
    }

    private void writeAssemblyOutput(final List<JunctionAssembly> assemblies)
    {
        if(mConfig.WriteTypes.contains(ASSEMBLIES) || mConfig.WriteTypes.contains(READS))
        {
            SV_LOGGER.debug("writing assembly TSV output");

            AssemblyReadWriter readWriter = mConfig.WriteTypes.contains(READS) ? mResultsWriter.readWriter() : null;
            AssemblyWriter assemblyWriter = mConfig.WriteTypes.contains(ASSEMBLIES) ? mResultsWriter.assemblyWriter() : null;

            for(JunctionAssembly assembly : assemblies)
            {
                if(assemblyWriter != null)
                    assemblyWriter.writeAssembly(assembly);

                if(readWriter != null)
                    readWriter.writeAssemblyReads(assembly);
            }
        }
    }

    private void writeAssemblyBam(final List<JunctionAssembly> assemblies)
    {
        // write BAM records
        if(mConfig.WriteTypes.contains(ASSEMBLY_BAM) && mResultsWriter.bamWriter().isValid())
        {
            SV_LOGGER.debug("writing assembly BAM");

            BamWriter bamWriter = mResultsWriter.bamWriter();

            for(JunctionAssembly assembly : assemblies)
            {
                bamWriter.writeAssembly(assembly);
            }
        }
    }

    private void writeVariantVcf(final List<JunctionAssembly> assemblies)
    {
        if(!mConfig.WriteTypes.contains(WriteType.VCF))
            return;

        SV_LOGGER.debug("writing variant VCF");

        VcfWriter vcfWriter = new VcfWriter(mConfig);

        for(JunctionAssembly assembly : assemblies)
        {
            vcfWriter.addVariant(assembly);
        }

        vcfWriter.close();
    }

    public void close()
    {
        mResultsWriter.close();
    }
}
