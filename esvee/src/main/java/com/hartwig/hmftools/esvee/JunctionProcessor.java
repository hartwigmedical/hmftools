package com.hartwig.hmftools.esvee;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.BAM_READ_JUNCTION_BUFFER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.DISCORDANT_FRAGMENT_LENGTH;
import static com.hartwig.hmftools.esvee.alignment.Alignment.skipJunctionAssembly;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.setAssemblyOutcome;
import static com.hartwig.hmftools.esvee.assembly.types.ThreadTask.mergePerfCounters;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.PREP_FRAG_LENGTH_FILE_ID;
import static com.hartwig.hmftools.esvee.assembly.types.JunctionGroup.buildJunctionGroups;
import static com.hartwig.hmftools.esvee.assembly.output.WriteType.JUNC_ASSEMBLY;
import static com.hartwig.hmftools.esvee.assembly.output.WriteType.ASSEMBLY_BAM;
import static com.hartwig.hmftools.esvee.assembly.output.WriteType.ASSEMBLY_READ;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.esvee.alignment.Alignment;
import com.hartwig.hmftools.esvee.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.alignment.Breakend;
import com.hartwig.hmftools.esvee.alignment.BwaAligner;
import com.hartwig.hmftools.esvee.assembly.output.BreakendWriter;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.PhaseSet;
import com.hartwig.hmftools.esvee.common.FragmentLengthBounds;
import com.hartwig.hmftools.esvee.assembly.output.AssemblyReadWriter;
import com.hartwig.hmftools.esvee.assembly.output.AssemblyWriter;
import com.hartwig.hmftools.esvee.assembly.output.BamWriter;
import com.hartwig.hmftools.esvee.assembly.phase.PhaseGroupBuilder;
import com.hartwig.hmftools.esvee.assembly.JunctionGroupAssembler;
import com.hartwig.hmftools.esvee.assembly.phase.PhaseSetTask;
import com.hartwig.hmftools.esvee.prep.FragmentSizeDistribution;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.JunctionGroup;
import com.hartwig.hmftools.esvee.assembly.types.PhaseGroup;
import com.hartwig.hmftools.esvee.assembly.types.ThreadTask;
import com.hartwig.hmftools.esvee.assembly.output.ResultsWriter;
import com.hartwig.hmftools.esvee.assembly.output.WriteType;
import com.hartwig.hmftools.esvee.assembly.read.BamReader;
import com.hartwig.hmftools.esvee.assembly.output.VcfWriter;
import com.hartwig.hmftools.esvee.assembly.read.ReadStats;

public class JunctionProcessor
{
    private final AssemblyConfig mConfig;
    private final ResultsWriter mResultsWriter;

    private final Map<String,List<Junction>> mChrJunctionsMap;
    private final Map<String,List<JunctionGroup>> mJunctionGroupMap;

    private final List<BamReader> mBamReaders;

    private final List<PerformanceCounter> mPerfCounters;

    public JunctionProcessor(final AssemblyConfig config)
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
                    junctionFile, mConfig.SpecificChrRegions, mConfig.ProcessDiscordant);

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

    private void loadFragmentLengthBounds()
    {
        if(mConfig.JunctionFiles.isEmpty())
            return;

        String fragLengthFilename = mConfig.formPrepInputFilename(mConfig.sampleId(), PREP_FRAG_LENGTH_FILE_ID);

        FragmentLengthBounds fragmentLengthBounds = FragmentSizeDistribution.loadFragmentLengthBounds(fragLengthFilename);

        if(fragmentLengthBounds.isValid())
        {
            SV_LOGGER.info("fragment length bounds(min={} max={})",
                    fragmentLengthBounds.LowerBound, fragmentLengthBounds.UpperBound);

            DISCORDANT_FRAGMENT_LENGTH = fragmentLengthBounds.UpperBound;
        }
    }

    public void run()
    {
        loadFragmentLengthBounds();

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

            List<JunctionAssembly> finalAssemblies = Lists.newArrayList();
            List<AssemblyAlignment> assemblyAlignments = Lists.newArrayList();

            gatherAssemblies(finalAssemblies, assemblyAlignments);

            runAlignment(assemblyAlignments);

            writeAssemblyOutput(finalAssemblies);

            writeVariants(assemblyAlignments);

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

        int assemblyCount = 0;
        int junctionReadCount = 0;
        int candidateReadCount = 0;

        for(JunctionGroup junctionGroup : junctionGroups)
        {
            for(JunctionAssembly assembly : junctionGroup.junctionAssemblies())
            {
                ++assemblyCount;
                junctionReadCount += assembly.supportCount();
                candidateReadCount += assembly.candidateSupport().size();
            }
        }

        ReadStats combinedReadStats = new ReadStats();
        primaryAssemblyTasks.forEach(x -> combinedReadStats.merge(x.readStats()));

        SV_LOGGER.info("created {} junction assemblies reads(junc={} candidate={})",
                assemblyCount, junctionReadCount, candidateReadCount);

        SV_LOGGER.info("extracted read stats: {}", combinedReadStats);

        addPerfCounters(primaryAssemblyTasks.stream().collect(Collectors.toList()));
    }

    private void formPhaseGroups()
    {
        PhaseGroupBuilder phaseGroupBuilder = new PhaseGroupBuilder(mConfig, mJunctionGroupMap, mResultsWriter.phaseGroupBuildWriter());

        phaseGroupBuilder.buildGroups(mPerfCounters);

        List<PhaseGroup> phaseGroups = phaseGroupBuilder.phaseGroups();

        if(phaseGroups.isEmpty())
            return;

        SV_LOGGER.info("building phase sets from {} phase groups", phaseGroups.size());

        List<Thread> threadTasks = new ArrayList<>();

        List<PhaseSetTask> phaseSetTasks = PhaseSetTask.createThreadTasks(
                mConfig, phaseGroupBuilder.phaseGroups(), mBamReaders, mConfig.Threads, threadTasks);

        if(!runThreadTasks(threadTasks))
            System.exit(1);

        addPerfCounters(phaseSetTasks.stream().collect(Collectors.toList()));

        int totalRemoteReadMatched = phaseSetTasks.stream().mapToInt(x -> x.totalRemoteReadsMatched()).sum();

        SV_LOGGER.info("created {} phase sets, remote reads extracted({})",
                phaseGroups.stream().mapToInt(x -> x.phaseSets().size()).sum(), totalRemoteReadMatched);
    }

    private void addPerfCounters(final List<ThreadTask> tasks)
    {
        mergePerfCounters(mPerfCounters, tasks);
    }

    private void runAlignment(final List<AssemblyAlignment> assemblyAlignments)
    {
        if(!mConfig.RunAlignment)
            return;

        boolean useCache = mConfig.AlignmentFile != null;
        Alignment alignment = new Alignment(mConfig, !useCache ? new BwaAligner(mConfig) : null);
        alignment.run(assemblyAlignments, mPerfCounters);
        alignment.close();
    }

    private void gatherAssemblies(final List<JunctionAssembly> allAssemblies, final List<AssemblyAlignment> assemblyAlignments)
    {
        int assemblyId = 0;

        Set<PhaseGroup> phaseGroups = Sets.newHashSet();

        for(List<JunctionGroup> junctionGroups : mJunctionGroupMap.values())
        {
            for(JunctionGroup junctionGroup : junctionGroups)
            {
                for(JunctionAssembly assembly : junctionGroup.junctionAssemblies())
                {
                    setAssemblyOutcome(assembly);

                    if(assembly.phaseGroup() != null)
                    {
                        phaseGroups.add(assembly.phaseGroup());
                    }

                    allAssemblies.add(assembly);
                }
            }
        }

        int assemblyAlignmentId = 0;
        for(PhaseGroup phaseGroup : phaseGroups)
        {
            allAssemblies.addAll(phaseGroup.derivedAssemblies());

            if(mConfig.RunAlignment)
            {
                for(PhaseSet phaseSet : phaseGroup.phaseSets())
                {
                    for(AssemblyLink assemblyLink : phaseSet.assemblyLinks())
                    {
                        assemblyAlignments.add(new AssemblyAlignment(assemblyAlignmentId++, assemblyLink));
                    }
                }

                for(JunctionAssembly assembly : phaseGroup.assemblies())
                {
                    if(assembly.phaseSet() == null && !skipJunctionAssembly(assembly))
                        assemblyAlignments.add(new AssemblyAlignment(assemblyAlignmentId++, assembly));
                }
            }
        }

        Collections.sort(allAssemblies, Comparator.comparing(x -> x.junction()));

        for(JunctionAssembly assembly : allAssemblies)
        {
            assembly.setId(assemblyId++);
        }
    }

    private void writeAssemblyOutput(final List<JunctionAssembly> assemblies)
    {
        if(mConfig.WriteTypes.contains(JUNC_ASSEMBLY) || mConfig.WriteTypes.contains(ASSEMBLY_READ))
        {
            SV_LOGGER.debug("writing assembly TSV output");

            AssemblyReadWriter readWriter = mConfig.WriteTypes.contains(ASSEMBLY_READ) ? mResultsWriter.readWriter() : null;
            AssemblyWriter assemblyWriter = mConfig.WriteTypes.contains(JUNC_ASSEMBLY) ? mResultsWriter.assemblyWriter() : null;

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

    private void writeVariants(final List<AssemblyAlignment> assemblyAlignments)
    {
        if(mConfig.WriteTypes.contains(WriteType.VCF))
        {
            List<Breakend> breakends = Lists.newArrayList();
            assemblyAlignments.forEach(x -> breakends.addAll(x.breakends()));
            Collections.sort(breakends);

            SV_LOGGER.debug("writing variant VCF with {} breakends");

            VcfWriter vcfWriter = new VcfWriter(mConfig);

            for(Breakend breakend : breakends)
            {
                vcfWriter.addBreakend(breakend);
            }

            vcfWriter.close();
        }

        if(mConfig.WriteTypes.contains(WriteType.BREAKEND))
        {
            BreakendWriter breakendWriter = mResultsWriter.breakendWriter();
            assemblyAlignments.forEach(x -> breakendWriter.writeBreakends((x)));
        }

    }

    public void close()
    {
        mResultsWriter.close();
    }
}
