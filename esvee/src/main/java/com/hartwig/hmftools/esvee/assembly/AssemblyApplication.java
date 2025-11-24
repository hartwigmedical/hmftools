package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.floor;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.perf.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.BAM_READ_JUNCTION_BUFFER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.DISC_RATE_DISC_ONLY_INCREMENT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.DISC_RATE_JUNC_INCREMENT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.MAX_OBSERVED_CONCORDANT_FRAG_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.alignment.Alignment.skipUnlinkedJunctionAssembly;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.setAssemblyOutcome;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.DECOY;
import static com.hartwig.hmftools.esvee.assembly.types.ThreadTask.mergePerfCounters;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;
import static com.hartwig.hmftools.esvee.assembly.types.JunctionGroup.buildJunctionGroups;
import static com.hartwig.hmftools.esvee.common.FileCommon.formDiscordantStatsFilename;
import static com.hartwig.hmftools.esvee.common.FileCommon.formFragmentLengthDistFilename;
import static com.hartwig.hmftools.esvee.common.WriteType.ASSEMBLY_BAM;
import static com.hartwig.hmftools.esvee.common.WriteType.ASSEMBLY_READ;
import static com.hartwig.hmftools.esvee.common.WriteType.JUNC_ASSEMBLY;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DISCORDANT_GROUP_MIN_FRAGMENTS_SHORT;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_HOTSPOT_JUNCTION_SUPPORT;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_JUNCTION_SUPPORT;
import static com.hartwig.hmftools.esvee.prep.types.DiscordantStats.loadDiscordantStats;

import java.nio.file.Files;
import java.nio.file.Paths;
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
import com.hartwig.hmftools.common.perf.PerformanceCounter;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.assembly.alignment.Alignment;
import com.hartwig.hmftools.esvee.assembly.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.assembly.alignment.Breakend;
import com.hartwig.hmftools.esvee.assembly.alignment.BwaAligner;
import com.hartwig.hmftools.esvee.assembly.output.BreakendWriter;
import com.hartwig.hmftools.esvee.assembly.types.PhaseSet;
import com.hartwig.hmftools.esvee.common.FragmentLengthBounds;
import com.hartwig.hmftools.esvee.assembly.output.AssemblyReadWriter;
import com.hartwig.hmftools.esvee.assembly.output.AssemblyWriter;
import com.hartwig.hmftools.esvee.assembly.output.BamWriter;
import com.hartwig.hmftools.esvee.assembly.phase.PhaseGroupBuilder;
import com.hartwig.hmftools.esvee.assembly.phase.PhaseSetTask;
import com.hartwig.hmftools.esvee.common.WriteType;
import com.hartwig.hmftools.esvee.prep.FragmentSizeDistribution;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.JunctionGroup;
import com.hartwig.hmftools.esvee.assembly.types.PhaseGroup;
import com.hartwig.hmftools.esvee.assembly.types.ThreadTask;
import com.hartwig.hmftools.esvee.assembly.output.ResultsWriter;
import com.hartwig.hmftools.esvee.assembly.read.BamReader;
import com.hartwig.hmftools.esvee.assembly.output.VcfWriter;
import com.hartwig.hmftools.esvee.assembly.read.ReadStats;
import com.hartwig.hmftools.esvee.prep.types.DiscordantStats;

public class AssemblyApplication
{
    private final AssemblyConfig mConfig;
    private final ResultsWriter mResultsWriter;

    private final Map<String,List<Junction>> mChrJunctionsMap;
    private final Map<String,List<JunctionGroup>> mJunctionGroupMap;

    private final List<BamReader> mBamReaders;

    private final List<PerformanceCounter> mPerfCounters;

    public AssemblyApplication(final ConfigBuilder configBuilder)
    {
        this(configBuilder, false);
    }

    public AssemblyApplication(final ConfigBuilder configBuilder, boolean asSubRoutine)
    {
        mConfig = new AssemblyConfig(configBuilder, asSubRoutine);

        mChrJunctionsMap = Maps.newHashMap();
        mJunctionGroupMap = Maps.newHashMap();
        mBamReaders = Lists.newArrayList();

        mResultsWriter = new ResultsWriter(mConfig);

        mPerfCounters = Lists.newArrayList();
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        SV_LOGGER.info("writing to output directory({}){}",
                mConfig.OutputDir, mConfig.OutputId != null ? format(" outputId(%s)", mConfig.OutputId) : "");

        loadFragmentLengthBounds();

        if(!loadJunctionFiles())
        {
            SV_LOGGER.error("failed to load junction file");
            System.exit(1);
        }

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

            List<Breakend> breakends = Lists.newArrayList();
            assemblyAlignments.forEach(x -> breakends.addAll(x.breakends()));
            Collections.sort(breakends);

            for(int i = 0; i < breakends.size(); ++i)
            {
                breakends.get(i).setId(i);
            }

            writeAssemblyOutput(finalAssemblies);

            writeVariants(assemblyAlignments, breakends);

            // note this is written after the VCF since writing reassigns the reads to the output BAM, there-by removing their association
            // with the BAM they were read from (ie as used in tumor vs ref counts)
            writeAssemblyBam(finalAssemblies);

            if(mConfig.PerfDebug)
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

        SV_LOGGER.info("Esvee assembly complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private boolean loadJunctionFiles()
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

        String discStatsFilename = formDiscordantStatsFilename(mConfig.PrepDir, mConfig.sampleId(), mConfig.OutputId);
        DiscordantStats discordantStats = loadDiscordantStats(discStatsFilename);

        int minJunctionFrags = MIN_JUNCTION_SUPPORT;
        int minHotspotFrags = MIN_HOTSPOT_JUNCTION_SUPPORT;
        int minDiscordantFrags = DISCORDANT_GROUP_MIN_FRAGMENTS_SHORT;

        double discordantRate = discordantStats.discordantRate();

        if(discordantRate >= mConfig.DiscordantRateIncrement)
        {
            minJunctionFrags += (int)floor(discordantRate / mConfig.DiscordantRateIncrement) * DISC_RATE_JUNC_INCREMENT;
            minHotspotFrags += (int)floor(discordantRate / mConfig.DiscordantRateIncrement) * DISC_RATE_JUNC_INCREMENT;
            minDiscordantFrags += (int)floor(discordantRate / mConfig.DiscordantRateIncrement) * DISC_RATE_DISC_ONLY_INCREMENT;

            SV_LOGGER.info("raised min fragments(hotspot={} junction={} disc-only={}) for discordantRate({})",
                    minHotspotFrags, minJunctionFrags, minDiscordantFrags, format("%.3f", discordantRate));
        }

        for(String junctionFile : mConfig.JunctionFiles)
        {
            Map<String,List<Junction>> newJunctionsMap = Junction.loadJunctions(
                    junctionFile, mConfig.SpecificChrRegions, minJunctionFrags, minHotspotFrags, minDiscordantFrags);

            if(newJunctionsMap == null)
                return false;

            Junction.mergeJunctions(mChrJunctionsMap, newJunctionsMap);
        }

        // if(AssemblyConfig.DevDebug && !validateJunctionMap(mChrJunctionsMap))
        //    System.exit(1);

        if(mConfig.JunctionFiles.size() > 1)
        {
            SV_LOGGER.debug("merged into {} junctions", mChrJunctionsMap.values().stream().mapToInt(x -> x.size()).sum());
        }

        return true;
    }

    private void loadFragmentLengthBounds()
    {
        String fragmentLengthFile = formFragmentLengthDistFilename(mConfig.PrepDir, mConfig.sampleId(), mConfig.OutputId);

        if(!Files.exists(Paths.get(fragmentLengthFile)))
        {
            SV_LOGGER.error("missing fragment length file: {}", fragmentLengthFile);
            System.exit(1);
        }

        FragmentLengthBounds fragmentLengthBounds = FragmentSizeDistribution.loadFragmentLengthBounds(fragmentLengthFile);

        if(fragmentLengthBounds.isValid())
        {
            SV_LOGGER.info("fragment length bounds(min={} max={})",
                    fragmentLengthBounds.LowerBound, fragmentLengthBounds.UpperBound);

            MAX_OBSERVED_CONCORDANT_FRAG_LENGTH = fragmentLengthBounds.UpperBound;
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

        if(mResultsWriter.assemblyWriter() != null)
        {
            for(JunctionGroupAssembler jgAssembler : primaryAssemblyTasks)
            {
                for(JunctionAssembly decoyAssembly : jgAssembler.decoyAssemblies())
                {
                    decoyAssembly.setOutcome(DECOY, true);
                    mResultsWriter.assemblyWriter().writeAssembly(decoyAssembly);
                }

                jgAssembler.decoyAssemblies().clear();
            }
        }

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
                candidateReadCount += assembly.unmappedReads().size();
            }
        }

        ReadStats combinedReadStats = new ReadStats();
        primaryAssemblyTasks.forEach(x -> combinedReadStats.merge(x.readStats()));

        SV_LOGGER.info("created {} junction assemblies from reads(junc={} candidate={})",
                assemblyCount, junctionReadCount, candidateReadCount);

        SV_LOGGER.info("extracted read stats: {}", combinedReadStats);

        if(AssemblyConfig.DevDebug && combinedReadStats.IdenticalSupplementaries > 0)
        {
            SV_LOGGER.info("filtered identical supplementaries({})", combinedReadStats.IdenticalSupplementaries);
        }

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

        int taskCount = mBamReaders.size();

        List<PhaseSetTask> phaseSetTasks = PhaseSetTask.createThreadTasks(
                mConfig, phaseGroupBuilder.phaseGroups(), mBamReaders, taskCount, threadTasks);

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
        Alignment alignment = new Alignment(mConfig, new BwaAligner(mConfig.RefGenomeImageFile));
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

        for(PhaseGroup phaseGroup : phaseGroups)
        {
            allAssemblies.addAll(phaseGroup.derivedAssemblies());

            // add link assemblies into the same assembly alignment
            for(PhaseSet phaseSet : phaseGroup.phaseSets())
            {
                if(phaseSet.hasValidAssemblyAlignment())
                    assemblyAlignments.add(phaseSet.assemblyAlignment());
            }

            // and then add any assemblies not in a phase set into their own for alignment
            for(JunctionAssembly assembly : phaseGroup.assemblies())
            {
                if(assembly.phaseSet() == null && !skipUnlinkedJunctionAssembly(assembly))
                {
                    AssemblyAlignment assemblyAlignment = new AssemblyAlignment(assembly);
                    assemblyAlignments.add(assemblyAlignment);
                }
            }
        }

        Collections.sort(allAssemblies, Comparator.comparing(x -> x.junction()));

        for(JunctionAssembly assembly : allAssemblies)
        {
            assembly.setId(assemblyId++);
        }

        int assemblyAlignmentId = 0;

        for(AssemblyAlignment assemblyAlignment : assemblyAlignments)
        {
            assemblyAlignment.setId(assemblyAlignmentId++);
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

    private void writeVariants(final List<AssemblyAlignment> assemblyAlignments, final List<Breakend> breakends)
    {
        if(mConfig.WriteTypes.contains(WriteType.VCF))
        {
            SV_LOGGER.debug("writing variant VCF with {} breakends", breakends.size());

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

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        AssemblyConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        AssemblyApplication assembly = new AssemblyApplication(configBuilder);
        assembly.run();
        assembly.close();
    }
}
