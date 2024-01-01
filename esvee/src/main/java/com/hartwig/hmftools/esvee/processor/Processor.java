package com.hartwig.hmftools.esvee.processor;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.BAM_READ_JUNCTION_BUFFER;
import static com.hartwig.hmftools.esvee.assembly.Aligner.mergeAlignedAssemblies;
import static com.hartwig.hmftools.esvee.assembly.JunctionGroupAssembler.mergePrimaryAssemblies;
import static com.hartwig.hmftools.esvee.assembly.PhasedMerger.createThreadTasks;
import static com.hartwig.hmftools.esvee.assembly.PhasedMerger.mergePhasedResults;
import static com.hartwig.hmftools.esvee.common.JunctionGroup.buildJunctionGroups;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.esvee.Context;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.assembly.Aligner;
import com.hartwig.hmftools.esvee.assembly.AssemblyMerger;
import com.hartwig.hmftools.esvee.assembly.JunctionGroupAssembler;
import com.hartwig.hmftools.esvee.assembly.PhasedMerger;
import com.hartwig.hmftools.esvee.assembly.SupportChecker;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.common.JunctionGroup;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.WriteType;
import com.hartwig.hmftools.esvee.assembly.AssemblyExtender;
import com.hartwig.hmftools.esvee.common.SampleSupport;
import com.hartwig.hmftools.esvee.common.ThreadTask;
import com.hartwig.hmftools.esvee.common.VariantCall;
import com.hartwig.hmftools.esvee.output.ResultsWriter;
import com.hartwig.hmftools.esvee.read.BamReader;
import com.hartwig.hmftools.esvee.sequence.AlignedAssembly;
import com.hartwig.hmftools.esvee.sequence.ExtendedAssembly;
import com.hartwig.hmftools.esvee.sequence.GappedAssembly;
import com.hartwig.hmftools.esvee.sequence.PrimaryAssembly;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.output.VcfWriter;
import com.hartwig.hmftools.esvee.html.SummaryPageGenerator;
import com.hartwig.hmftools.esvee.html.VariantCallPageGenerator;
import com.hartwig.hmftools.esvee.sequence.ReadSupport;
import com.hartwig.hmftools.esvee.sequence.Sequence;
import com.hartwig.hmftools.esvee.sequence.SupportedAssembly;

public class Processor
{
    private final SvConfig mConfig;
    private final Context mContext;
    private final ResultsWriter mResultsWriter;
    private final SupportChecker mSupportChecker;

    private final Map<String,List<Junction>> mChrJunctionsMap;

    private final List<PerformanceCounter> mPerfCounters;

    private final OverallCounters mCounters;

    public Processor(final SvConfig config, final Context context, final ResultsWriter resultsWriter)
    {
        mConfig = config;
        mContext = context;
        mResultsWriter = resultsWriter;
        mSupportChecker = new SupportChecker();
        mChrJunctionsMap = Maps.newHashMap();
        mCounters = new OverallCounters();
        mPerfCounters = Lists.newArrayList();
    }

    public boolean loadJunctionFiles()
    {
        for(String junctionFile : mConfig.JunctionFiles)
        {
            Map<String,List<Junction>> newJunctionsMap = Junction.loadJunctions(junctionFile);

            if(newJunctionsMap == null)
                return false;

            Junction.mergeJunctions(mChrJunctionsMap, newJunctionsMap);
        }

        if(mConfig.JunctionFiles.size() > 1)
        {
            SV_LOGGER.debug("merged into {} junctions", mChrJunctionsMap.values().stream().mapToInt(x -> x.size()).sum());
        }

        return true;
    }

    public List<VariantCall> run()
    {
        try
        {
            Map<String,List<JunctionGroup>> junctionGroupMap = Maps.newHashMap();

            for(Map.Entry<String,List<Junction>> entry : mChrJunctionsMap.entrySet())
            {
                junctionGroupMap.put(entry.getKey(), buildJunctionGroups(entry.getValue(), BAM_READ_JUNCTION_BUFFER));
            }

            List<Junction> allJunctions = Lists.newArrayList();
            mChrJunctionsMap.values().forEach(x -> allJunctions.addAll(x));

            List<VariantCall> variantCalls = run(allJunctions, junctionGroupMap);

            return variantCalls;
        }
        catch(Exception e)
        {
            SV_LOGGER.error("process run error: {}", e.toString());
            e.printStackTrace();
            System.exit(1);
        }

        return null;
    }

    public List<VariantCall> run(final List<Junction> junctions, final Map<String,List<JunctionGroup>> junctionGroupMap)
    {
        List<JunctionGroup> junctionGroups = Lists.newArrayList();
        junctionGroupMap.values().forEach(x -> junctionGroups.addAll(x));

        mCounters.JunctionsProcessed.add(junctions.size());

        Collections.sort(junctionGroups);

        List<BamReader> bamReaders = Lists.newArrayList();
        List<Thread> threadTasks = new ArrayList<>();

        int taskCount = min(junctionGroups.size(), mConfig.Threads);

        for(int i = 0; i < taskCount; ++i)
        {
            BamReader bamReader = new BamReader(mConfig);
            bamReaders.add(bamReader);
        }

        // Primary Junction Assembly
        List<JunctionGroupAssembler> primaryAssemblyTasks = JunctionGroupAssembler.createThreadTasks(
                junctionGroups, bamReaders, mConfig, taskCount, threadTasks);

        if(!runThreadTasks(threadTasks))
            System.exit(1);

        List<PrimaryAssembly> primaryAssemblies = mergePrimaryAssemblies(primaryAssemblyTasks);

        mPerfCounters.add(ThreadTask.mergePerfCounters(primaryAssemblyTasks.stream().collect(Collectors.toList())));
        primaryAssemblyTasks.clear();
        threadTasks.clear();

        SV_LOGGER.info("created {} primary assemblies", primaryAssemblies.size());

        int totalCachedReads = junctionGroups.stream().mapToInt(x -> x.candidateReadCount()).sum();
        SV_LOGGER.debug("cached read count({}) from {} junction groups", totalCachedReads, junctionGroups.size());

        // FIXME: is this necessary? should clear when no longer required, and now factor in junction groups
        if(!mConfig.OtherDebug)
        {
            junctions.clear();
        }

        // Inter-junction deduplication
        // FIXME: surely can just be within junction groups or at least same chromosome?
        AssemblyMerger assemblyMerger = new AssemblyMerger(); // only instantiated because of SupportChecker

        List<PrimaryAssembly> mergedPrimaryAssemblies = assemblyMerger.consolidatePrimaryAssemblies(primaryAssemblies);

        SV_LOGGER.info("reduced to {} assemblies", mergedPrimaryAssemblies.size());

        if(!mConfig.OtherDebug)
            primaryAssemblies.clear();

        // CHECK
        // if(mConfig.dropGermline())
        //    primaryAssemblies.removeIf(assembly -> assembly.getSupportRecords().stream().anyMatch(Record::isGermline));

        // assembly extension
        List<AssemblyExtender> assemblyExtenderTasks = AssemblyExtender.createThreadTasks(
                junctionGroupMap, mergedPrimaryAssemblies, mConfig, taskCount, threadTasks);

        if(!runThreadTasks(threadTasks))
            System.exit(1);

        List<ExtendedAssembly> extendedAssemblies = Lists.newArrayList();
        assemblyExtenderTasks.forEach(x -> extendedAssemblies.addAll(x.extendedAssemblies()));

        mPerfCounters.add(ThreadTask.mergePerfCounters(assemblyExtenderTasks.stream().collect(Collectors.toList())));
        assemblyExtenderTasks.clear();
        threadTasks.clear();

        if(!mConfig.OtherDebug)
        {
            primaryAssemblies.clear();
        }

        // FIXME: now clear all cached reads not assigned as support
        junctionGroups.forEach(x -> x.clearCandidateReads());

        SV_LOGGER.info("created {} extended assemblies", extendedAssemblies.size());

        // Primary phasing
        List<Set<ExtendedAssembly>> primaryPhaseSets = PrimaryPhasing.run(extendedAssemblies);

        SV_LOGGER.info("created {} primary phase sets", primaryPhaseSets.size());

        if(!mConfig.OtherDebug)
            extendedAssemblies.clear();

        // Phased assembly merging
        List<PhasedMerger> phasedMergerTasks = createThreadTasks(primaryPhaseSets, mConfig, taskCount, threadTasks);

        if(!runThreadTasks(threadTasks))
            System.exit(1);

        List<Set<ExtendedAssembly>> mergedPhaseSets = mergePhasedResults(phasedMergerTasks);

        mPerfCounters.add(ThreadTask.mergePerfCounters(phasedMergerTasks.stream().collect(Collectors.toList())));
        phasedMergerTasks.clear();
        threadTasks.clear();

        //List<Set<ExtendedAssembly>> mergedPhaseSets = ParallelMapper.mapWithProgress(
        //        mContext.Executor, primaryPhaseSets, assemblyMerger::primaryPhasedMerging);

        SV_LOGGER.info("merged to {} primary phase sets", mergedPhaseSets.size());

        // secondary phasing, not multi-threaded but could it be?
        List<Set<ExtendedAssembly>> secondaryPhaseSets = SecondaryPhasing.run(mergedPhaseSets);
        SV_LOGGER.info("created {} secondary phase sets", secondaryPhaseSets.size());

        // Secondary merging
        List<GappedAssembly> mergedSecondaries = Lists.newArrayList();

        // FIXME: Gapped assemblies
        //for(int i = 0; i < secondaryPhaseSets.size(); i++)
        //    mergedSecondaries.add(createGapped(secondaryPhaseSets.get(i), i));
        for(int i = 0; i < secondaryPhaseSets.size(); i++)
        {
            Set<ExtendedAssembly> phaseSet = secondaryPhaseSets.get(i);
            int j = 0;
            for(ExtendedAssembly assembly : phaseSet)
            {
                GappedAssembly newAssembly = new GappedAssembly(String.format("Assembly%s-%s", i, j++), List.of(assembly));
                newAssembly.addErrata(assembly.getAllErrata());

                assembly.readSupport().forEach(x -> newAssembly.addEvidenceAt(x.Read, x.Index));

                mergedSecondaries.add(newAssembly);
            }
        }

        SV_LOGGER.info("merged to {} secondaries", mergedSecondaries.size());

        // alignment
        // List<AlignedAssembly> aligned = ParallelMapper.mapWithProgress(mContext.Executor, mergedSecondaries, mContext.Aligner::align);

        // now also handles homology sliding
        // ParallelMapper.mapWithProgress(mContext.Executor, homologised, supportScanner::rescanAssemblySupport);

        List<Aligner> alignerTasks = Aligner.createThreadTasks(mergedSecondaries, mConfig, taskCount, threadTasks);

        if(!runThreadTasks(threadTasks))
            System.exit(1);

        List<AlignedAssembly> alignedAssemblies = mergeAlignedAssemblies(alignerTasks);

        mPerfCounters.add(ThreadTask.mergePerfCounters(alignerTasks.stream().collect(Collectors.toList())));
        alignerTasks.clear();
        threadTasks.clear();

        SV_LOGGER.info("created {} alignments", alignedAssemblies.size());

        // Left sliding (we will find mid-points after calling + de-duping)
        //List<AlignedAssembly> homologised = ParallelMapper.mapWithProgress( mContext.Executor, aligned, mHomologySlider::slideHomology);
        // SV_LOGGER.info("processed homology");

        // Support scanning
        // FIXME: consider not rescanning any aligned assembly that matches exactly or closely to the original
        List<SupportScanner> supportScannerTasks = SupportScanner.createThreadTasks(
                alignedAssemblies, bamReaders, mConfig, taskCount, threadTasks);

        if(!runThreadTasks(threadTasks))
            System.exit(1);

        // existing alignment assemblies are maintained, just with any additional support
        // List<AlignedAssembly> rescannedAssemblies = SupportScanner.mergeRescannedAssemblies(supportScannerTasks);

        mPerfCounters.add(ThreadTask.mergePerfCounters(alignerTasks.stream().collect(Collectors.toList())));
        int rescannedAddedReadCount = supportScannerTasks.stream().mapToInt(x -> x.addedReadCount()).sum();
        alignerTasks.clear();
        threadTasks.clear();

        SV_LOGGER.info("rescanned support, adding {} new reads", rescannedAddedReadCount);


        // Calling
        List<VariantCall> variants = new VariantCaller(mContext.Executor).callVariants(alignedAssemblies);
        SV_LOGGER.info("called {} variants", variants.size());

        // variant deduplication
        VariantDeduplication deduplicator = new VariantDeduplication(mContext, mCounters.VariantDeduplicationCounters);
        List<VariantCall> deduplicated = deduplicator.deduplicate(variants);

        SV_LOGGER.info("{} variants remaining after deduplication", deduplicated.size());
        deduplicated.removeIf(variant -> variant.supportingFragments().isEmpty());
        SV_LOGGER.info("{} variants remaining after removing unsubstantiated", deduplicated.size());

        long lowQualityVariants = deduplicated.stream()
                .filter(variant -> variant.quality() < SvConstants.VCFLOWQUALITYTHRESHOLD)
                .count();
        SV_LOGGER.info("{} low-quality variants found", lowQualityVariants);

        long lowSupportVariants = deduplicated.stream()
                .filter(variant -> variant.quality() >= SvConstants.VCFLOWQUALITYTHRESHOLD)
                .filter(variant -> variant.supportingFragments().size() < SvConstants.MIN_READS_SUPPORT_ASSEMBLY)
                .count();
        SV_LOGGER.info("{} low-support variants found (excl low-quality)", lowSupportVariants);

        /* CHECK
        if(mConfig.dropGermline())
        {
            deduplicated.removeIf(VariantCall::isGermline);
            SV_LOGGER.info("{} variants remaining after dropping those with germline support", deduplicated.size());
        }
        */

        if(!mContext.Problems.isEmpty())
        {
            SV_LOGGER.warn("encountered {} problems", mContext.Problems.size());

            if(mContext.Problems.size() < 50)
            {
                for(Problem problem : mContext.Problems)
                {
                    SV_LOGGER.warn("{}", problem);
                }
            }
        }

        if(mConfig.writeHtmlFiles())
            writeHTMLSummaries(deduplicated);

        writeVCF(deduplicated);

        if(mConfig.WriteTypes.contains(WriteType.BREAKEND_TSV))
        {
            deduplicated.forEach(x -> mResultsWriter.writeVariant(x));
        }

        mResultsWriter.writeVariantAssemblyBamRecords(deduplicated);

        mPerfCounters.forEach(x -> x.logStats());

        return deduplicated;
    }

    private void writeHTMLSummaries(final List<VariantCall> variants)
    {
        int summariesWritten = 0;
        for(VariantCall call : variants)
        {
            if(summariesWritten++ > SvConstants.MAX_HTML_SUMMARIES)
            {
                SV_LOGGER.warn("Not writing further HTML summaries -- limit reached. Increase -max_html_summaries to see more.");
                break;
            }

            try
            {
                VariantCallPageGenerator.generatePage(mConfig.HtmlOutputDir, mContext.ReferenceGenome, mSupportChecker, call);
            }
            catch(Exception ex)
            {
                SV_LOGGER.error("Failed to generate HTML for {}", call, ex);
            }
        }

        try
        {
            SummaryPageGenerator.generatePage(mConfig.HtmlOutputDir, mCounters, variants);
        }
        catch(Exception ex)
        {
            SV_LOGGER.error("Failure while generating summary HTML", ex);
        }
    }

    private void writeVCF(final List<VariantCall> variants)
    {
        List<String> sampleNames = variants.stream()
                .flatMap(call -> call.sampleSupport().stream().map(SampleSupport::sampleName))
                .distinct()
                .sorted()
                .collect(Collectors.toList());

        VcfWriter writer = new VcfWriter(mContext, sampleNames);
        for(VariantCall call : variants)
        {
            try
            {
                writer.append(call);
            }
            catch(Exception ex)
            {
                SV_LOGGER.error("Failure while appending to call VCF: {}", call, ex);
            }
        }
        writer.close();
    }

    // CHECK: these aren't called in original SvAssembly code, what is their purpose?
    public GappedAssembly createGapped(final Collection<ExtendedAssembly> assemblies, final int index)
    {
        GappedAssembly gappedAssembly = new GappedAssembly("Assembly" + index, orderExtendedAssemblies(assemblies));

        for(ExtendedAssembly assembly : assemblies)
        {
            for(Read support : assembly.supportingReads())
            {
                if(!gappedAssembly.tryAddSupport(mSupportChecker, support))
                {
                    SV_LOGGER.info("Failed to add support for assembly {}: {}", gappedAssembly.Name, support.getName());
                }
            }
        }

        return gappedAssembly;
    }

    private List<ExtendedAssembly> orderExtendedAssemblies(final Collection<ExtendedAssembly> assemblies)
    {
        // FIXME: Correctly order these
        if(assemblies.size() > 1)
            SV_LOGGER.warn("Found more than 1 assembly ({}) while creating gapped ({})", assemblies.size(),
                    assemblies.stream().map(assembly -> assembly.Name).collect(Collectors.toList()));

        Map<ExtendedAssembly, Map<ExtendedAssembly, Long>> leftWise = new IdentityHashMap<>();
        Map<ExtendedAssembly, Map<ExtendedAssembly, Long>> rightWise = new IdentityHashMap<>();
        for(ExtendedAssembly first : assemblies)
            for(ExtendedAssembly second : assemblies)
            {
                if(first == second)
                    continue;

            }

        return new ArrayList<>(assemblies);
    }

    public AlignedAssembly mergeAlignedAssembly(
            final AlignedAssembly left, final AlignedAssembly right, final int supportIndex,
            final HomologySlider homologySlider, final Aligner aligner)
    {
        final Sequence mergedSequence = SequenceMerger.merge(left, right, supportIndex);

        final ExtendedAssembly merged = new ExtendedAssembly(left.Name, mergedSequence.getBasesString(), left.Source);
        left.Source.Sources.get(0).Diagrams.forEach(merged::addDiagrams);

        final GappedAssembly gapped = new GappedAssembly(merged.Name, List.of(merged));
        reAddSupport(gapped, left);
        reAddSupport(gapped, right);

        return homologySlider.slideHomology(aligner.align(gapped));
    }

    private void reAddSupport(final SupportedAssembly merged, final SupportedAssembly old)
    {
        final int offset = merged.Assembly.indexOf(old.Assembly);

        for(ReadSupport readSupport : old.readSupport())
        {
            Read potentialSupport = readSupport.Read;

            if(offset != -1)
            {
                int oldSupportIndex = readSupport.Index;
                if(mSupportChecker.AssemblySupport.supportsAt(merged, potentialSupport, oldSupportIndex + offset))
                {
                    merged.addEvidenceAt(potentialSupport, oldSupportIndex + offset);
                    continue;
                }
            }
            merged.tryAddSupport(mSupportChecker, potentialSupport);
        }
    }
}
