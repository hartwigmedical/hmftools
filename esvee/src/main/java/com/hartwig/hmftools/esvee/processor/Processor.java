package com.hartwig.hmftools.esvee.processor;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.BAM_READ_JUNCTION_BUFFER;
import static com.hartwig.hmftools.esvee.assembly.JunctionGroupAssembler.mergePrimaryAssemblies;
import static com.hartwig.hmftools.esvee.common.JunctionGroup.buildJunctionGroups;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
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
import com.hartwig.hmftools.esvee.assembly.AssemblyMerger;
import com.hartwig.hmftools.esvee.assembly.JunctionGroupAssembler;
import com.hartwig.hmftools.esvee.assembly.SupportChecker;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.common.JunctionGroup;
import com.hartwig.hmftools.esvee.common.RegionOfInterest;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.WriteType;
import com.hartwig.hmftools.esvee.assembly.AssemblyExtender;
import com.hartwig.hmftools.esvee.common.SampleSupport;
import com.hartwig.hmftools.esvee.common.VariantCall;
import com.hartwig.hmftools.esvee.output.ResultsWriter;
import com.hartwig.hmftools.esvee.read.BamReader;
import com.hartwig.hmftools.esvee.sequence.AlignedAssembly;
import com.hartwig.hmftools.esvee.sequence.ExtendedAssembly;
import com.hartwig.hmftools.esvee.sequence.GappedAssembly;
import com.hartwig.hmftools.esvee.sequence.PrimaryAssembly;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.sequence.Sequence;
import com.hartwig.hmftools.esvee.sequence.SupportedAssembly;
import com.hartwig.hmftools.esvee.output.VcfWriter;
import com.hartwig.hmftools.esvee.html.SummaryPageGenerator;
import com.hartwig.hmftools.esvee.html.VariantCallPageGenerator;
import com.hartwig.hmftools.esvee.util.ParallelMapper;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

public class Processor
{
    private final SvConfig mConfig;
    private final Context mContext;
    private final ResultsWriter mResultsWriter;
    private final HomologySlider mHomologySlider;
    private final SupportChecker mSupportChecker;

    private final Map<String,List<Junction>> mChrJunctionsMap;

    private final List<PerformanceCounter> mPerfCounters;

    private final OverallCounters mCounters;

    public Processor(final SvConfig config, final Context context, final ResultsWriter resultsWriter)
    {
        mConfig = config;
        mContext = context;
        mResultsWriter = resultsWriter;
        mHomologySlider = new HomologySlider(mContext.ReferenceGenome);
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
        catch(final Exception e)
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

        mPerfCounters.add(ThreadTask.mergePerfCounters(primaryAssemblyTasks.stream().collect(Collectors.toList())));
        threadTasks.clear();

        List<PrimaryAssembly> primaryAssemblies = mergePrimaryAssemblies(primaryAssemblyTasks);

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

        if(!mConfig.OtherDebug)
        {
            primaryAssemblies.clear();
        }

        // FIXME: now clear all cached reads not assigned as support
        junctionGroups.forEach(x -> x.clearCandidateReads());

        SV_LOGGER.info("created {} extended assemblies", extendedAssemblies.size());

        // Primary phasing
        final List<Set<ExtendedAssembly>> primaryPhaseSets = PrimaryPhasing.run(extendedAssemblies);

        SV_LOGGER.info("created {} primary phase sets", primaryPhaseSets.size());

        if(!mConfig.OtherDebug)
            extendedAssemblies.clear();

        // Phased assembly merging
        final List<Set<ExtendedAssembly>> mergedPhaseSets = ParallelMapper.mapWithProgress(
                mContext.Executor, primaryPhaseSets, this::primaryPhasedMerging); // FIXME: switch to use AssemblyMerger

        SV_LOGGER.info("merged primary phase sets");

        // Secondary phasing
        final List<Set<ExtendedAssembly>> secondaryPhaseSets = SecondaryPhasing.run(mergedPhaseSets);
        SV_LOGGER.info("created {} secondary phase sets", secondaryPhaseSets.size());

        // Secondary merging
        final List<GappedAssembly> mergedSecondaries = new ArrayList<>();

        // FIXME: Gapped assemblies
        //for(int i = 0; i < secondaryPhaseSets.size(); i++)
        //    mergedSecondaries.add(createGapped(secondaryPhaseSets.get(i), i));
        for(int i = 0; i < secondaryPhaseSets.size(); i++)
        {
            final Set<ExtendedAssembly> phaseSet = secondaryPhaseSets.get(i);
            int j = 0;
            for(ExtendedAssembly assembly : phaseSet)
            {
                final GappedAssembly newAssembly = new GappedAssembly(String.format("Assembly%s-%s", i, j++), List.of(assembly));
                newAssembly.addErrata(assembly.getAllErrata());
                for(Map.Entry<Read, Integer> entry : assembly.getSupport())
                    newAssembly.addEvidenceAt(entry.getKey(), entry.getValue());

                mergedSecondaries.add(newAssembly);
            }
        }

        SV_LOGGER.info("merged secondaries");

        // alignment
        final List<AlignedAssembly> aligned = ParallelMapper.mapWithProgress(
                mContext.Executor, mergedSecondaries, mContext.Aligner::align);

        SV_LOGGER.info("created {} alignments", aligned.size());

        // Left sliding (we will find mid-points after calling + de-duping)
        final List<AlignedAssembly> homologised = ParallelMapper.mapWithProgress(
                mContext.Executor, aligned, mHomologySlider::slideHomology);

        SV_LOGGER.info("processed homology");

        // Support scanning
        // FIXME: consider not rescanning any aligned assembly that matches exactly or closely to the original
        /*
        SupportScanner supportScanner = new SupportScanner(mContext, Counters.ExtraScannedSupport);

        ParallelMapper.mapWithProgress(
                mContext.Executor, homologised, supportScanner::tryRescanSupport);

        SV_LOGGER.info("rescanned support, adding {} new reads", Counters.ExtraScannedSupport.formatValue());
        */


        // Calling
        final List<VariantCall> variants = new VariantCaller(mContext.Executor).callVariants(homologised);
        SV_LOGGER.info("called {} variants", variants.size());

        // variant deduplication
        final VariantDeduplication deduplicator = new VariantDeduplication(mContext, mCounters.VariantDeduplicationCounters);
        final List<VariantCall> deduplicated = deduplicator.deduplicate(variants);

        SV_LOGGER.info("{} variants remaining after deduplication", deduplicated.size());
        deduplicated.removeIf(variant -> variant.supportingFragments().isEmpty());
        SV_LOGGER.info("{} variants remaining after removing unsubstantiated", deduplicated.size());

        final long lowQualityVariants = deduplicated.stream()
                .filter(variant -> variant.quality() < SvConstants.VCFLOWQUALITYTHRESHOLD)
                .count();
        SV_LOGGER.info("{} low-quality variants found", lowQualityVariants);

        final long lowSupportVariants = deduplicated.stream()
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
            catch(final Exception ex)
            {
                SV_LOGGER.error("Failed to generate HTML for {}", call, ex);
            }
        }

        try
        {
            SummaryPageGenerator.generatePage(mConfig.HtmlOutputDir, mCounters, variants);
        }
        catch(final Exception ex)
        {
            SV_LOGGER.error("Failure while generating summary HTML", ex);
        }
    }

    private void writeVCF(final List<VariantCall> variants)
    {
        final List<String> sampleNames = variants.stream()
                .flatMap(call -> call.sampleSupport().stream().map(SampleSupport::sampleName))
                .distinct()
                .sorted()
                .collect(Collectors.toList());

        final VcfWriter writer = new VcfWriter(mContext, sampleNames);
        for(VariantCall call : variants)
        {
            try
            {
                writer.append(call);
            }
            catch(final Exception ex)
            {
                SV_LOGGER.error("Failure while appending to call VCF: {}", call, ex);
            }
        }
        writer.close();
    }

    private Set<ExtendedAssembly> primaryPhasedMerging(final Set<ExtendedAssembly> primaryPhaseSet)
    {
        try
        {
            final Set<ExtendedAssembly> result = new HashSet<>(primaryPhaseSet);
            final Set<Pair<ExtendedAssembly, ExtendedAssembly>> checked = new HashSet<>();
            boolean merged = true;
            while(merged)
            {
                merged = false;

                loopHead:
                for(ExtendedAssembly left : result)
                    for(ExtendedAssembly right : result)
                    {
                        if(left == right)
                            continue;
                        if(!checked.add(Pair.of(left, right)))
                            continue;

                        final int minOverlap = Math.min(30, Math.min(left.getLength(), right.getLength()));
                        @Nullable
                        Integer index = mSupportChecker.AssemblySupport.supportIndex(left, right, minOverlap);
                        if(index != null)
                            index = mSupportChecker.AssemblySupport.bestSupportIndex(left, right, minOverlap);
                        final ExtendedAssembly mergedAssembly;
                        if(index == null)
                        {
                            final ExtendedAssembly flippedRight = right.flipStrand();
                            index = mSupportChecker.AssemblySupport.supportIndex(left, flippedRight, minOverlap);
                            if(index != null)
                                index = mSupportChecker.AssemblySupport.bestSupportIndex(left, right, minOverlap);
                            if(index == null)
                                continue;

                            mergedAssembly = merge(left, flippedRight, index);
                        }
                        else
                            mergedAssembly = merge(left, right, index);

                        result.remove(left);
                        result.remove(right);
                        result.add(mergedAssembly);

                        merged = true;
                        break loopHead;
                    }
            }

            return result;
        }
        catch(final Throwable throwable)
        {
            SV_LOGGER.warn("Failure during phased assembly merging with group of size {}", primaryPhaseSet.size(), throwable);
            SV_LOGGER.warn("{}", RegionOfInterest.tryMerge(
                    primaryPhaseSet.stream()
                            .flatMap(assembly -> assembly.getSupport().stream())
                            .map(Map.Entry::getKey)
                            .filter(record -> !record.isUnmapped())
                            .map(record -> new RegionOfInterest(record.getChromosome(), record.getAlignmentStart(), record.getAlignmentEnd()))
                            .collect(Collectors.toList())
            ));
            return null;
        }
    }

    private ExtendedAssembly merge(final ExtendedAssembly left, final ExtendedAssembly right, final int supportIndex)
    {
        left.markDecompositionStale();
        right.markDecompositionStale();
        final Sequence mergedSequence = SequenceMerger.merge(left, right, supportIndex);

        final ExtendedAssembly merged = new ExtendedAssembly(left.Name, mergedSequence.getBasesString(), left.Source);
        left.Diagrams.forEach(merged::addDiagrams);

        left.getSupportRecords().forEach(support -> merged.tryAddSupport(mSupportChecker, support));
        right.getSupportRecords().forEach(support -> merged.tryAddSupport(mSupportChecker, support));

        merged.addErrata(left.getAllErrata());
        merged.addErrata(right.getAllErrata());

        return merged;
    }

    private void reAddSupport(final SupportedAssembly merged, final SupportedAssembly old)
    {
        final int offset = merged.Assembly.indexOf(old.Assembly);
        for(Map.Entry<Read, Integer> entry : old.getSupport())
        {
            final Read potentialSupport = entry.getKey();
            if(offset != -1)
            {
                final int oldSupportIndex = entry.getValue();
                if(mSupportChecker.AssemblySupport.supportsAt(merged, potentialSupport, oldSupportIndex + offset))
                {
                    merged.addEvidenceAt(potentialSupport, oldSupportIndex + offset);
                    continue;
                }
            }
            merged.tryAddSupport(mSupportChecker, potentialSupport);
        }
    }

    private List<ExtendedAssembly> order(final Collection<ExtendedAssembly> assemblies)
    {
        // FIXME: Correctly order these
        if(assemblies.size() > 1)
            SV_LOGGER.warn("Found more than 1 assembly ({}) while creating gapped ({})", assemblies.size(),
                    assemblies.stream().map(assembly -> assembly.Name).collect(Collectors.toList()));

        final Map<ExtendedAssembly, Map<ExtendedAssembly, Long>> leftWise = new IdentityHashMap<>();
        final Map<ExtendedAssembly, Map<ExtendedAssembly, Long>> rightWise = new IdentityHashMap<>();
        for(ExtendedAssembly first : assemblies)
            for(ExtendedAssembly second : assemblies)
            {
                if(first == second)
                    continue;



            }

        return new ArrayList<>(assemblies);
    }

    public GappedAssembly createGapped(final Collection<ExtendedAssembly> assemblies, final int index)
    {
        final GappedAssembly gappedAssembly = new GappedAssembly("Assembly" + index, order(assemblies));

        for(ExtendedAssembly assembly : assemblies)
            for(Read support : assembly.getSupportRecords())
                if(!gappedAssembly.tryAddSupport(mSupportChecker, support))
                    SV_LOGGER.info("Failed to add support for assembly {}: {}", gappedAssembly.Name, support.getName());

        return gappedAssembly;
    }
}
