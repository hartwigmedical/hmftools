package com.hartwig.hmftools.esvee.processor;

import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BiFunction;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.Context;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.RegionOfInterest;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.WriteType;
import com.hartwig.hmftools.esvee.assembly.AssemblyExtender;
import com.hartwig.hmftools.esvee.assembly.PrimaryAssembler;
import com.hartwig.hmftools.esvee.common.SampleSupport;
import com.hartwig.hmftools.esvee.common.VariantCall;
import com.hartwig.hmftools.esvee.output.ResultsWriter;
import com.hartwig.hmftools.esvee.models.AlignedAssembly;
import com.hartwig.hmftools.esvee.models.ExtendedAssembly;
import com.hartwig.hmftools.esvee.models.GappedAssembly;
import com.hartwig.hmftools.esvee.models.PrimaryAssembly;
import com.hartwig.hmftools.esvee.models.Record;
import com.hartwig.hmftools.esvee.models.Sequence;
import com.hartwig.hmftools.esvee.models.SupportedAssembly;
import com.hartwig.hmftools.esvee.output.VcfWriter;
import com.hartwig.hmftools.esvee.html.SummaryPageGenerator;
import com.hartwig.hmftools.esvee.html.VariantCallPageGenerator;
import com.hartwig.hmftools.esvee.util.NaturalSortComparator;
import com.hartwig.hmftools.esvee.util.ParallelMapper;
import com.hartwig.hmftools.esvee.util.CommonUtils;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

public class Processor
{
    private final SvConfig mConfig;
    private final Context mContext;
    private final ResultsWriter mResultsWriter;
    private final HomologySlider mHomologySlider;

    private final Map<String,List<Junction>> mChrJunctionsMap;

    public final OverallCounters Counters = new OverallCounters();

    public Processor(final SvConfig config, final Context context, final ResultsWriter resultsWriter)
    {
        mConfig = config;
        mContext = context;
        mResultsWriter = resultsWriter;
        mHomologySlider = new HomologySlider(mContext.ReferenceGenome);
        mChrJunctionsMap = Maps.newHashMap();
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
            // for now merge all junctions into a single list - alternatives would be combined by proximity or chromosome
            List<Junction> allJunctions = Lists.newArrayList();
            mChrJunctionsMap.values().forEach(x -> allJunctions.addAll(x));

            List<VariantCall> variantCalls = run(allJunctions);

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

    public List<VariantCall> run(final List<Junction> junctions)
    {
        Counters.JunctionsProcessed.add(junctions.size());

        SV_LOGGER.info("starting primary assembly on {} junctions", junctions.size());

        // Primary Junction Assembly
        List<PrimaryAssemblyResult> primaryAssemblyResults = ParallelMapper.mapWithProgress(
                mContext.Executor, junctions,
                junction -> PrimaryAssembler.process(mContext, junction, Counters.PrimaryAssemblerCounters));

        final int primaryAssemblyCount = primaryAssemblyResults.stream().mapToInt(r -> r.Assemblies.size()).sum();

        SV_LOGGER.info("created {} primary assemblies", primaryAssemblyCount);

        if(!mContext.Config.OtherDebug)
        {
            junctions.clear();
        }

        // Inter-junction deduplication
        List<PrimaryAssembly> primaryAssemblies = consolidateNearbyAssemblies(flatten(primaryAssemblyResults));

        SV_LOGGER.info("reduced to {} assemblies", primaryAssemblies.size());

        if(!mContext.Config.OtherDebug)
            primaryAssemblyResults.clear();

        // CHECK
        // if(mContext.Config.dropGermline())
        //    primaryAssemblies.removeIf(assembly -> assembly.getSupportRecords().stream().anyMatch(Record::isGermline));

        // Assembly Extension
        List<ExtendedAssembly> extendedAssemblies = ParallelMapper.mapWithProgress(
                        mContext.Executor, primaryAssemblies,
                        assembly -> AssemblyExtender.process(mContext, assembly, Counters.AssemblyExtenderCounters)).stream()
                .flatMap(Collection::stream)
                .collect(Collectors.toList());

        if(!mContext.Config.OtherDebug)
        {
            primaryAssemblies.clear();
        }

        SV_LOGGER.info("created {} extended assemblies", extendedAssemblies.size());

        // Primary phasing
        final List<Set<ExtendedAssembly>> primaryPhaseSets = PrimaryPhasing.run(extendedAssemblies);

        SV_LOGGER.info("created {} primary phase sets", primaryPhaseSets.size());

        if(!mContext.Config.OtherDebug)
            extendedAssemblies.clear();

        // Phased assembly merging
        final List<Set<ExtendedAssembly>> mergedPhaseSets = ParallelMapper.mapWithProgress(
                mContext.Executor, primaryPhaseSets, this::primaryPhasedMerging);

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
            for(final ExtendedAssembly assembly : phaseSet)
            {
                final GappedAssembly newAssembly = new GappedAssembly(String.format("Assembly%s-%s", i, j++), List.of(assembly));
                newAssembly.addErrata(assembly.getAllErrata());
                for(final Map.Entry<Record, Integer> entry : assembly.getSupport())
                    newAssembly.addEvidenceAt(entry.getKey(), entry.getValue());

                mergedSecondaries.add(newAssembly);
            }
        }

        SV_LOGGER.info("merged secondaries");

        // Alignment
        final List<AlignedAssembly> aligned = ParallelMapper.mapWithProgress(
                mContext.Executor, mergedSecondaries, mContext.Aligner::align);

        SV_LOGGER.info("created {} alignments", aligned.size());

        // Left sliding (we will find mid-points after calling + de-duping)
        final List<AlignedAssembly> homologised = ParallelMapper.mapWithProgress(
                mContext.Executor, aligned, mHomologySlider::slideHomology);

        SV_LOGGER.info("processed homology");

        // Support scanning
        final SupportScanner supportScanner = new SupportScanner(mContext, Counters.ExtraScannedSupport);

        ParallelMapper.mapWithProgress(
                mContext.Executor, homologised, supportScanner::tryRescanSupport);

        SV_LOGGER.info("rescanned support, adding {} new reads", Counters.ExtraScannedSupport.formatValue());

        // Calling
        final List<VariantCall> variants = new VariantCaller(mContext.Executor).callVariants(homologised);
        SV_LOGGER.info("called {} variants", variants.size());

        // Variant Deduplication
        final VariantDeduplication deduplicator = new VariantDeduplication(mContext, Counters.VariantDeduplicationCounters);
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
                .filter(variant -> variant.supportingFragments().size() < SvConstants.MINREADSTOSUPPORTASSEMBLY)
                .count();
        SV_LOGGER.info("{} low-support variants found (excl low-quality)", lowSupportVariants);

        /* CHECK
        if(mContext.Config.dropGermline())
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
                for(final Problem problem : mContext.Problems)
                    SV_LOGGER.warn("{}", problem);
            }
        }

        if(mContext.Config.writeHtmlFiles())
            writeHTMLSummaries(deduplicated);

        writeVCF(deduplicated);

        if(mContext.Config.WriteTypes.contains(WriteType.BREAKEND_TSV))
        {
            deduplicated.forEach(x -> mResultsWriter.writeVariant(x));
        }

        mResultsWriter.writeVariantAssemblyBamRecords(deduplicated);

        return deduplicated;
    }

    private void writeHTMLSummaries(final List<VariantCall> variants)
    {
        int summariesWritten = 0;
        for(final VariantCall call : variants)
        {
            if(summariesWritten++ > SvConstants.MAX_HTML_SUMMARIES)
            {
                SV_LOGGER.warn("Not writing further HTML summaries -- limit reached. Increase -max_html_summaries to see more.");
                break;
            }

            try
            {
                VariantCallPageGenerator.generatePage(mContext.Config.HtmlOutputDir, mContext.ReferenceGenome, mContext.SupportChecker, call);
            }
            catch(final Exception ex)
            {
                SV_LOGGER.error("Failed to generate HTML for {}", call, ex);
            }
        }

        try
        {
            SummaryPageGenerator.generatePage(mContext.Config.HtmlOutputDir, Counters, variants);
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
        for(final VariantCall call : variants)
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
                for(final ExtendedAssembly left : result)
                    for(final ExtendedAssembly right : result)
                    {
                        if(left == right)
                            continue;
                        if(!checked.add(Pair.of(left, right)))
                            continue;

                        final int minOverlap = Math.min(30, Math.min(left.getLength(), right.getLength()));
                        @Nullable
                        Integer index = mContext.SupportChecker.AssemblySupport.supportIndex(left, right, minOverlap);
                        if (index != null)
                            index = mContext.SupportChecker.AssemblySupport.bestSupportIndex(left, right, minOverlap);
                        final ExtendedAssembly mergedAssembly;
                        if(index == null)
                        {
                            final ExtendedAssembly flippedRight = right.flipStrand();
                            index = mContext.SupportChecker.AssemblySupport.supportIndex(left, flippedRight, minOverlap);
                            if (index != null)
                                index = mContext.SupportChecker.AssemblySupport.bestSupportIndex(left, right, minOverlap);
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

        left.getSupportRecords().forEach(support -> merged.tryAddSupport(mContext.SupportChecker, support));
        right.getSupportRecords().forEach(support -> merged.tryAddSupport(mContext.SupportChecker, support));

        merged.addErrata(left.getAllErrata());
        merged.addErrata(right.getAllErrata());

        return merged;
    }

    private AlignedAssembly merge(final AlignedAssembly left, final AlignedAssembly right, final int supportIndex)
    {
        final Sequence mergedSequence = SequenceMerger.merge(left, right, supportIndex);

        final ExtendedAssembly merged = new ExtendedAssembly(left.Name, mergedSequence.getBasesString(), left.Source);
        left.Source.Sources.get(0).Diagrams.forEach(merged::addDiagrams);

        final GappedAssembly gapped = new GappedAssembly(merged.Name, List.of(merged));
        reAddSupport(gapped, left);
        reAddSupport(gapped, right);

        return mHomologySlider.slideHomology(mContext.Aligner.align(gapped));
    }

    private void reAddSupport(final SupportedAssembly merged, final SupportedAssembly old)
    {
        final int offset = merged.Assembly.indexOf(old.Assembly);
        for(final Map.Entry<Record, Integer> entry : old.getSupport())
        {
            final Record potentialSupport = entry.getKey();
            if(offset != -1)
            {
                final int oldSupportIndex = entry.getValue();
                if(mContext.SupportChecker.AssemblySupport.supportsAt(merged, potentialSupport, oldSupportIndex + offset))
                {
                    merged.addEvidenceAt(potentialSupport, oldSupportIndex + offset);
                    continue;
                }
            }
            merged.tryAddSupport(mContext.SupportChecker, potentialSupport);
        }
    }

    public List<PrimaryAssembly> flatten(final List<PrimaryAssemblyResult> results)
    {
        return results.stream()
                .flatMap(r -> r.Assemblies.stream())
                .collect(Collectors.toList());
    }

    public List<PrimaryAssembly> consolidateNearbyAssemblies(final List<PrimaryAssembly> results)
    {
        final List<PrimaryAssembly> firstPass = consolidateNearbyAssemblies(results, this::tryMerge);
        return consolidateNearbyAssemblies(firstPass, (left, right) ->
        {
            final Set<String> leftOnly = new HashSet<>(left.getSupportFragments());
            leftOnly.removeAll(right.getSupportFragments());

            final Set<String> rightOnly = new HashSet<>(right.getSupportFragments());
            rightOnly.removeAll(left.getSupportFragments());

            if (leftOnly.size() < 2 && rightOnly.size() < 2)
            {
                final boolean returnRight;
                if(leftOnly.size() == rightOnly.size())
                {
                    if(left.getAverageBaseQuality() == right.getAverageBaseQuality())
                        returnRight = left.getLength() < right.getLength();
                    else
                        returnRight = left.getAverageBaseQuality() < right.getAverageBaseQuality();
                }
                else
                    returnRight = leftOnly.size() < rightOnly.size();

                if(returnRight)
                {
                    right.addErrata(left.getAllErrata());
                    return right;
                }
                else
                {
                    left.addErrata(right.getAllErrata());
                    return left;
                }
            }
            if (leftOnly.size() < 2)
            {
                right.addErrata(left.getAllErrata());
                return right;
            }
            else if (rightOnly.size() < 2)
            {
                left.addErrata(right.getAllErrata());
                return left;
            }
            else
                return null;
        });
    }

    public List<PrimaryAssembly> consolidateNearbyAssemblies(final List<PrimaryAssembly> results,
            final BiFunction<PrimaryAssembly, PrimaryAssembly, PrimaryAssembly> merger)
    {
        results.sort(NaturalSortComparator.<PrimaryAssembly>of(r -> r.AnchorChromosome)
                .thenComparing(r -> r.AnchorPosition));

        final List<PrimaryAssembly> assemblies = new ArrayList<>();
        for(int i = 0; i < results.size(); i++)
        {
            @Nullable
            PrimaryAssembly current = results.get(i);
            if(current == null)
                continue;

            final int maxDedupeDistance = SvConstants.MAXDISTANCETODEDUPEASSEMBLIES;
            final int maxToCheck = Math.min(results.size() - i - 1, maxDedupeDistance * 2);
            for(int j = 0; j < maxToCheck; j++)
            {
                final PrimaryAssembly next = results.get(i + j + 1);
                if(next == null)
                    continue;
                if(!current.AnchorChromosome.equals(next.AnchorChromosome))
                    break;

                final int currentStart = current.AnchorPosition - current.AnchorPositionInAssembly;
                final int nextStart = next.AnchorPosition - next.AnchorPositionInAssembly;

                if(!CommonUtils.overlaps(
                        currentStart - maxDedupeDistance, currentStart + current.Assembly.length() + maxDedupeDistance,
                        nextStart - maxDedupeDistance, nextStart + next.Assembly.length() + maxDedupeDistance))
                {
                    break;
                }

                try
                {
                    @Nullable
                    final PrimaryAssembly merged = merger.apply(current, next);
                    if(merged != null)
                    {
                        current = merged;
                        results.set(i + j + 1, null); // Null out next
                    }
                }
                catch(final Throwable throwable)
                {
                    mContext.Problems.add(new Problem(String.format("Problem merging %s and %s", current.Name, next.Name),
                            throwable, current));
                }
            }

            assemblies.add(current);
        }
        return assemblies;
    }

    @Nullable
    public PrimaryAssembly tryMerge(final PrimaryAssembly left, final PrimaryAssembly right)
    {
        @Nullable
        final Integer mergeIndex = mContext.SupportChecker.AssemblySupport.supportIndex(left, right, 100);
        if(mergeIndex == null)
            return null;

        return merge(left, right, mergeIndex);
    }

    private PrimaryAssembly merge(final PrimaryAssembly left, final PrimaryAssembly right, final int supportIndex)
    {
        final Sequence mergedSequence = SequenceMerger.merge(left, right, supportIndex);

        final var merged = new PrimaryAssembly(left.Name, mergedSequence.getBasesString(),
                "?", 0, 0, left);

        final int leftDelta = supportIndex > 0 ? 0 : -supportIndex;
        for(final Map.Entry<Record, Integer> entry : left.getSupport())
            merged.tryAddSupport(mContext.SupportChecker, entry.getKey(), entry.getValue() + leftDelta);

        final int rightDelta = Math.max(supportIndex, 0);
        for(final Map.Entry<Record, Integer> entry : right.getSupport())
            merged.tryAddSupport(mContext.SupportChecker, entry.getKey(), entry.getValue() + rightDelta);

        merged.addErrata(left.getAllErrata());
        merged.addErrata(right.getAllErrata());

        return merged;
    }

    private List<ExtendedAssembly> order(final Collection<ExtendedAssembly> assemblies)
    {
        // FIXME: Correctly order these
        if(assemblies.size() > 1)
            SV_LOGGER.warn("Found more than 1 assembly ({}) while creating gapped ({})", assemblies.size(),
                    assemblies.stream().map(assembly -> assembly.Name).collect(Collectors.toList()));

        final Map<ExtendedAssembly, Map<ExtendedAssembly, Long>> leftWise = new IdentityHashMap<>();
        final Map<ExtendedAssembly, Map<ExtendedAssembly, Long>> rightWise = new IdentityHashMap<>();
        for (final ExtendedAssembly first : assemblies)
            for (final ExtendedAssembly second : assemblies)
            {
                if (first == second)
                    continue;



            }

        return new ArrayList<>(assemblies);
    }

    public GappedAssembly createGapped(final Collection<ExtendedAssembly> assemblies, final int index)
    {
        final GappedAssembly gappedAssembly = new GappedAssembly("Assembly" + index, order(assemblies));

        for(final ExtendedAssembly assembly : assemblies)
            for(final Record support : assembly.getSupportRecords())
                if(!gappedAssembly.tryAddSupport(mContext.SupportChecker, support))
                    SV_LOGGER.info("Failed to add support for assembly {}: {}", gappedAssembly.Name, support.getName());

        return gappedAssembly;
    }
}
