package com.hartwig.hmftools.esvee.assembly;

import static java.lang.String.format;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.common.SampleSupport;
import com.hartwig.hmftools.esvee.common.VariantAssembly;
import com.hartwig.hmftools.esvee.common.VariantCall;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.sequence.ReadSupport;

import org.apache.commons.lang3.tuple.Pair;

public class VariantDeduplication
{
    private final VariantDeduplicationCounters mCounters;

    public VariantDeduplication()
    {
        mCounters = new VariantDeduplicationCounters(); // previously took a global set
    }

    public List<VariantCall> deduplicate(final List<VariantCall> variants)
    {
        List<VariantCall> matchedGroups = deduplicateMatchedGroups(variants);
        List<VariantCall> cleaned = removeIdenticalAssemblies(matchedGroups);
        List<VariantCall> withSGLsMerged = mergeSGLsWithNearbyDoubles(cleaned);
        List<VariantCall> deduped = deduplicateNearbyVariants(withSGLsMerged);

        mCounters.VariantsRemoved.add(variants.size() - deduped.size());

        return deduped;
    }

    private List<VariantCall> deduplicateMatchedGroups(final List<VariantCall> variants)
    {
        // finds variants that have identical descriptors, and de-duplicate them
        Map<String,List<VariantCall>> variantsByPairLocation = Maps.newHashMap();

        for(VariantCall variant : variants)
        {
            String variantKey = formVariantKey(variant);
            variantsByPairLocation.computeIfAbsent(variantKey, k -> new ArrayList<>()).add(variant);
        }

        List<List<VariantCall>> groups = new ArrayList<>(variantsByPairLocation.values());

        // FIXME: was multi-threaded - is it beneficial to be?
        // return ParallelMapper.mapWithCounter(mContext.Executor, groups, this::deduplicateMatchedGroup, null);

        List<VariantCall> dedupedVariants = groups.stream().map(x -> deduplicateMatchedGroup(x)).collect(Collectors.toList());
        return dedupedVariants;
    }

    private static String formVariantKey(final VariantCall variant)
    {
        if(variant.isSingleSided())
        {
            return variant.LeftDescriptor != null ? variant.LeftDescriptor : variant.RightDescriptor;
        }
        else
        {
            return format("%s:%s-%s:%s (%s/%s)",
                    variant.LeftChromosome, variant.LeftPosition, variant.RightChromosome, variant.RightPosition,
                    variant.LeftDescriptor, variant.RightDescriptor);
        }
    }

    private VariantCall deduplicateMatchedGroup(List<VariantCall> variants)
    {
        if(variants.size() == 1)
            return variants.get(0);

        mCounters.MatchedGroupRemoved.add(variants.size());

        Set<Integer> newPhaseSets = new HashSet<>();
        Set<VariantAssembly> newAssemblies = new HashSet<>();
        Map<String, List<SampleSupport>> sampleSupport = new LinkedHashMap<>();
        int leftMapQ = 0, rightMapQ = 0;
        for(VariantCall call : variants)
        {
            newPhaseSets.addAll(call.PhaseSets);
            newAssemblies.addAll(call.variantAssemblies());

            for(SampleSupport support : call.sampleSupport())
            {
                sampleSupport.computeIfAbsent(support.sampleName(), __ -> new ArrayList<>())
                        .add(support);
            }

            leftMapQ = Math.max(leftMapQ, call.LeftMappingQuality);
            rightMapQ = Math.max(rightMapQ, call.RightMappingQuality);
        }

        List<SampleSupport> newSampleSupport = new ArrayList<>();
        for(List<SampleSupport> sampleSupportList : sampleSupport.values())
        {
            String sampleName = sampleSupportList.get(0).sampleName();
            int quality = sampleSupportList.stream().mapToInt(SampleSupport::quality).max().orElseThrow();

            Set<Read> splitReads = sampleSupportList.stream()
                    .flatMap(s -> s.splitReads().stream())
                    .collect(Collectors.toSet());
            Set<Read> discordantReads = sampleSupportList.stream()
                    .flatMap(s -> s.discordantReads().stream())
                    .filter(record -> !splitReads.contains(record))
                    .collect(Collectors.toSet());

            newSampleSupport.add(new SampleSupport(sampleName, quality, splitReads, discordantReads));
        }

        VariantCall existing = variants.get(0);
        return VariantCall.create(existing.LeftChromosome, existing.LeftPosition, existing.RightChromosome, existing.RightPosition,
                existing.LeftDescriptor, existing.RightDescriptor, newPhaseSets, newAssemblies,
                leftMapQ, rightMapQ, newSampleSupport, existing.Classification);
    }

    private List<VariantCall> removeIdenticalAssemblies(final List<VariantCall> variants)
    {
        // FIXME: was multi-threaded - how beneficial
        // return ParallelMapper.mapWithCounter(mContext.Executor, variants, this::removeIdenticalAssemblies,null);

        List<VariantCall> adjustedVariants = variants.stream().map(x -> removeIdenticalAssemblies(x)).collect(Collectors.toList());
        return adjustedVariants;
    }

    private VariantCall removeIdenticalAssemblies(final VariantCall variant)
    {
        if(variant.associatedAssemblies().size() == 1)
            return variant;

        Map<Pair<String, Integer>, VariantAssembly> assemblies = variant.variantAssemblies().stream()
                .collect(Collectors.toMap(assembly -> Pair.of(assembly.Assembly.Assembly, assembly.LeftPosition),
                        assembly -> assembly, (left, right) ->
                        {
                            if(left.Assembly.getSupportReadNames().equals(right.Assembly.getSupportReadNames()))
                                return left;

                            // Ensure that both are equally supported
                            for(ReadSupport readSupport : right.Assembly.readSupport())
                            {
                                left.Assembly.addEvidenceAt(readSupport.Read, readSupport.Index);
                            }

                            return left;
                        }, LinkedHashMap::new));

        Set<VariantAssembly> newAssemblies = new LinkedHashSet<>(assemblies.values());
        int assembliesRemoved = variant.associatedAssemblies().size() - newAssemblies.size();

        if(assembliesRemoved == 0)
            return variant;

        mCounters.AssembliesRemoved.add(assembliesRemoved);

        return VariantCall.create(variant.LeftChromosome, variant.LeftPosition,
                variant.RightChromosome, variant.RightPosition,
                variant.LeftDescriptor, variant.RightDescriptor, variant.PhaseSets,
                newAssemblies, variant.LeftMappingQuality, variant.RightMappingQuality,
                variant.sampleSupport(), variant.Classification);
    }

    private List<VariantCall> mergeSGLsWithNearbyDoubles(final List<VariantCall> variants)
    {
        // Index double-ended variants
        Map<String, TreeMap<Integer, List<VariantCall>>> indexedVariants = new HashMap<>();
        for(VariantCall call : variants)
        {
            if(call.isSingleSided())
                continue;

            indexedVariants.computeIfAbsent(call.LeftChromosome, __ -> new TreeMap<>())
                    .computeIfAbsent(call.LeftPosition, __ -> new ArrayList<>())
                    .add(call);
            indexedVariants.computeIfAbsent(call.RightChromosome, __ -> new TreeMap<>())
                    .computeIfAbsent(call.RightPosition, __ -> new ArrayList<>())
                    .add(call);
        }
        for(VariantCall call : variants)
        {
            if(!call.isSingleSided())
                continue;


        }

        // FIXME: This
        return variants;
    }

    private List<VariantCall> deduplicateNearbyVariants(final List<VariantCall> variants)
    {
        // FIXME: This
        return variants;
    }
}
