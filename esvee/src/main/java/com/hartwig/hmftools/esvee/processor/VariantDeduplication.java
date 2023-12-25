package com.hartwig.hmftools.esvee.processor;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.Context;
import com.hartwig.hmftools.esvee.models.Record;
import com.hartwig.hmftools.esvee.util.ParallelMapper;

import org.apache.commons.lang3.tuple.Pair;

public class VariantDeduplication
{
    private final Context mContext;
    private final VariantDeduplicationCounters mCounters;

    public VariantDeduplication(final Context context, final VariantDeduplicationCounters counters)
    {
        mContext = context;
        mCounters = counters;
    }

    public List<VariantCall> deduplicate(final List<VariantCall> variants)
    {
        final List<VariantCall> matchedGroups = deduplicateMatchedGroups(variants);
        final List<VariantCall> cleaned = removeIdenticalAssemblies(matchedGroups);
        final List<VariantCall> withSGLsMerged = mergeSGLsWithNearbyDoubles(cleaned);
        final List<VariantCall> deduped = deduplicateNearbyVariants(withSGLsMerged);

        mCounters.VariantsRemoved.add(variants.size() - deduped.size());

        return deduped;
    }

    /** Finds variants that have identical descriptors, and de-duplicates them */
    private List<VariantCall> deduplicateMatchedGroups(final List<VariantCall> variants)
    {
        final Map<String, List<VariantCall>> callsByPairLocation = new HashMap<>();
        for(final VariantCall call : variants)
        {
            if (call.isSingleSided())
            {
                final String key = Objects.requireNonNullElse(call.LeftDescriptor, call.RightDescriptor);
                callsByPairLocation.computeIfAbsent(key, k -> new ArrayList<>()).add(call);
                continue;
            }

            final var key = String.format("%s:%s-%s:%s (%s/%s)",
                    call.LeftChromosome, call.LeftPosition,
                    call.RightChromosome, call.RightPosition,
                    call.LeftDescriptor, call.RightDescriptor);

            callsByPairLocation.computeIfAbsent(key, k -> new ArrayList<>()).add(call);
        }

        final List<List<VariantCall>> groups = new ArrayList<>(callsByPairLocation.values());

        return ParallelMapper.mapWithCounter(mContext.Executor, groups, this::deduplicateMatchedGroup, null);
    }

    private VariantCall deduplicateMatchedGroup(final List<VariantCall> variants)
    {
        if(variants.size() == 1)
            return variants.get(0);
        mCounters.MatchedGroupRemoved.add(variants.size());

        final Set<Integer> newPhaseSets = new HashSet<>();
        final Set<VariantCall.VariantAssembly> newAssemblies = new HashSet<>();
        final Map<String, List<VariantCall.SampleSupport>> sampleSupport = new LinkedHashMap<>();
        int leftMapQ = 0, rightMapQ = 0;
        for(final VariantCall call : variants)
        {
            newPhaseSets.addAll(call.PhaseSets);
            newAssemblies.addAll(call.variantAssemblies());

            for(final VariantCall.SampleSupport support : call.sampleSupport())
            {
                sampleSupport.computeIfAbsent(support.sampleName(), __ -> new ArrayList<>())
                        .add(support);
            }

            leftMapQ = Math.max(leftMapQ, call.LeftMappingQuality);
            rightMapQ = Math.max(rightMapQ, call.RightMappingQuality);
        }

        final List<VariantCall.SampleSupport> newSampleSupport = new ArrayList<>();
        for(final List<VariantCall.SampleSupport> sampleSupportList : sampleSupport.values())
        {
            final String sampleName = sampleSupportList.get(0).sampleName();
            final boolean isGermline = sampleSupportList.get(0).isGermline();
            final int quality = sampleSupportList.stream().mapToInt(VariantCall.SampleSupport::quality).max().orElseThrow();

            final Set<Record> splitReads = sampleSupportList.stream()
                    .flatMap(s -> s.splitReads().stream())
                    .collect(Collectors.toSet());
            final Set<Record> discordantReads = sampleSupportList.stream()
                    .flatMap(s -> s.discordantReads().stream())
                    .filter(record -> !splitReads.contains(record))
                    .collect(Collectors.toSet());

            newSampleSupport.add(new VariantCall.SampleSupport(sampleName, isGermline, quality, splitReads, discordantReads));
        }

        final VariantCall existing = variants.get(0);
        return VariantCall.create(existing.LeftChromosome, existing.LeftPosition, existing.RightChromosome, existing.RightPosition,
                existing.LeftDescriptor, existing.RightDescriptor, newPhaseSets, newAssemblies,
                leftMapQ, rightMapQ, newSampleSupport, existing.Classification);
    }

    /** For each assembly,  */
    private List<VariantCall> removeIdenticalAssemblies(final List<VariantCall> variants)
    {
        return ParallelMapper.mapWithCounter(mContext.Executor, variants, this::removeIdenticalAssemblies,null);
    }

    private VariantCall removeIdenticalAssemblies(final VariantCall variant)
    {
        if(variant.associatedAssemblies().size() == 1)
            return variant;

        final Map<Pair<String, Integer>, VariantCall.VariantAssembly> assemblies = variant.variantAssemblies().stream()
                .collect(Collectors.toMap(assembly -> Pair.of(assembly.Assembly.Assembly, assembly.LeftPosition),
                        assembly -> assembly, (left, right) ->
                        {
                            if (left.Assembly.getSupportFragments().equals(right.Assembly.getSupportFragments()))
                                return left;

                            // Ensure that both are equally supported
                            for(final Map.Entry<Record, Integer> entry : right.Assembly.getSupport())
                                left.Assembly.addEvidenceAt(entry.getKey(), entry.getValue());

                            return left;
                        }, LinkedHashMap::new));

        final Set<VariantCall.VariantAssembly> newAssemblies = new LinkedHashSet<>(assemblies.values());
        final int assembliesRemoved = variant.associatedAssemblies().size() - newAssemblies.size();
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
        final Map<String, TreeMap<Integer, List<VariantCall>>> indexedVariants = new HashMap<>();
        for (final VariantCall call : variants)
        {
            if (call.isSingleSided())
                continue;

            indexedVariants.computeIfAbsent(call.LeftChromosome, __ -> new TreeMap<>())
                    .computeIfAbsent(call.LeftPosition, __ -> new ArrayList<>())
                    .add(call);
            indexedVariants.computeIfAbsent(call.RightChromosome, __ -> new TreeMap<>())
                    .computeIfAbsent(call.RightPosition, __ -> new ArrayList<>())
                    .add(call);
        }
        for (final VariantCall call : variants)
        {
            if (!call.isSingleSided())
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
