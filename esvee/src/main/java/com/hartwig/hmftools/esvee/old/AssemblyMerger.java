package com.hartwig.hmftools.esvee.old;

import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.BiFunction;

import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.util.CommonUtils;
import com.hartwig.hmftools.esvee.util.NaturalSortComparator;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

public class AssemblyMerger
{
    private final SupportChecker mSupportChecker;

    public AssemblyMerger()
    {
        mSupportChecker = new SupportChecker();
    }

    public List<PrimaryAssembly> consolidatePrimaryAssemblies(final List<PrimaryAssembly> primaryAssemblies) // was consolidateNearbyAssemblies
    {
        final List<PrimaryAssembly> firstPass = consolidateNearbyAssemblies(primaryAssemblies, this::tryMergePrimaryAssemblies);
        return consolidateNearbyAssemblies(firstPass, (left, right) ->
        {
            final Set<String> leftOnly = new HashSet<>(left.getSupportReadNames());
            leftOnly.removeAll(right.getSupportReadNames());

            final Set<String> rightOnly = new HashSet<>(right.getSupportReadNames());
            rightOnly.removeAll(left.getSupportReadNames());

            if(leftOnly.size() < 2 && rightOnly.size() < 2)
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
                    return right;
                else
                    return left;
            }
            if(leftOnly.size() < 2)
            {
                return right;
            }
            else if(rightOnly.size() < 2)
            {
                return left;
            }
            else
                return null;
        });
    }

    private List<PrimaryAssembly> consolidateNearbyAssemblies(
            final List<PrimaryAssembly> results, final BiFunction<PrimaryAssembly, PrimaryAssembly, PrimaryAssembly> merger)
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

            final int maxDedupeDistance = SvConstants.MAX_DISTANCE_EDUPE_ASSEMBLIES;
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
                    SV_LOGGER.error("error consolidating near-by assemblies: {} and {}", current.Name, next.Name);
                }
            }

            assemblies.add(current);
        }
        return assemblies;
    }

    @Nullable
    private PrimaryAssembly tryMergePrimaryAssemblies(final PrimaryAssembly left, final PrimaryAssembly right)
    {
        @Nullable
        final Integer mergeIndex = mSupportChecker.AssemblySupport.supportIndex(left, right, 100);
        if(mergeIndex == null)
            return null;

        return mergePrimaryAssemblies(left, right, mergeIndex);
    }

    private PrimaryAssembly mergePrimaryAssemblies(final PrimaryAssembly left, final PrimaryAssembly right, final int supportIndex)
    {
        final Sequence mergedSequence = SequenceMerger.merge(left, right, supportIndex);

        final var merged = new PrimaryAssembly(left.Name, mergedSequence.getBasesString(), left.OriginalJunction,
                "?", 0, 0, left);

        final int leftDelta = supportIndex > 0 ? 0 : -supportIndex;

        for(ReadSupport readSupport : left.readSupport())
        {
            merged.tryAddSupport(mSupportChecker, readSupport.Read, readSupport.Index + leftDelta);
        }

        final int rightDelta = Math.max(supportIndex, 0);

        for(ReadSupport readSupport : right.readSupport())
        {
            merged.tryAddSupport(mSupportChecker, readSupport.Read, readSupport.Index + rightDelta);
        }

        return merged;
    }

    public Set<ExtendedAssembly> primaryPhasedMerging(final Set<ExtendedAssembly> primaryPhaseSet)
    {
        // CHECK: use of break goto
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

                            mergedAssembly = mergeExtendedAssembly(left, flippedRight, index);
                        }
                        else
                            mergedAssembly = mergeExtendedAssembly(left, right, index);

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

            /* FIXME:
            SV_LOGGER.warn("{}", RegionOfInterest.tryMerge(
                    primaryPhaseSet.stream()
                            .flatMap(assembly -> assembly.getSupport().stream())
                            .map(Map.Entry::getKey)
                            .filter(record -> !record.isUnmapped())
                            .map(record -> new RegionOfInterest(record.getChromosome(), record.getAlignmentStart(), record.getAlignmentEnd()))
                            .collect(Collectors.toList())
            ));

             */
            return null;
        }
    }

    private ExtendedAssembly mergeExtendedAssembly(final ExtendedAssembly left, final ExtendedAssembly right, final int supportIndex)
    {
        left.markDecompositionStale();
        right.markDecompositionStale();
        final Sequence mergedSequence = SequenceMerger.merge(left, right, supportIndex);

        final ExtendedAssembly merged = new ExtendedAssembly(left.Name, mergedSequence.getBasesString(), left.Source);

        left.readSupport().forEach(x -> merged.tryAddSupport(mSupportChecker, x.Read));
        right.readSupport().forEach(x -> merged.tryAddSupport(mSupportChecker, x.Read));

        return merged;
    }
}
