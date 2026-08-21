package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_SHIFT_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.REGION_UNCOVERED_MAX;
import static com.hartwig.hmftools.panelbuilder.ProbeSelector.selectBestProbe;
import static com.hartwig.hmftools.panelbuilder.ProbeTiling.calculateProbeTiling;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.maxProbeEndOverlapping;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.maxProbeEndWithoutGap;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.minProbeStartOverlapping;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.minProbeStartWithoutGap;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionCenteredAt;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionEndingAt;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionStartingAt;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionCentre;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionCentreStartOffset;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionIntersection;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionNegatedIntersection;
import static com.hartwig.hmftools.panelbuilder.Utils.outwardMovingOffsets;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.function.Consumer;
import java.util.function.IntPredicate;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.panelbuilder.probequality.ProbeQualityModel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

// Encapsulates all functionality for creating probes.
public class ProbeGenerator
{
    private final Map<String, Integer> mChromosomeLengths;
    private final ProbeEvaluator mProbeEvaluator;
    // Hook to catch all candidate probes for output.
    @Nullable
    private final Consumer<Probe> mCandidateCallback;

    private static final Logger LOGGER = LogManager.getLogger(ProbeGenerator.class);

    public ProbeGenerator(final Map<String, Integer> chromosomeLengths, final ProbeEvaluator probeEvaluator,
            final @Nullable Consumer<Probe> candidateCallback)
    {
        mChromosomeLengths = chromosomeLengths;
        mProbeEvaluator = probeEvaluator;
        mCandidateCallback = candidateCallback;
    }

    public static ProbeGenerator construct(final RefGenomeInterface refGenome, final ProbeQualityProfile probeQualityProfile,
            final ProbeQualityModel probeQualityModel, final Consumer<Probe> candidateCallback)
    {
        ProbeQualityScorer probeQualityScorer = new ProbeQualityScorer(probeQualityProfile, probeQualityModel);
        ProbeEvaluator probeEvaluator = new ProbeEvaluator(refGenome, probeQualityScorer);
        return new ProbeGenerator(refGenome.chromosomeLengths(), probeEvaluator, candidateCallback);
    }

    // Batches probe generation so the expensive path - scoring a novel (non-reference-contiguous) sequence via BWA alignment - can align many
    // candidates in one call, where BWA parallelises across queries (~7x speedup). Only spec types that produce novel sequences are pooled:
    // variant/SV probes span a breakend so are pooled across all specs; RNA exon ranges cross splice junctions so are pooled per exon (which
    // already saturates the aligner). Reference-contiguous probes are profile-scored with nothing to parallelise, so they run one at a time.
    public class Batch
    {
        private final List<ProbeGenerationSpec> mGenericSpecs = new ArrayList<>();
        private final List<ProbeGenerationSpec.SingleProbe> mSingleProbeSpecs = new ArrayList<>();
        private final List<ProbeGenerationSpec.CoverExonRange> mExonRangeSpecs = new ArrayList<>();

        public void add(final ProbeGenerationSpec spec)
        {
            if(spec instanceof ProbeGenerationSpec.SingleProbe singleProbeSpec)
            {
                mSingleProbeSpecs.add(singleProbeSpec);
            }
            else if(spec instanceof ProbeGenerationSpec.CoverExonRange exonRangeSpec)
            {
                mExonRangeSpecs.add(exonRangeSpec);
            }
            else
            {
                mGenericSpecs.add(spec);
            }
        }

        public void addAll(Stream<ProbeGenerationSpec> specs)
        {
            specs.forEachOrdered(this::add);
        }

        public ProbeGenerationResult generate(PanelBuffer resultBuffer)
        {
            ProbeGenerationResult result = generateGenericSpecs(resultBuffer, resultBuffer)
                    .add(generateSingleProbeSpecs(resultBuffer, resultBuffer))
                    .add(generateExonRangeSpecs(mExonRangeSpecs, resultBuffer));
            clear();
            return result;
        }

        private ProbeGenerationResult generateGenericSpecs(final PanelCoverage coverage, PanelStore resultStore)
        {
            return mGenericSpecs.stream().map(spec ->
                    {
                        ProbeGenerationResult result = generateGenericSpec(spec, coverage);
                        resultStore.addResult(result);
                        return result;
                    })
                    .reduce(new ProbeGenerationResult(), ProbeGenerationResult::add);
        }

        private ProbeGenerationResult generateSingleProbeSpecs(final PanelCoverage coverage, PanelStore resultStore)
        {
            return singleProbeBatch(mSingleProbeSpecs, coverage, resultStore);
        }

        private void clear()
        {
            mGenericSpecs.clear();
            mSingleProbeSpecs.clear();
            mExonRangeSpecs.clear();
        }
    }

    public ProbeGenerationResult generateBatch(Stream<ProbeGenerationSpec> specs, PanelBuffer resultBuffer)
    {
        Batch batch = new Batch();
        batch.addAll(specs);
        return batch.generate(resultBuffer);
    }

    private ProbeGenerationResult generateGenericSpec(final ProbeGenerationSpec spec, final PanelCoverage coverage)
    {
        if(spec instanceof ProbeGenerationSpec.CoverRegion specObj)
        {
            return coverRegion(
                    specObj.region(), specObj.metadata(), specObj.evalCriteria(), specObj.localSelectStrategy(), coverage);
        }
        else if(spec instanceof ProbeGenerationSpec.CoverOneSubregion specObj)
        {
            return coverOneSubregion(
                    specObj.region(), specObj.metadata(), specObj.evalCriteria(), specObj.selectStrategy(), coverage);
        }
        else if(spec instanceof ProbeGenerationSpec.CoverOnePosition specObj)
        {
            return coverOnePosition(
                    specObj.positions(), specObj.metadata(), specObj.evalCriteria(), specObj.selectStrategy(), coverage);
        }
        else
        {
            throw new IllegalArgumentException("Unhandled ProbeGenerationSpec");
        }
    }

    // General purpose method for generating the best acceptable probes to cover an entire region.
    // Behaviour:
    //   - Prefer a probe set which is symmetrical and centered on the target region, unless this would violate constraints.
    //   - Avoid placing probes in already covered regions, if `coverage` is not null. But not guaranteed!
    //   - The edges of the region may be slightly uncovered (since the probes will capture a slightly larger region).
    //   - Probes may overlap and extend outside the target region.
    //   - Probes are shifted slightly to optimise for the selection criteria.
    private ProbeGenerationResult coverRegion(final ChrBaseRegion region, final TargetMetadata metadata,
            final ProbeEvaluator.Criteria evalCriteria, final ProbeSelector.Strategy localSelect, final PanelCoverage coverage)
    {
        List<ChrBaseRegion> subregions;
        // Split the region into uncovered subregions to avoid overlap with regions already covered by probes.
        subregions = regionNegatedIntersection(
                region.baseRegion(), coverage.coveredRegions(region.chromosome())
                        .map(ChrBaseRegion::baseRegion))
                .stream()
                .map(baseRegion -> ChrBaseRegion.from(region.chromosome(), baseRegion))
                .toList();
        if(subregions.size() > 1)
        {
            subregions.forEach(subregion -> LOGGER.debug("Split region into uncovered subregion: {}", subregion));
        }

        // Identity (whole-chromosome) mapping with no pinned edges, so this reproduces the plain genome-space tiling exactly.
        RegionMapping mapping = RegionMapping.wholeChromosome(region.chromosome(), mChromosomeLengths.get(region.chromosome()));
        IntPredicate noPinning = position -> false;
        ProbeGenerationResult result = subregions.stream()
                .map(subregion -> coverMappedRange(
                        mapping, subregion.start() - 1, subregion.end(), metadata, evalCriteria, localSelect, noPinning))
                .reduce(new ProbeGenerationResult(), ProbeGenerationResult::add);

        // Add in the candidate target region that's not added by coverSubregion().
        TargetRegion candidateTargetRegion = new TargetRegion(region, metadata);
        result = result.add(new ProbeGenerationResult(emptyList(), List.of(candidateTargetRegion), emptyList()));

        return result;
    }

    private record CoverAcceptableSubregionResult(
            List<Probe> probes,
            // Region where the tiling algorithm decided it was optimal to not place probes, but shouldn't be marked as rejected.
            List<BaseRegion> permittedUncoveredRegions
    )
    {
    }

    // Generates exon-aware (RNA) probes for a batch of exon-range specs, combining the results. Tiling edges are pinned flush to exon
    // boundaries for splice-junction coverage; a window near an edge becomes a spliced multi-region probe.
    private ProbeGenerationResult generateExonRangeSpecs(final List<ProbeGenerationSpec.CoverExonRange> specs, final PanelStore resultStore)
    {
        ProbeGenerationResult total = new ProbeGenerationResult();
        for(ProbeGenerationSpec.CoverExonRange spec : specs)
        {
            ProbeGenerationResult result = coverMappedRange(
                    spec.mapping(), spec.rangeStart(), spec.rangeEnd(), spec.metadata(), spec.evalCriteria(), spec.localSelectStrategy(),
                    spec.mapping()::isRegionBoundary);
            // The coverage step does not add the candidate target region; add the exon region here.
            ChrBaseRegion exonRegion = spec.mapping().toGenomeRegions(spec.rangeStart(), spec.rangeEnd()).get(0);
            result = result.add(
                    new ProbeGenerationResult(emptyList(), List.of(new TargetRegion(exonRegion, spec.metadata())), emptyList()));
            resultStore.addResult(result);
            total = total.add(result);
        }
        return total;
    }

    // Covers a probe-space target range [rangeStart, rangeEnd) using the region mapping (identity for DNA, ordered exons for RNA). Candidate
    // windows map back to genome regions - one for DNA, possibly spliced for RNA. Probe-space positions are wrapped as 1-indexed regions so
    // the shared region maths apply. The caller adds candidate target regions, not this method.
    private ProbeGenerationResult coverMappedRange(final RegionMapping mapping, int rangeStart, int rangeEnd, final TargetMetadata metadata,
            final ProbeEvaluator.Criteria evalCriteria, final ProbeSelector.Strategy localSelect, final IntPredicate pinBoundary)
    {
        int minStart = max(0, rangeStart - PROBE_LENGTH + 1);
        int maxStart = min(mapping.length() - PROBE_LENGTH, rangeEnd - 1);
        BaseRegion targetRegion = spaceRegion(rangeStart, rangeEnd);
        // TODO? When the whole mapping is shorter than a probe, no window fits and the target gets no coverage - the short-exon padding limit.
        //  Probes are fixed length throughout, so a fix means a pipeline-wide change (e.g. shorter probes). Accepted for now; only tiny
        //  non-coding RNAs are this short, and the uncovered target is reported in the rejection output.

        // Accepted candidate windows indexed by probe-space start (as a 1-indexed coordinate), merged into maximal acceptable sub-ranges.
        Map<Integer, Probe> acceptableProbes = new HashMap<>();

        // Pass 1: cheap within-exon (single-region) candidates only. The expensive spliced candidates are deferred to pass 2 and built only
        // where these leave the target uncovered. For the DNA identity mapping every candidate is single-region, so pass 2 is inert.
        List<Probe> singleRegionCandidates = IntStream.rangeClosed(minStart, maxStart)
                .mapToObj(start -> probeAtStart(mapping, start, rangeStart, rangeEnd, metadata))
                .filter(probe -> probe.definition().regions().size() == 1)
                .toList();
        evaluateProbes(singleRegionCandidates.stream(), evalCriteria)
                .filter(Probe::accepted)
                .forEachOrdered(probe -> acceptableProbes.put(spaceWindow(probeSpaceStart(mapping, probe)).start(), probe));

        // Pass 2: spliced candidates only where no acceptable single-region probe covers the target (short exons, or sub-ranges abutting a
        // junction). Bounded to the junction neighbourhood, so it is empty for a covered exon interior.
        List<BaseRegion> uncoveredBySingleRegion = regionNegatedIntersection(
                targetRegion, acceptableProbes.values().stream().map(probe -> spaceWindow(probeSpaceStart(mapping, probe))));
        if(!uncoveredBySingleRegion.isEmpty())
        {
            SortedSet<Integer> splicedStarts = allSplicedStartsOverlapping(uncoveredBySingleRegion, minStart, maxStart);

            // Short-exon fast path: with no single-region candidate the exon is covered by one centred padded probe placeable only near the
            // centred position, so evaluate just that neighbourhood, falling back to the full sweep only if nothing lands there. Exact.
            SortedSet<Integer> toEvaluate = splicedStarts;
            List<Integer> centreStarts = List.of();
            if(singleRegionCandidates.isEmpty())
            {
                BaseRegion tilingBounds = new BaseRegion(minStart + 1, maxStart + PROBE_LENGTH);
                List<Integer> tiling = calculateProbeTiling(targetRegion, tilingBounds, false, false);
                if(!tiling.isEmpty())
                {
                    centreStarts = tiling.stream().map(probeStart -> probeStart - 1).toList();
                    SortedSet<Integer> nearCentre = new TreeSet<>();
                    for(int centreStart : centreStarts)
                    {
                        for(int start = max(minStart, centreStart - PROBE_SHIFT_MAX);
                                start <= min(maxStart, centreStart + PROBE_SHIFT_MAX); ++start)
                        {
                            nearCentre.add(start);
                        }
                    }
                    toEvaluate = nearCentre;
                }
            }

            evaluateSplicedInto(mapping, toEvaluate, rangeStart, rangeEnd, metadata, evalCriteria, acceptableProbes);

            // If the fast path placed no probe at a centred position, the real placement may shift; evaluate the rest to stay exact.
            boolean centredPlaced = centreStarts.stream().anyMatch(start -> acceptableProbes.containsKey(start + 1));
            if(!centreStarts.isEmpty() && !centredPlaced)
            {
                splicedStarts.removeAll(toEvaluate);
                evaluateSplicedInto(mapping, splicedStarts, rangeStart, rangeEnd, metadata, evalCriteria, acceptableProbes);
            }
        }

        // Merge all accepted candidate windows (ascending) into maximal acceptable sub-ranges.
        List<BaseRegion> acceptableSubregions = new ArrayList<>();
        acceptableProbes.keySet().stream().sorted().forEach(windowStart ->
        {
            BaseRegion window = new BaseRegion(windowStart, windowStart + PROBE_LENGTH - 1);
            BaseRegion prev = acceptableSubregions.isEmpty() ? null : acceptableSubregions.get(acceptableSubregions.size() - 1);
            if(prev != null && window.start() <= prev.end() + 1)
            {
                if(window.end() > prev.end())
                {
                    acceptableSubregions.set(acceptableSubregions.size() - 1, new BaseRegion(prev.start(), window.end()));
                }
            }
            else
            {
                acceptableSubregions.add(window);
            }
        });

        List<Probe> probes = new ArrayList<>();
        List<BaseRegion> permittedUncoveredRegions = new ArrayList<>();
        for(BaseRegion acceptableSubregion : acceptableSubregions)
        {
            CoverAcceptableSubregionResult result =
                    coverAcceptableMappedSubrange(mapping, acceptableSubregion, targetRegion, acceptableProbes, localSelect, pinBoundary);
            probes.addAll(result.probes());
            permittedUncoveredRegions.addAll(result.permittedUncoveredRegions());
        }

        Stream<BaseRegion> probeWindows = probes.stream().map(probe -> spaceWindow(probeSpaceStart(mapping, probe)));
        Stream<BaseRegion> unrejectedRegions = Stream.concat(probeWindows, permittedUncoveredRegions.stream());
        List<RejectedFeature> rejectedFeatures = regionNegatedIntersection(targetRegion, unrejectedRegions).stream()
                .flatMap(uncovered -> mapping.toGenomeRegions(uncovered.start() - 1, uncovered.end()).stream())
                .map(region -> RejectedFeature.fromRegion(region, metadata))
                .toList();

        return new ProbeGenerationResult(probes, emptyList(), rejectedFeatures);
    }

    // Tiles and shifts probes over an acceptable sub-range. A tiling edge coinciding with an exon boundary is pinned flush; otherwise the
    // usual shift constraints apply (probes may extend past the target). Pinning never triggers for the DNA identity mapping.
    private CoverAcceptableSubregionResult coverAcceptableMappedSubrange(final RegionMapping mapping, final BaseRegion acceptableSubregion,
            final BaseRegion targetRegion, final Map<Integer, Probe> acceptableProbes, final ProbeSelector.Strategy localSelect,
            final IntPredicate pinBoundary)
    {
        BaseRegion tilingTarget = regionIntersection(acceptableSubregion, targetRegion).orElseThrow();
        BaseRegion probeBounds = acceptableSubregion;

        // Pinning the outermost probes flush to boundaries is only meaningful with multiple probes; a range fitting a single probe is centred
        // instead. Inert for DNA.
        boolean singleProbe = tilingTarget.baseLength() <= PROBE_LENGTH + REGION_UNCOVERED_MAX;
        boolean pinStart = !singleProbe && pinBoundary.test(tilingTarget.start() - 1);
        boolean pinEnd = !singleProbe && pinBoundary.test(tilingTarget.end());

        List<BaseRegion> tiling = calculateProbeTiling(tilingTarget, probeBounds, pinStart, pinEnd).stream()
                .map(ProbeUtils::probeRegionStartingAt)
                .toList();

        List<Probe> finalProbes = new ArrayList<>();
        List<Integer> finalStarts = new ArrayList<>();
        boolean[] couldPlaceProbe = new boolean[tiling.size()];
        for(int i = 0; i < tiling.size(); ++i)
        {
            BaseRegion originalProbe = tiling.get(i);
            BaseRegion prevProbe = finalStarts.isEmpty() ? null : probeRegionStartingAt(finalStarts.get(finalStarts.size() - 1));
            BaseRegion nextProbe = i + 1 < tiling.size() ? tiling.get(i + 1) : null;

            int minStart = max(probeBounds.start(), minProbeStartOverlapping(targetRegion));
            int maxEnd = min(probeBounds.end(), maxProbeEndOverlapping(targetRegion));
            if(prevProbe != null)
            {
                maxEnd = min(maxEnd, maxProbeEndWithoutGap(prevProbe));
            }
            if(nextProbe != null)
            {
                minStart = max(minStart, minProbeStartWithoutGap(nextProbe));
            }

            if(i == 0)
            {
                if(originalProbe.start() < targetRegion.start())
                {
                    minStart = originalProbe.start();
                    maxEnd = min(maxEnd, probeRegionStartingAt(targetRegion.start()).end());
                }
                else
                {
                    minStart = max(minStart, targetRegion.start());
                    maxEnd = originalProbe.end();
                }
            }
            if(i == tiling.size() - 1)
            {
                if(originalProbe.end() > targetRegion.end())
                {
                    minStart = max(minStart, probeRegionEndingAt(targetRegion.end()).start());
                    maxEnd = originalProbe.end();
                }
                else
                {
                    minStart = originalProbe.start();
                    maxEnd = min(maxEnd, targetRegion.end());
                }
            }

            // A pinned outer edge is fixed flush to the boundary (no shift).
            boolean pinnedEdge = (i == 0 && pinStart) || (i == tiling.size() - 1 && pinEnd);
            if(pinnedEdge)
            {
                minStart = originalProbe.start();
                maxEnd = originalProbe.end();
            }
            else
            {
                minStart = max(minStart, originalProbe.start() - PROBE_SHIFT_MAX);
                maxEnd = min(maxEnd, originalProbe.end() + PROBE_SHIFT_MAX);
            }

            int maxStart = probeRegionEndingAt(maxEnd).start();
            if(minStart > maxStart)
            {
                minStart = maxStart = originalProbe.start();
            }

            int minOffset = minStart - originalProbe.start();
            int maxOffset = maxStart - originalProbe.start();
            Stream<Probe> candidates = outwardMovingOffsets(minOffset, maxOffset)
                    .map(offset -> originalProbe.start() + offset)
                    .mapToObj(acceptableProbes::get)
                    .filter(Objects::nonNull);
            Optional<Probe> bestProbe = selectBestProbe(candidates, localSelect);
            if(bestProbe.isPresent())
            {
                finalProbes.add(bestProbe.get());
                finalStarts.add(spaceWindow(probeSpaceStart(mapping, bestProbe.get())).start());
            }
            couldPlaceProbe[i] = bestProbe.isPresent();
        }

        List<BaseRegion> permittedUncoveredRegions = new ArrayList<>();
        if(!tiling.isEmpty() && couldPlaceProbe[0])
        {
            int tilingStart = finalStarts.get(0);
            if(tilingStart > acceptableSubregion.start())
            {
                permittedUncoveredRegions.add(new BaseRegion(acceptableSubregion.start(), tilingStart));
            }
        }
        if(!tiling.isEmpty() && couldPlaceProbe[tiling.size() - 1])
        {
            int tilingEnd = finalStarts.get(finalStarts.size() - 1) + PROBE_LENGTH - 1;
            if(tilingEnd < acceptableSubregion.end())
            {
                permittedUncoveredRegions.add(new BaseRegion(tilingEnd, acceptableSubregion.end()));
            }
        }

        return new CoverAcceptableSubregionResult(finalProbes, permittedUncoveredRegions);
    }

    // Builds the candidate probe whose probe-length window starts at the given probe-space position, mapping the window to genome regions (a
    // single region within one exon, or spliced across exon boundaries).
    private static Probe probeAtStart(final RegionMapping mapping, int start, int rangeStart, int rangeEnd, final TargetMetadata metadata)
    {
        List<ChrBaseRegion> regions = mapping.toGenomeRegions(start, start + PROBE_LENGTH);
        SequenceDefinition definition = SequenceDefinition.spliced(regions);
        TargetedRange targetedRange = new TargetedRange(max(0, rangeStart - start), min(PROBE_LENGTH, rangeEnd - start));
        return new Probe(definition, targetedRange, metadata);
    }

    // Probe-space starts whose probe-length window overlaps any of the given (1-indexed) ranges, within the mapping start bounds.
    private static SortedSet<Integer> allSplicedStartsOverlapping(final List<BaseRegion> ranges, int minStart, int maxStart)
    {
        SortedSet<Integer> starts = new TreeSet<>();
        for(BaseRegion range : ranges)
        {
            int lo = max(minStart, range.start() - PROBE_LENGTH);
            int hi = min(maxStart, range.end() - 1);
            for(int start = lo; start <= hi; ++start)
            {
                starts.add(start);
            }
        }
        return starts;
    }

    // Evaluates the spliced (multi-region) candidate at each given probe-space start and adds the accepted ones to the acceptable-probes map
    // (keyed by 1-indexed window start). Single-region starts are ignored (handled by pass 1).
    private void evaluateSplicedInto(final RegionMapping mapping, final Collection<Integer> starts, int rangeStart, int rangeEnd,
            final TargetMetadata metadata, final ProbeEvaluator.Criteria evalCriteria, final Map<Integer, Probe> acceptableProbes)
    {
        Stream<Probe> candidates = starts.stream()
                .map(start -> probeAtStart(mapping, start, rangeStart, rangeEnd, metadata))
                .filter(probe -> probe.definition().regions().size() > 1);
        evaluateProbes(candidates, evalCriteria)
                .filter(Probe::accepted)
                .forEachOrdered(probe -> acceptableProbes.put(spaceWindow(probeSpaceStart(mapping, probe)).start(), probe));
    }

    private static int probeSpaceStart(final RegionMapping mapping, final Probe probe)
    {
        ChrBaseRegion first = probe.definition().regions().get(0);
        return mapping.toProbeSpacePosition(first.chromosome(), first.start()).orElseThrow();
    }

    // A probe-length window at the given probe-space start, as a 1-indexed BaseRegion.
    private static BaseRegion spaceWindow(int spaceStart)
    {
        return new BaseRegion(spaceStart + 1, spaceStart + PROBE_LENGTH);
    }

    // A probe-space range [spaceStart, spaceEnd) as a 1-indexed inclusive BaseRegion.
    private static BaseRegion spaceRegion(int spaceStart, int spaceEnd)
    {
        return new BaseRegion(spaceStart + 1, spaceEnd);
    }

    // Generates all probes overlapping a region, in order from left to right.
    // Use with care. Generally, you would want to choose probes with a more careful approach.
    protected Stream<Probe> allOverlappingProbes(final ChrBaseRegion region, final TargetMetadata metadata)
    {
        int minProbeStart = max(minProbeStartOverlapping(region.baseRegion()), 1);
        int maxProbeEnd = min(maxProbeEndOverlapping(region.baseRegion()), mChromosomeLengths.get(region.chromosome()));
        int maxProbeStart = probeRegionEndingAt(maxProbeEnd).start();
        return IntStream.rangeClosed(minProbeStart, maxProbeStart)
                .mapToObj(start ->
                {
                    ChrBaseRegion probeRegion = probeRegionStartingAt(region.chromosome(), start);
                    SequenceDefinition definition = SequenceDefinition.singleRegion(probeRegion);
                    TargetedRange targetedRange = TargetedRange.fromRegions(region, probeRegion);
                    return new Probe(definition, targetedRange, metadata);
                });
    }

    // Generates the one best acceptable probe that is contained within the specified region.
    private ProbeGenerationResult coverOneSubregion(final ChrBaseRegion region, final TargetMetadata metadata,
            final ProbeEvaluator.Criteria evalCriteria, final ProbeSelector.Strategy selectStrategy, final PanelCoverage coverage)
    {
        TargetRegion candidateTargetRegion = new TargetRegion(region, metadata);
        Stream<Probe> candidates = coverOneSubregionCandidates(region, metadata);
        return selectBestCandidate(candidates, evalCriteria, selectStrategy)
                .map(probe ->
                {
                    if(coverage.isCovered(probe.definition(), probe.targetedRange()))
                    {
                        return ProbeGenerationResult.alreadyCoveredTargets(List.of(candidateTargetRegion));
                    }
                    else
                    {
                        return new ProbeGenerationResult(List.of(probe), List.of(candidateTargetRegion), emptyList());
                    }
                })
                .orElseGet(() -> ProbeGenerationResult.rejectTargets(List.of(candidateTargetRegion)));
    }

    // Generates candidate probes covering any subregion within the target region.
    // Probes do not extend outside the target region.
    protected Stream<Probe> coverOneSubregionCandidates(final ChrBaseRegion region, final TargetMetadata metadata)
    {
        if(region.baseLength() < PROBE_LENGTH)
        {
            // This method is designed to find probes contained within the region, which requires the region fit at least one probe.
            throw new IllegalArgumentException("region must be larger than a probe");
        }

        BasePosition initialPosition = regionCentre(region);
        // Stop once the probes go outside the target region.
        int minProbeStart = region.start();
        int maxProbeEnd = region.end();
        return outwardMovingCenterAlignedProbes(initialPosition, minProbeStart, maxProbeEnd)
                .map(probeRegion ->
                {
                    SequenceDefinition definition = SequenceDefinition.singleRegion(probeRegion);
                    TargetedRange targetedRange = TargetedRange.fromRegions(region, probeRegion);
                    return new Probe(definition, targetedRange, metadata);
                });
    }

    // Returns candidate probe regions shifting left and right with offsets: 0, 1, -1, 2, -2, 3, -3, ...
    // Probes are aligned to their centre position.
    // Probe bounds:
    //   - Can't start before start of chromosome
    //   - Can't start before `minProbeStart`
    //   - Can't end after `maxProbeEnd`
    //   - Can't end after end of chromosome
    // Useful because we prefer to select probes which are closest to the target position or centre of a region.
    private Stream<ChrBaseRegion> outwardMovingCenterAlignedProbes(final BasePosition initialPosition, int minProbeStart, int maxProbeEnd)
    {
        if(maxProbeEnd - minProbeStart + 1 < PROBE_LENGTH)
        {
            // Probably indicates a bug in the calling code.
            throw new IllegalArgumentException("minProbeStart and maxProbeEnd forbid all possible probes");
        }

        // Must be consistent with probeRegionCenteredAt().
        int centreOffset = regionCentreStartOffset(PROBE_LENGTH);

        // minProbeStart = initialPosition + offset + centreOffset
        int minOffset = minProbeStart - initialPosition.Position + centreOffset;

        // maxProbeEnd = initialPosition + offset + centreOffset + PROBE_LENGTH - 1
        int maxOffset = maxProbeEnd - initialPosition.Position + centreOffset - PROBE_LENGTH + 1;

        return outwardMovingOffsets(minOffset, maxOffset)
                .mapToObj(offset -> probeRegionCenteredAt(initialPosition.Chromosome, initialPosition.Position + offset));
    }

    // Generates the one best acceptable probe that is centered on one of the given positions.
    private ProbeGenerationResult coverOnePosition(List<BasePosition> positions, final TargetMetadata metadata,
            final ProbeEvaluator.Criteria evalCriteria, final ProbeSelector.Strategy selectStrategy, final PanelCoverage coverage)
    {
        Map<ChrBaseRegion, BasePosition> probeToPosition = new HashMap<>();
        List<Probe> candidateProbes = positions.stream()
                .map(position ->
                {
                    ChrBaseRegion probeRegion = probeRegionCenteredAt(position);
                    SequenceDefinition definition = SequenceDefinition.singleRegion(probeRegion);
                    TargetedRange targetedRange = TargetedRange.singleBase(regionCentreStartOffset(definition.baseLength()));
                    Probe probe = new Probe(definition, targetedRange, metadata);
                    // Store the target position so it can be retrieved later once the final probe is selected.
                    // Also needed to compute the rejected regions.
                    probeToPosition.put(probeRegion, position);
                    return probe;
                })
                .toList();

        Optional<Probe> bestCandidate = selectBestCandidate(candidateProbes.stream(), evalCriteria, selectStrategy);

        // Always include the full list of candidate positions in the result.
        // Perhaps surprising, but it's consistent with coverOneSubregion().
        // Also nothing should be precisely relying on the candidate targets output.
        List<TargetRegion> candidateTargetRegions = probeToPosition.values().stream()
                .map(position -> new TargetRegion(ChrBaseRegion.from(position), metadata))
                .toList();

        return bestCandidate
                .map(probe ->
                {
                    if(coverage.isCovered(probe.definition(), probe.targetedRange()))
                    {
                        return ProbeGenerationResult.alreadyCoveredTargets(candidateTargetRegions);
                    }
                    else
                    {
                        return new ProbeGenerationResult(List.of(probe), candidateTargetRegions, emptyList());
                    }
                })
                .orElseGet(() -> ProbeGenerationResult.rejectTargets(candidateTargetRegions));
    }

    private ProbeGenerationResult singleProbeBatch(List<ProbeGenerationSpec.SingleProbe> specs, final PanelCoverage coverage,
            PanelStore resultStore)
    {
        Stream<Probe> candidateProbes = specs.stream().map(spec ->
                new Probe(spec.sequenceDefinition(), spec.targetedRange(), spec.metadata()).withEvaluationCriteria(spec.evalCriteria()));
        Stream<Probe> evaluatedProbes = evaluateProbes(candidateProbes);
        return evaluatedProbes.map(probe ->
                {
                    ProbeGenerationResult result;
                    if(coverage.isCovered(probe.definition(), probe.targetedRange()))
                    {
                        result = ProbeGenerationResult.alreadyCoveredProbe(probe);
                    }
                    else
                    {
                        if(probe.accepted())
                        {
                            result = ProbeGenerationResult.acceptProbe(probe);
                        }
                        else
                        {
                            result = ProbeGenerationResult.rejectProbe(probe);
                        }
                    }
                    resultStore.addResult(result);
                    return result;
                })
                .reduce(new ProbeGenerationResult(), ProbeGenerationResult::add);
    }

    // Evaluates candidate probes and select the one best probe.
    private Optional<Probe> selectBestCandidate(Stream<Probe> candidates, final ProbeEvaluator.Criteria evalCriteria,
            final ProbeSelector.Strategy selectStrategy)
    {
        Stream<Probe> evaluatedCandidates = evaluateProbes(candidates, evalCriteria);
        return selectBestProbe(evaluatedCandidates, selectStrategy);
    }

    private Stream<Probe> evaluateProbes(Stream<Probe> probes, final ProbeEvaluator.Criteria criteria)
    {
        return evaluateProbes(probes.map(probe -> probe.withEvaluationCriteria(criteria)));
    }

    // Probes must have the evaluation criteria already set.
    private Stream<Probe> evaluateProbes(Stream<Probe> probes)
    {
        return mProbeEvaluator.evaluateProbes(probes).map(this::logCandidateProbe);
    }

    private Probe logCandidateProbe(final Probe probe)
    {
        if(mCandidateCallback != null)
        {
            mCandidateCallback.accept(probe);
        }
        return probe;
    }
}
