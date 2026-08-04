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
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
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

    // Object to encapsulate batch generation of probes.
    // Generating in batches is more efficient because of computing the probe quality model (alignment).
    public class Batch
    {
        private final List<ProbeGenerationSpec> mGenericSpecs = new ArrayList<>();
        private final List<ProbeGenerationSpec.SingleProbe> mSingleProbeSpecs = new ArrayList<>();

        public void add(final ProbeGenerationSpec spec)
        {
            if(spec instanceof ProbeGenerationSpec.SingleProbe singleProbeSpec)
            {
                mSingleProbeSpecs.add(singleProbeSpec);
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
                    .add(generateSingleProbeSpecs(resultBuffer, resultBuffer));
            clear();
            return result;
        }

        private ProbeGenerationResult generateGenericSpecs(final PanelCoverage coverage, PanelStore resultStore)
        {
            // TODO: batch probe evaluation
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

        // DNA runs the shared probe-space coverage over an identity (whole-chromosome) mapping: probe-space position = genome position - 1.
        // No tiling edges are pinned (chromosome ends and target edges are arbitrary, not feature boundaries), so it reproduces the
        // genome-space tiling exactly.
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

    // Covers a probe-space target range within a single exon of a transcript's exon mapping (RNA). Tiling edges that coincide with an exon
    // boundary are pinned flush (good splice-junction coverage); a probe window near an exon edge maps to a spliced (multi-region) probe.
    // The range must lie within a single mapped exon. Candidate target regions are not added here; the caller adds them.
    public ProbeGenerationResult coverExonRange(final RegionMapping mapping, int rangeStart, int rangeEnd, final TargetMetadata metadata,
            final ProbeEvaluator.Criteria evalCriteria, final ProbeSelector.Strategy localSelect)
    {
        if(!(rangeStart >= 0 && rangeStart < rangeEnd && rangeEnd <= mapping.length()))
        {
            throw new IllegalArgumentException("Invalid probe-space range");
        }
        if(mapping.toGenomeRegions(rangeStart, rangeEnd).size() != 1)
        {
            throw new IllegalArgumentException("Target range must lie within a single mapped region (exon)");
        }
        return coverMappedRange(mapping, rangeStart, rangeEnd, metadata, evalCriteria, localSelect, mapping::isRegionBoundary);
    }

    // Probe-space analogue of coverUncoveredRegion, driven off a RegionMapping (identity mapping for DNA, exon mapping for RNA). The target
    // range [rangeStart, rangeEnd) is in probe-space. Candidate windows map back to genome via the mapping (a single region for DNA, possibly
    // spliced for RNA). Candidate target regions are not added here; the caller adds them. Probe-space positions are wrapped as 1-indexed
    // BaseRegions (position + 1) so the shared BaseRegion maths apply directly.
    private ProbeGenerationResult coverMappedRange(final RegionMapping mapping, int rangeStart, int rangeEnd, final TargetMetadata metadata,
            final ProbeEvaluator.Criteria evalCriteria, final ProbeSelector.Strategy localSelect, final IntPredicate pinBoundary)
    {
        // Accepted candidate windows indexed by probe-space start (as a 1-indexed coordinate), merged into maximal acceptable sub-ranges.
        Map<Integer, Probe> acceptableProbes = new HashMap<>();
        List<BaseRegion> acceptableSubregions = new ArrayList<>();
        evaluateProbes(mappedOverlappingProbes(mapping, rangeStart, rangeEnd, metadata), evalCriteria)
                .filter(Probe::accepted)
                .forEachOrdered(probe ->
                {
                    BaseRegion window = spaceWindow(probeSpaceStart(mapping, probe));
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
                    acceptableProbes.put(window.start(), probe);
                });

        BaseRegion targetRegion = spaceRegion(rangeStart, rangeEnd);

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

    // Probe-space analogue of coverAcceptableSubregion, pin-aware. A tiling edge that coincides with a mapping boundary (exon boundary for
    // RNA; never for the DNA identity mapping) is pinned flush; otherwise the DNA shift constraints apply (probes may extend past the target).
    private CoverAcceptableSubregionResult coverAcceptableMappedSubrange(final RegionMapping mapping, final BaseRegion acceptableSubregion,
            final BaseRegion targetRegion, final Map<Integer, Probe> acceptableProbes, final ProbeSelector.Strategy localSelect,
            final IntPredicate pinBoundary)
    {
        BaseRegion tilingTarget = regionIntersection(acceptableSubregion, targetRegion).orElseThrow();
        BaseRegion probeBounds = acceptableSubregion;

        // Pinning flushes the outermost probes to feature boundaries, which is only meaningful when tiling multiple probes. A range that fits
        // in a single probe (within the uncovered budget) is centred instead - it cannot be flush to both boundaries, and a short range is
        // centred and padded across the boundary. This is inert for DNA (never pins).
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

    // Generates every candidate probe whose window overlaps the probe-space target range and lies within the mapping bounds, ascending.
    private static Stream<Probe> mappedOverlappingProbes(final RegionMapping mapping, int rangeStart, int rangeEnd,
            final TargetMetadata metadata)
    {
        int minStart = max(0, rangeStart - PROBE_LENGTH + 1);
        int maxStart = min(mapping.length() - PROBE_LENGTH, rangeEnd - 1);
        return IntStream.rangeClosed(minStart, maxStart)
                .mapToObj(start ->
                {
                    List<ChrBaseRegion> regions = mapping.toGenomeRegions(start, start + PROBE_LENGTH);
                    SequenceDefinition definition = SequenceDefinition.spliced(regions);
                    TargetedRange targetedRange = new TargetedRange(max(0, rangeStart - start), min(PROBE_LENGTH, rangeEnd - start));
                    return new Probe(definition, targetedRange, metadata);
                });
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
