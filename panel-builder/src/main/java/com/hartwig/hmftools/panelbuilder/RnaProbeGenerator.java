package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_SHIFT_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.RNA_EXON_SINGLE_PROBE_SLACK;
import static com.hartwig.hmftools.panelbuilder.ProbeSelector.selectBestProbe;
import static com.hartwig.hmftools.panelbuilder.ProbeTiling.calculateContainedTiling;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionNegatedIntersection;
import static com.hartwig.hmftools.panelbuilder.Utils.outwardMovingOffsets;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.function.Consumer;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

// Exon-aware probe generation for RNA transcripts, driven off a RegionMapping in probe-space (the exons concatenated as a contiguous
// coordinate system). Placing probes in probe-space means a window near an exon edge maps naturally to a spliced (multi-region) probe, so no
// intronic bases are ever covered. Structurally the DNA probe-space analogue of ProbeGenerator.coverUncoveredRegion / coverAcceptableSubregion:
// generate all candidate windows overlapping the target, split into acceptable sub-ranges, then tile each by the exon size rules.
// Kept separate from the DNA generator for now; see RNA_DESIGN_NOTES for the eventual merge.
public class RnaProbeGenerator
{
    private final ProbeEvaluator mProbeEvaluator;
    @Nullable
    private final Consumer<Probe> mCandidateCallback;

    private static final Logger LOGGER = LogManager.getLogger(RnaProbeGenerator.class);

    public RnaProbeGenerator(final ProbeEvaluator probeEvaluator, final @Nullable Consumer<Probe> candidateCallback)
    {
        mProbeEvaluator = probeEvaluator;
        mCandidateCallback = candidateCallback;
    }

    // Covers a probe-space target range with the best acceptable probes, applying the exon size rules:
    //   - Range shorter than a probe: one probe padded out across the adjacent exon junction (spliced).
    //   - Range up to a probe plus the single-probe slack: one centred probe.
    //   - Longer: multiple probes tiled evenly, pinned flush to any edge that is a true exon boundary.
    // The range must lie within a single mapped region (exon); the caller enforces this. Junction crossing only occurs in the short-range
    // padding branch, never within-exon tiling. Candidate target regions are not added here (mirrors coverUncoveredRegion); the caller adds them.
    public ProbeGenerationResult coverExonRange(final RegionMapping mapping, int rangeStart, int rangeEnd, final TargetMetadata metadata,
            final ProbeEvaluator.Criteria evalCriteria, final ProbeSelector.Strategy selectStrategy)
    {
        if(!(rangeStart >= 0 && rangeStart < rangeEnd && rangeEnd <= mapping.length()))
        {
            throw new IllegalArgumentException("Invalid probe-space range");
        }
        if(mapping.toGenomeRegions(rangeStart, rangeEnd).size() != 1)
        {
            throw new IllegalArgumentException("Target range must lie within a single mapped region (exon)");
        }

        // Index all acceptable candidate windows by their probe-space start, and merge them into the maximal sub-ranges in which acceptable
        // probes can be placed. A rejected (bad QS/GC) window splits the target into separate acceptable sub-ranges.
        Map<Integer, Probe> acceptableProbes = new HashMap<>();
        List<BaseRegion> acceptableSubranges = new ArrayList<>();
        evaluate(allOverlappingProbes(mapping, rangeStart, rangeEnd, metadata), evalCriteria)
                .filter(Probe::accepted)
                .forEachOrdered(probe ->
                {
                    int start = probeSpaceStart(mapping, probe);
                    acceptableProbes.put(start, probe);
                    // Probe window in 1-indexed probe-space so the shared region maths apply directly.
                    BaseRegion window = toProbeSpaceRegion(start, start + PROBE_LENGTH);
                    BaseRegion prev = acceptableSubranges.isEmpty() ? null : acceptableSubranges.get(acceptableSubranges.size() - 1);
                    if(prev != null && window.start() <= prev.end() + 1)
                    {
                        if(window.end() > prev.end())
                        {
                            acceptableSubranges.set(acceptableSubranges.size() - 1, new BaseRegion(prev.start(), window.end()));
                        }
                    }
                    else
                    {
                        acceptableSubranges.add(window);
                    }
                });

        if(acceptableSubranges.size() > 1)
        {
            LOGGER.trace("Exon target range split into {} acceptable sub-ranges", acceptableSubranges.size());
        }

        List<Probe> probes = new ArrayList<>();
        List<BaseRegion> unrejected = new ArrayList<>();
        for(BaseRegion subrange : acceptableSubranges)
        {
            // Bound the sub-range (a union of probe windows, so wider than the target) back to the target.
            int tilingStart = max(fromProbeSpacePosition(subrange.start()), rangeStart);
            int tilingEnd = min(subrange.end(), rangeEnd);
            coverAcceptableSubrange(mapping, tilingStart, tilingEnd, acceptableProbes, selectStrategy, probes, unrejected);
        }

        List<RejectedFeature> rejectedFeatures = regionNegatedIntersection(
                toProbeSpaceRegion(rangeStart, rangeEnd), unrejected.stream())
                .stream()
                .flatMap(rejected -> mapping.toGenomeRegions(fromProbeSpacePosition(rejected.start()), rejected.end()).stream())
                .map(region -> RejectedFeature.fromRegion(region, metadata))
                .toList();

        return new ProbeGenerationResult(probes, emptyList(), rejectedFeatures);
    }

    // Tiles one acceptable sub-range (already bounded to the target) and appends the placed probes plus the ranges that should not be marked
    // rejected (probe coverage plus permitted uncovered edge slack).
    private void coverAcceptableSubrange(final RegionMapping mapping, int tilingStart, int tilingEnd,
            final Map<Integer, Probe> acceptableProbes, final ProbeSelector.Strategy selectStrategy, final List<Probe> probes,
            final List<BaseRegion> unrejected)
    {
        int length = tilingEnd - tilingStart;

        if(length < PROBE_LENGTH)
        {
            // Short range: pad a single probe out across an exon boundary to cover the whole range.
            selectPaddingProbe(mapping, tilingStart, tilingEnd, acceptableProbes, selectStrategy).ifPresent(probe ->
            {
                probes.add(probe);
                // The padding probe covers the entire target range (extending beyond the exon), so nothing here is rejected.
                unrejected.add(toProbeSpaceRegion(tilingStart, tilingEnd));
            });
            return;
        }

        // Pin an edge flush only when it coincides with a true exon boundary. An interior edge left by a rejected-region split is not a
        // boundary, so it stays a normal (unpinned) tiling edge with a small permitted uncovered gap.
        boolean pinStart = mapping.isRegionBoundary(tilingStart);
        boolean pinEnd = mapping.isRegionBoundary(tilingEnd);

        List<Integer> tiling;
        if(length <= PROBE_LENGTH + RNA_EXON_SINGLE_PROBE_SLACK)
        {
            // Single probe centred in the range (a single probe cannot pin both boundaries).
            tiling = List.of(tilingStart + (int) round((length - PROBE_LENGTH) / 2.0));
            pinStart = false;
            pinEnd = false;
        }
        else
        {
            tiling = calculateContainedTiling(toProbeSpaceRegion(tilingStart, tilingEnd), pinStart, pinEnd).stream()
                    .map(RnaProbeGenerator::fromProbeSpacePosition)
                    .toList();
        }

        List<Probe> placed = shiftAndSelect(tiling, tilingStart, tilingEnd, pinStart, pinEnd, acceptableProbes, selectStrategy);
        placed.forEach(probe ->
        {
            probes.add(probe);
            int start = probeSpaceStart(mapping, probe);
            unrejected.add(toProbeSpaceRegion(max(start, tilingStart), min(start + PROBE_LENGTH, tilingEnd)));
        });

        // Slack between an outermost probe and an unpinned edge is permitted uncovered (not rejected); pinned edges are flush so have none.
        if(!placed.isEmpty())
        {
            int firstStart = probeSpaceStart(mapping, placed.get(0));
            if(firstStart > tilingStart)
            {
                unrejected.add(toProbeSpaceRegion(tilingStart, firstStart));
            }
            int lastEnd = probeSpaceStart(mapping, placed.get(placed.size() - 1)) + PROBE_LENGTH;
            if(lastEnd < tilingEnd)
            {
                unrejected.add(toProbeSpaceRegion(lastEnd, tilingEnd));
            }
        }
    }

    // For each tiled position, shifts the probe within its allowed window (staying contained, keeping pinned edges flush, and not creating
    // gaps with neighbours) and selects the best acceptable probe. Mirrors ProbeGenerator.coverAcceptableSubregion but contained and
    // pin-aware. Positions with no acceptable probe in range are skipped.
    private static List<Probe> shiftAndSelect(final List<Integer> tiling, int tilingStart, int tilingEnd, boolean pinStart, boolean pinEnd,
            final Map<Integer, Probe> acceptableProbes, final ProbeSelector.Strategy selectStrategy)
    {
        List<Probe> placed = new ArrayList<>();
        Integer prevStart = null;
        for(int i = 0; i < tiling.size(); ++i)
        {
            int original = tiling.get(i);

            // Keep the probe contained within the range.
            int minStart = tilingStart;
            int maxStart = tilingEnd - PROBE_LENGTH;

            // Don't create a gap with the previously placed probe or the next probe's original tiled position.
            if(prevStart != null)
            {
                maxStart = min(maxStart, prevStart + PROBE_LENGTH);
            }
            if(i + 1 < tiling.size())
            {
                minStart = max(minStart, tiling.get(i + 1) - PROBE_LENGTH);
            }

            // Pinned outer edges are fixed flush to the boundary.
            if(i == 0 && pinStart)
            {
                minStart = maxStart = tilingStart;
            }
            if(i == tiling.size() - 1 && pinEnd)
            {
                minStart = maxStart = tilingEnd - PROBE_LENGTH;
            }

            // Limit the shift distance.
            if(!(i == 0 && pinStart) && !(i == tiling.size() - 1 && pinEnd))
            {
                minStart = max(minStart, original - PROBE_SHIFT_MAX);
                maxStart = min(maxStart, original + PROBE_SHIFT_MAX);
            }

            // If the constraints are contradictory, fall back to the original tiled position.
            if(minStart > maxStart)
            {
                minStart = maxStart = original;
            }

            Optional<Probe> best = selectBestFromPositions(original, minStart, maxStart, acceptableProbes, selectStrategy);
            if(best.isPresent())
            {
                placed.add(best.get());
                prevStart = original;
            }
        }
        return placed;
    }

    // Selects the best acceptable probe fully covering the short target range, preferring one centred over it (extending across the exon
    // boundary as needed). Empty if no acceptable padding probe exists.
    private Optional<Probe> selectPaddingProbe(final RegionMapping mapping, int tilingStart, int tilingEnd,
            final Map<Integer, Probe> acceptableProbes, final ProbeSelector.Strategy selectStrategy)
    {
        int length = tilingEnd - tilingStart;
        // Windows fully containing the range, within the mapping bounds.
        int minStart = max(0, tilingEnd - PROBE_LENGTH);
        int maxStart = min(mapping.length() - PROBE_LENGTH, tilingStart);
        int centred = tilingStart + (int) round((length - PROBE_LENGTH) / 2.0);
        return selectBestFromPositions(centred, minStart, maxStart, acceptableProbes, selectStrategy);
    }

    // Selects the best acceptable probe among the given probe-space start positions, visited outward from a preferred position so ties
    // favour the position closest to it.
    private static Optional<Probe> selectBestFromPositions(int preferredStart, int minStart, int maxStart,
            final Map<Integer, Probe> acceptableProbes, final ProbeSelector.Strategy selectStrategy)
    {
        Stream<Probe> candidates = outwardMovingOffsets(minStart - preferredStart, maxStart - preferredStart)
                .map(offset -> preferredStart + offset)
                .mapToObj(acceptableProbes::get)
                .filter(Objects::nonNull);
        return selectBestProbe(candidates, selectStrategy);
    }

    // Generates every candidate probe whose window overlaps the target range and lies within the mapping bounds, in ascending order.
    // Windows crossing an exon junction become spliced (multi-region) probes; windows within one exon stay single-region.
    private static Stream<Probe> allOverlappingProbes(final RegionMapping mapping, int rangeStart, int rangeEnd,
            final TargetMetadata metadata)
    {
        int minStart = max(0, rangeStart - PROBE_LENGTH + 1);
        int maxStart = min(mapping.length() - PROBE_LENGTH, rangeEnd - 1);
        return IntStream.rangeClosed(minStart, maxStart)
                .mapToObj(start ->
                {
                    List<ChrBaseRegion> regions = mapping.toGenomeRegions(start, start + PROBE_LENGTH);
                    SequenceDefinition definition = SequenceDefinition.spliced(regions);
                    TargetedRange targetedRange = new TargetedRange(
                            max(0, rangeStart - start), min(PROBE_LENGTH, rangeEnd - start));
                    return new Probe(definition, targetedRange, metadata);
                });
    }

    // Recovers a probe's probe-space start from its first (lowest, genome-forward) region.
    private static int probeSpaceStart(final RegionMapping mapping, final Probe probe)
    {
        ChrBaseRegion first = probe.definition().regions().get(0);
        return mapping.toProbeSpacePosition(first.chromosome(), first.start()).orElseThrow();
    }

    // Probe-space is 0-indexed half-open; the shared region maths use 1-indexed inclusive BaseRegions.
    private static BaseRegion toProbeSpaceRegion(int spaceStart, int spaceEnd)
    {
        return new BaseRegion(spaceStart + 1, spaceEnd);
    }

    private static int fromProbeSpacePosition(int baseRegionPosition)
    {
        return baseRegionPosition - 1;
    }

    private Stream<Probe> evaluate(Stream<Probe> probes, final ProbeEvaluator.Criteria criteria)
    {
        return mProbeEvaluator.evaluateProbes(probes.map(probe -> probe.withEvaluationCriteria(criteria)))
                .map(this::logCandidate);
    }

    private Probe logCandidate(final Probe probe)
    {
        if(mCandidateCallback != null)
        {
            mCandidateCallback.accept(probe);
        }
        return probe;
    }
}
