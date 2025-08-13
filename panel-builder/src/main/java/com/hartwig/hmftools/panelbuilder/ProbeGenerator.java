package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.ceil;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.util.Collections.emptyList;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_OVERLAP_EXTENSION_BALANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_SHIFT_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.REGION_UNCOVERED_MAX;
import static com.hartwig.hmftools.panelbuilder.ProbeSelector.selectBestProbe;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.maxProbeEndOverlapping;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.maxProbeEndWithoutGap;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.minProbeStartOverlapping;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.minProbeStartWithoutGap;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionCenteredAt;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionEndingAt;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionStartingAt;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.computeUncoveredRegions;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.mergeOverlapAndAdjacentRegions;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionCentreFloat;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionIntersection;
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
    private final ProbeFactory mProbeFactory;
    private final CandidateProbeGenerator mCandidateGenerator;
    private final ProbeEvaluator mProbeEvaluator;

    private static final Logger LOGGER = LogManager.getLogger(ProbeGenerator.class);

    public ProbeGenerator(final ProbeFactory probeFactory, final CandidateProbeGenerator candidateGenerator,
            final ProbeEvaluator probeEvaluator)
    {
        mProbeFactory = probeFactory;
        mCandidateGenerator = candidateGenerator;
        mProbeEvaluator = probeEvaluator;
    }

    public static ProbeGenerator construct(final RefGenomeInterface refGenome, final ProbeQualityProfile probeQualityProfile,
            final ProbeQualityModel probeQualityModel, final Consumer<Probe> candidateCallback)
    {
        ProbeFactory probeFactory = new ProbeFactory(refGenome, probeQualityProfile, probeQualityModel);
        return new ProbeGenerator(
                probeFactory,
                new CandidateProbeGenerator(probeFactory, refGenome.chromosomeLengths()),
                new ProbeEvaluator(candidateCallback)
        );
    }

    // General purpose method for generating the best acceptable probes to cover an entire region.
    // Behaviour:
    //   - Prefer a probe set which is symmetrical and centered on the target region, unless this would violate constraints.
    //   - Avoid placing probes in already covered regions, if `coverage` is not null.
    //   - The edges of the region may be slightly uncovered (since the probes will capture a slightly larger region).
    //   - Probes may overlap and extend outside the target region.
    //   - Probes are shifted slightly to optimise for the selection criteria.
    public ProbeGenerationResult coverRegion(final ChrBaseRegion region, final TargetMetadata metadata,
            final ProbeEvaluator.Criteria evalCriteria, final ProbeSelector.Strategy localSelect, @Nullable final PanelCoverage coverage)
    {
        List<ChrBaseRegion> subregions;
        if(coverage == null)
        {
            subregions = List.of(region);
        }
        else
        {
            // Split the region into subregions to avoid overlap with regions already covered by probes.
            subregions = computeUncoveredRegions(region.baseRegion(), coverage.coveredRegions().map(ChrBaseRegion::baseRegion))
                    .stream()
                    .map(baseRegion -> ChrBaseRegion.from(region.chromosome(), baseRegion))
                    .toList();
            if(subregions.size() > 1)
            {
                subregions.forEach(subregion -> LOGGER.debug("Split region into uncovered subregion: {}", subregion));
            }
        }

        ProbeGenerationResult result = subregions.stream()
                .map(subregion -> coverUncoveredRegion(subregion, metadata, evalCriteria, localSelect))
                .reduce(new ProbeGenerationResult(), ProbeGenerationResult::add);

        // Add in the candidate target region that's not added by coverSubregion().
        TargetRegion candidateTarget = new TargetRegion(region, metadata);
        result = result.add(new ProbeGenerationResult(emptyList(), List.of(candidateTarget), emptyList(), emptyList()));

        return result;
    }

    private ProbeGenerationResult coverUncoveredRegion(final ChrBaseRegion region, final TargetMetadata metadata,
            final ProbeEvaluator.Criteria evalCriteria, final ProbeSelector.Strategy localSelect)
    {
        // Methodology:
        //   1. Generate all possible probes and check which are acceptable.
        //   2. Based on acceptable probes, generate regions in which acceptable probes can be tiled.
        //   3. Within each acceptable region, tile probes according to the ideal tiling algorithm.
        //   4. For each probe, try shifting it left and right slightly and pick the local best probe.

        String chromosome = region.chromosome();
        BaseRegion baseRegion = region.baseRegion();

        Stream<Probe> allPlausibleProbes = mCandidateGenerator.allOverlapping(region, metadata)
                .map(candidate -> mProbeEvaluator.evaluateProbe(candidate, evalCriteria))
                .filter(Probe::accepted);

        // These are the subregions in which probes can be placed.
        // Probes are guaranteed to be rejected if they overlap subregions between the acceptable subregions.
        List<BaseRegion> acceptableSubregions = new ArrayList<>();
        // Map from start position to probe.
        Map<Integer, Probe> acceptableProbes = new HashMap<>();
        // This requires sorted by position, but it's already in that order.
        allPlausibleProbes.forEach(probe ->
        {
            BaseRegion probeRegion = requireNonNull(probe.region()).baseRegion();
            BaseRegion prev = acceptableSubregions.isEmpty() ? null : acceptableSubregions.get(acceptableSubregions.size() - 1);
            if(prev != null && probeRegion.start() <= prev.end() + 1)
            {
                if(probeRegion.end() > prev.end())
                {
                    // Don't mutate in place because we borrowed the object from the probe.
                    acceptableSubregions.set(acceptableSubregions.size() - 1, new BaseRegion(prev.start(), probeRegion.end()));
                }
            }
            else
            {
                acceptableSubregions.add(probeRegion);
            }
            acceptableProbes.put(probeRegion.start(), probe);
        });

        if(acceptableSubregions.size() > 1)
        {
            acceptableSubregions.forEach(subregion -> LOGGER.trace("Split region into acceptable subregion: {}", subregion));
        }

        List<Probe> probes = new ArrayList<>();
        String rejectionReason = "No probe covering region, producing valid tiling, and meeting criteria " + evalCriteria;
        List<RejectedRegion> rejectedRegions = new ArrayList<>();
        acceptableSubregions.forEach(acceptableSubregion ->
        {
            CoverAcceptableSubregionResult result =
                    coverAcceptableSubregion(acceptableSubregion, baseRegion, acceptableProbes, localSelect);
            probes.addAll(result.probes());
            // Compute rejected regions based on what has been covered by the probes.
            // Not that the reference region here is not necessarily the target region, but the subregion of the target which the tiling
            // algorithm found to be optimal to cover. This may exclude a few bases on the edge which are uncovered. However, we want to
            // count those as covered and not mark them as rejected.
            Stream<BaseRegion> probeRegions = result.probes().stream().map(probe -> requireNonNull(probe.region()).baseRegion());
            computeUncoveredRegions(result.tilingIntendedCoverage(), probeRegions).stream()
                    .map(uncovered -> new RejectedRegion(ChrBaseRegion.from(chromosome, uncovered), metadata, rejectionReason))
                    .forEach(rejectedRegions::add);
        });

        // Compute covered target regions by merging all probe regions and intersecting with the desired target region.
        List<TargetRegion> coveredTargets = mergeOverlapAndAdjacentRegions(probes.stream().map(Probe::region)).stream()
                .map(covered -> new TargetRegion(regionIntersection(covered, region).orElseThrow(), metadata))
                .toList();

        // Candidate target is not added here because it would be added multiple times if there are multiple calls to this function for one
        // target region. Candidate target must be added by the caller.
        return new ProbeGenerationResult(probes, emptyList(), coveredTargets, rejectedRegions);
    }

    private record CoverAcceptableSubregionResult(
            List<Probe> probes,
            // Region that the tiling algorithm found was optimal to cover, taking into account the edges which may be left uncovered.
            BaseRegion tilingIntendedCoverage
    )
    {
    }

    private static CoverAcceptableSubregionResult coverAcceptableSubregion(final BaseRegion acceptableSubregion,
            final BaseRegion targetRegion, final Map<Integer, Probe> acceptableProbes, final ProbeSelector.Strategy localSelect)
    {
        BaseRegion subregion = acceptableSubregion;

        // Bound the subregion to the target region to prevent producing too many probes (our bounds for generating candidates
        // earlier were the least strict possible).
        // calculateOptimalProbeTiling() will produce probes extending past the target region if that's allowed and optimal.
        BaseRegion tilingTarget = regionIntersection(subregion, targetRegion).orElseThrow();

        // The acceptable subregions are maximal because we checked all probes which overlap the target region.
        // Thus, the subregion is bounded on both sides by unacceptable regions or completely off-target regions, and probes
        // cannot be placed outside it.
        BaseRegion probeBounds = subregion;

        List<BaseRegion> tiling = calculateOptimalProbeTiling(tilingTarget, probeBounds).stream()
                .map(ProbeUtils::probeRegionStartingAt)
                .toList();

        // For each probe, if there is space around the probe within the tiling, try shifting the probe around to find the best.
        List<Probe> finalProbes = new ArrayList<>();
        for(int i = 0; i < tiling.size(); ++i)
        {
            BaseRegion originalProbe = tiling.get(i);
            BaseRegion prevProbe = finalProbes.isEmpty()
                    ? null
                    : requireNonNull(finalProbes.get(finalProbes.size() - 1).region()).baseRegion();
            // We don't know exactly what the next probe will be but allow at least its original tiled position to be valid.
            BaseRegion nextProbe = i + 1 < tiling.size() ? tiling.get(i + 1) : null;

            // Shift must ensure:
            //  - Can't extend outside the hard bounds;
            //  - Can't reduce coverage of the target region;
            //  - Can't shift further than a configured amount.

            // Apply hard bounds and min coverage constraints.
            int minStart = max(probeBounds.start(), minProbeStartOverlapping(targetRegion));
            int maxEnd = min(probeBounds.end(), maxProbeEndOverlapping(targetRegion));

            // If the probe is adjacent to a probe, don't shift it such that it creates a gap.
            if(prevProbe != null)
            {
                maxEnd = min(maxEnd, maxProbeEndWithoutGap(prevProbe));
            }
            if(nextProbe != null)
            {
                minStart = max(minStart, minProbeStartWithoutGap(nextProbe));
            }

            // If the probe is on the edge, don't shift it such that it reduces the target coverage.
            if(i == 0)
            {
                if(originalProbe.start() < targetRegion.start())
                {
                    // First probe, extending before the target region:
                    //   - Can't shift left at all
                    //   - Can't shift right past the target start
                    minStart = originalProbe.start();
                    maxEnd = min(maxEnd, probeRegionStartingAt(targetRegion.start()).end());
                }
                else
                {
                    // First probe, starting after the target start:
                    //   - Can't shift left before the target start
                    //   - Can't shift right at all
                    minStart = max(minStart, targetRegion.start());
                    maxEnd = originalProbe.end();
                }
            }
            if(i == tiling.size() - 1)
            {
                if(originalProbe.end() > targetRegion.end())
                {
                    // Last probe, extending after the target region:
                    //   - Can't shift left before target end
                    //   - Can't shift right at all
                    minStart = max(minStart, probeRegionEndingAt(targetRegion.end()).start());
                    maxEnd = originalProbe.end();
                }
                else
                {
                    // Last probe, ending before the target end:
                    //   - Can't shift left at all
                    //   - Can't shift right past target end
                    minStart = originalProbe.start();
                    maxEnd = min(maxEnd, targetRegion.end());
                }
            }

            // Max shift amount constraint.
            minStart = max(minStart, originalProbe.start() - PROBE_SHIFT_MAX);
            maxEnd = min(maxEnd, originalProbe.end() + PROBE_SHIFT_MAX);

            int maxStart = probeRegionEndingAt(maxEnd).start();

            // In case our additional checks here fail, just use the original tiling position as we have no better alternative.
            if(minStart > maxStart)
            {
                minStart = maxStart = originalProbe.start();
            }

            // Iterate in an outward moving pattern so in the case of a tie, prefer probes closer to the original tiling position.
            int minOffset = minStart - originalProbe.start();
            int maxOffset = maxStart - originalProbe.start();
            Stream<Probe> candidates = outwardMovingOffsets(minOffset, maxOffset)
                    .map(offset -> originalProbe.start() + offset)
                    .mapToObj(acceptableProbes::get)
                    .filter(Objects::nonNull);
            Optional<Probe> bestProbe = selectBestProbe(candidates, localSelect);
            bestProbe.ifPresent(finalProbes::add);
        }

        // The region where tiling found was acceptable to cover.
        // This may be equal to the target region or smaller if the tiling found it was optimal to leave some bases uncovered at the edges.
        // We want the uncovered edges to count as covered though (don't want to mark as rejected) so we pass this back to the caller.
        BaseRegion intendedRegion = tiling.isEmpty()
                ? targetRegion
                : regionIntersection(targetRegion, new BaseRegion(tiling.get(0).start(), tiling.get(tiling.size() - 1).end())).orElse(null);

        return new CoverAcceptableSubregionResult(finalProbes, intendedRegion);
    }

    // Calculates the best probe tiling of a region.
    // Objectives:
    //   - Cover most of the region, with possibly a few bases on the edge uncovered;
    //   - Centre the tiling on the region;
    //   - Equally spaced probes;
    //   - Probes cannot extend outside `probeBounds`.
    // Returns the start positions of the probes.
    private static List<Integer> calculateOptimalProbeTiling(final BaseRegion region, final BaseRegion probeBounds)
    {
        if(!region.hasValidPositions() || !probeBounds.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
        if(!probeBounds.containsRegion(region))
        {
            // Probably a bug in the caller.
            throw new IllegalArgumentException("probeBounds forbids all possible tilings");
        }

        int maxProbeOffTarget = PROBE_LENGTH - 1;
        // Hard bounds on the region in which probes can be placed.
        BaseRegion tilingBounds = new BaseRegion(
                max(probeBounds.start(), region.start() - maxProbeOffTarget),
                min(probeBounds.end(), region.end() + maxProbeOffTarget));

        double centre = regionCentreFloat(region);

        // Lower bound is number of probes which completely cover the region, possibly excluding the max allowed uncovered bases.
        int minProbes = (int) ceil(max(1.0, region.baseLength() - REGION_UNCOVERED_MAX) / PROBE_LENGTH);

        // Upper bound is maximally overlapped and maximally extending outside the target region.
        int maxProbeOverlap = PROBE_LENGTH - 1;
        // maxProbes * PROBE_LENGTH - (maxProbes-1) * maxProbeOverlap <= regionSize + startExtend + endExtend
        int maxProbes = (tilingBounds.baseLength() - maxProbeOverlap) / (PROBE_LENGTH - maxProbeOverlap);

        if(minProbes > maxProbes)
        {
            // No tiling is possible given the constraints.
            return emptyList();
        }

        // The optimal number of probes is always the minimum possible, since we guaranteed the minimum acceptably covers the region,
        // and adding more probes will only increase overlap or extension for no gain.
        int probeCount = minProbes;
        // How many probe bases "left over" from covering the region?
        int extra = probeCount * PROBE_LENGTH - region.baseLength();

        if(probeCount <= 0)
        {
            // Can't even place a single probe with the given constraints.
            return emptyList();
        }

        // Calculate the start of the tiling and the space between each probe start position.
        // Then just tile the probes regularly from those parameters (since the constraints were checked earlier).
        double tilingStart;
        double probeStartSpacing;
        if(extra == 0)
        {
            // Perfect tiling.
            tilingStart = region.start();
            probeStartSpacing = PROBE_LENGTH;
        }
        else
        {
            // General case.

            // If extra is negative (i.e. some bases uncovered) then reduce the desired tiling region.
            // If extra is positive (i.e. some overlap) then distribute that between overlap and extension based on the constant.
            double tilingSpace = min(
                    region.baseLength() + min(0, extra) + max(0, extra * PROBE_OVERLAP_EXTENSION_BALANCE),
                    tilingBounds.baseLength());
            probeStartSpacing = probeCount <= 1
                    ? 0.0
                    : max(
                            (tilingSpace - PROBE_LENGTH) / (probeCount - 1.0),
                            PROBE_LENGTH - maxProbeOverlap);
            double tilingLength = (probeCount - 1) * probeStartSpacing + PROBE_LENGTH;
            tilingStart = centre - tilingLength / 2;

            // Adjust the tiling alignment to adhere to the hard bounds.
            tilingStart = max(tilingStart, tilingBounds.start());
            double end = tilingStart + tilingLength - 1;
            if(end > tilingBounds.end())
            {
                tilingStart -= end - tilingBounds.end();
            }
        }

        double tilingStartCopy = tilingStart;
        List<Integer> probes = IntStream.range(0, probeCount)
                .map(i -> (int) round(tilingStartCopy + i * probeStartSpacing))
                .boxed().toList();
        return probes;
    }

    // Generates the one best acceptable probe that is contained within the specified region.
    public ProbeGenerationResult coverOneSubregion(final ChrBaseRegion region, final TargetMetadata metadata,
            final ProbeEvaluator.Criteria evalCriteria, final ProbeSelector.Strategy selectStrategy, @Nullable final PanelCoverage coverage)
    {
        // TODO: is it correct for the target to be the whole region?
        TargetRegion target = new TargetRegion(region, metadata);
        Stream<Probe> candidates = mCandidateGenerator.coverOneSubregion(region, metadata);
        return selectBestCandidate(candidates, evalCriteria, selectStrategy)
                .map(probe ->
                {
                    if(coverage != null && coverage.isCovered(requireNonNull(probe.region())))
                    {
                        return ProbeGenerationResult.alreadyCoveredTarget(target);
                    }
                    else
                    {
                        return ProbeGenerationResult.coveredTarget(target, probe);
                    }
                })
                .orElseGet(() ->
                {
                    String rejectionReason = "No probe in region meeting criteria " + evalCriteria;
                    return ProbeGenerationResult.rejectTarget(target, rejectionReason);
                });
    }

    // Generates the one best acceptable probe that centered on one of the given positions.
    public ProbeGenerationResult coverOnePosition(Stream<BasePosition> positions, final TargetMetadata metadata,
            final ProbeEvaluator.Criteria evalCriteria, final ProbeSelector.Strategy selectStrategy, final PanelCoverage coverage)
    {
        Map<ChrBaseRegion, BasePosition> probeToPosition = new HashMap<>();
        Stream<Probe> candidates = positions
                .flatMap(position ->
                {
                    ChrBaseRegion probeRegion = probeRegionCenteredAt(position);
                    Optional<Probe> probe = mProbeFactory.createProbeFromRegion(probeRegion, metadata);
                    if(probe.isPresent())
                    {
                        // Store the target position so it can be retrieved later once the final probe is selected.
                        // Also needed to compute the rejected regions.
                        probeToPosition.put(probeRegion, position);
                    }
                    return probe.stream();
                });

        return selectBestCandidate(candidates, evalCriteria, selectStrategy)
                .map(probe ->
                {
                    BasePosition position = probeToPosition.get(probe.region());
                    TargetRegion target = new TargetRegion(ChrBaseRegion.from(position), probe.metadata());
                    // TODO: is it correct to have only the probe as the candidate target?
                    if(coverage.isCovered(target.region()))
                    {
                        return ProbeGenerationResult.alreadyCoveredTarget(target);
                    }
                    else
                    {
                        return ProbeGenerationResult.coveredTarget(target, probe);
                    }
                })
                .orElseGet(() ->
                {
                    // Produce a rejected region for every position. Not sure a better way to handle it. No rejected region seems wrong.
                    List<TargetRegion> targets = probeToPosition.values().stream()
                            .map(position -> new TargetRegion(ChrBaseRegion.from(position), metadata))
                            .toList();
                    String rejectionReason = "Probe at position does not meet criteria " + evalCriteria;
                    return ProbeGenerationResult.rejectTargets(targets, rejectionReason);
                });
    }

    // Generates a probe with the given sequence.
    // `regions` are the target regions which the probe (partially) hits, if any.
    // Coverage check must be done by the caller since it's not necessarily known if this custom sequence is already covered by a probe.
    public ProbeGenerationResult probeWithSequence(final String sequence, final List<ChrBaseRegion> regions,
            final TargetMetadata metadata, final ProbeEvaluator.Criteria evalCriteria, boolean covered)
    {
        List<TargetRegion> targetRegions = regions.stream().map(region -> new TargetRegion(region, metadata)).toList();

        if(covered)
        {
            return ProbeGenerationResult.alreadyCoveredTargets(targetRegions);
        }
        else
        {
            return mProbeFactory.createProbeFromSequence(sequence, metadata)
                    .map(probe ->
                    {
                        probe = mProbeEvaluator.evaluateProbe(probe, evalCriteria);
                        if(probe.accepted())
                        {
                            return new ProbeGenerationResult(List.of(probe), targetRegions, targetRegions, emptyList());
                        }
                        else
                        {
                            String rejectionReason = "Probe does not meet criteria " + evalCriteria;
                            return ProbeGenerationResult.rejectTargets(targetRegions, rejectionReason);
                        }
                    })
                    .orElseGet(ProbeGenerationResult::new);
        }
    }

    // Evaluates candidate probes and select the one best probe.
    public Optional<Probe> selectBestCandidate(Stream<Probe> candidates, final ProbeEvaluator.Criteria evalCriteria,
            final ProbeSelector.Strategy selectStrategy)
    {
        Stream<Probe> evaluatedCandidates = mProbeEvaluator.evaluateProbes(candidates, evalCriteria);
        return selectBestProbe(evaluatedCandidates, selectStrategy);
    }
}
