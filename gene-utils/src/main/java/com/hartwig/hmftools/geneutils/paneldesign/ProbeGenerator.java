package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_SHIFT_MAX;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeSelector.selectBestProbe;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.calculateOptimalProbeTiling;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.maxProbeEndPartiallyCovering;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.maxProbeStartOverlapping;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.minProbeStartPartiallyCovering;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.minProbeStartOverlapping;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.nextProbeStartPosition;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.prevProbeEndPosition;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.probeRegionEndingAt;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.probeRegionStartingAt;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.computeUncoveredRegions;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.outwardMovingOffsets;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.regionIntersection;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

// TODO: unit test

public class ProbeGenerator
{
    public final CandidateProbeGenerator mCandidateGenerator;
    public final ProbeEvaluator mProbeEvaluator;

    private static final Logger LOGGER = LogManager.getLogger(ProbeGenerator.class);

    public ProbeGenerator(final CandidateProbeGenerator candidateGenerator, final ProbeEvaluator probeEvaluator)
    {
        mCandidateGenerator = candidateGenerator;
        mProbeEvaluator = probeEvaluator;
    }

    // Generates the best acceptable probes to cover an entire region. The probes may overlap and extend outside the target region.
    // If `coverage` is not null, avoid placing probes in already covered regions.
    public ProbeGenerationResult coverRegion(final ChrBaseRegion region, final CandidateProbeContext context,
            final ProbeSelectCriteria criteria, @Nullable final PanelCoverage coverage)
    {
        List<ChrBaseRegion> subregions;
        if(coverage == null)
        {
            subregions = List.of(region);
        }
        else
        {
            // Split the region into subregions to avoid overlap with regions already covered by probes.
            subregions =
                    computeUncoveredRegions(region.baseRegion(), coverage.coveredRegions().map(ChrBaseRegion::baseRegion))
                            .stream()
                            .map(baseRegion -> ChrBaseRegion.from(region.chromosome(), baseRegion))
                            .toList();
            if(subregions.size() > 1)
            {
                subregions.forEach(subregion -> LOGGER.debug("Split region into uncovered subregion: {}", subregion));
            }
        }

        ProbeGenerationResult result = subregions.stream()
                .map(subregion -> coverUncoveredRegion(subregion, context, criteria))
                .reduce(new ProbeGenerationResult(), ProbeGenerationResult::add);
        // Add in the target region that's not added by coverSubregion().
        result = result.add(new ProbeGenerationResult(List.of(context.targetRegion()), emptyList(), emptyList()));
        return result;
    }

    private ProbeGenerationResult coverUncoveredRegion(final ChrBaseRegion region, final CandidateProbeContext context,
            final ProbeSelectCriteria criteria)
    {
        // Methodology:
        //   1. Generate all possible probes and check which are acceptable.
        //   2. Based on acceptable probes, generate regions in which acceptable probes can be tiled.
        //   3. Within each acceptable region, tile probes according to the ideal tiling algorithm.
        //   4. For each probe, try shifting it left and right slightly and pick the local best probe.

        String chromosome = region.chromosome();
        BaseRegion baseRegion = region.baseRegion();

        int minPlausibleProbeStart = minProbeStartOverlapping(baseRegion);
        int maxPlausibleProbeStart = maxProbeStartOverlapping(baseRegion);
        Stream<EvaluatedProbe> allPlausibleProbes = IntStream.rangeClosed(minPlausibleProbeStart, maxPlausibleProbeStart)
                .mapToObj(start -> context.createProbe(probeRegionStartingAt(start)))
                .map(candidate -> mProbeEvaluator.evaluateCandidate(candidate, criteria.eval()))
                .filter(EvaluatedProbe::accepted);

        // These are the subregions in which probes can be placed.
        // Probes are guaranteed to be rejected if they overlap subregions between the acceptable subregions.
        List<BaseRegion> acceptableSubregions = new ArrayList<>();
        Map<Integer, EvaluatedProbe> acceptableProbes = new HashMap<>();
        // This requires sorted by position, but it's already in that order.
        allPlausibleProbes.forEach(probe ->
        {
            BaseRegion probeRegion = probe.candidate().probeRegion().baseRegion();
            BaseRegion prev = acceptableSubregions.isEmpty() ? null : acceptableSubregions.get(acceptableSubregions.size() - 1);
            if(prev != null && probeRegion.start() <= prev.end() + 1)
            {
                prev.setEnd(probeRegion.end());
            }
            else
            {
                acceptableSubregions.add(probeRegion);
            }
            acceptableProbes.put(probeRegion.start(), probe);
        });

        List<EvaluatedProbe> probes = acceptableSubregions.stream()
                .flatMap(acceptableSubregion ->
                {
                    BaseRegion subregion = acceptableSubregion;

                    // Bound the subregion to the target region to prevent producing too many probes (our bounds for generating candidates
                    // earlier were the least strict possible).
                    // calculateOptimalProbeTiling() will produce probes extending past the target region if that's allowed and optimal.
                    BaseRegion tilingTarget = regionIntersection(subregion, baseRegion);

                    // The acceptable subregions are maximal because we checked all probes which overlap the target region.
                    // Thus, the subregion is bounded on both sides by unacceptable regions or completely off-target regions, and probes
                    // cannot be placed outside it.
                    BaseRegion probeBounds = subregion;

                    // For each probe, if there is space around the probe within the tiling, try shifting the probe around to find the best.
                    List<Integer> tiling = calculateOptimalProbeTiling(tilingTarget, probeBounds);
                    List<EvaluatedProbe> tiledProbes = new ArrayList<>();
                    for(int i = 0; i < tiling.size(); ++i)
                    {
                        BaseRegion originalProbe = probeRegionStartingAt(tiling.get(i));

                        // Shift must ensure:
                        //  - Can't extend outside the hard bounds;
                        //  - Can't overlap the adjacent probes too much;
                        //  - Can't cover too few target bases;
                        //  - Can't reduce coverage of the target region;
                        //  - Can't shift further than a configured amount.

                        // Apply hard bounds and min coverage constraints.
                        int minStart = max(probeBounds.start(), minProbeStartPartiallyCovering(baseRegion));
                        int maxEnd = min(probeBounds.end(), maxProbeEndPartiallyCovering(baseRegion));

                        // TODO: possible to avoid increasing the total overlap?
                        // Apply probe overlap constraint.
                        // For the previous adjacent probe, use the actual probe selected.
                        // For the next adjacent probe, use the original tiling position to ensure at least 1 option for the next iteration.
                        if(!tiledProbes.isEmpty())
                        {
                            minStart = max(minStart,
                                    nextProbeStartPosition(tiledProbes.get(tiledProbes.size() - 1).candidate().probeRegion().end()));
                        }
                        if(i + 1 < tiling.size())
                        {
                            maxEnd = min(maxEnd, prevProbeEndPosition(tiling.get(i + 1)));
                        }

                        // If the probe is on the edge, don't shift it such that it reduces the target coverage.
                        if(i == 0)
                        {
                            if(originalProbe.start() < baseRegion.start())
                            {
                                // First probe, extending before the target region:
                                //   - Can't shift left at all
                                //   - Can't shift right past the target start
                                minStart = originalProbe.start();
                                maxEnd = min(maxEnd, probeRegionStartingAt(baseRegion.start()).end());
                            }
                            else
                            {
                                // First probe, starting after the target start:
                                //   - Can't shift left before the target start
                                //   - Can't shift right at all
                                minStart = max(minStart, baseRegion.start());
                                maxEnd = originalProbe.end();
                            }
                        }
                        if(i == tiling.size() - 1)
                        {
                            // TODO: technically could shift right more if there was overlap with the previous probe
                            if(originalProbe.end() > baseRegion.end())
                            {
                                // Last probe, extending after the target region:
                                //   - Can't shift left before target end
                                //   - Can't shift right at all
                                minStart = max(minStart, probeRegionEndingAt(baseRegion.end()).start());
                                maxEnd = originalProbe.end();
                            }
                            else
                            {
                                // Last probe, ending before the target end:
                                //   - Can't shift left at all
                                //   - Can't shift right past target end
                                minStart = originalProbe.start();
                                maxEnd = min(maxEnd, baseRegion.end());
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
                        Stream<EvaluatedProbe> candidates = outwardMovingOffsets(minOffset, maxOffset)
                                .map(offset -> originalProbe.start() + offset)
                                .mapToObj(acceptableProbes::get)
                                .filter(Objects::nonNull);
                        Optional<EvaluatedProbe> bestProbe = selectBestProbe(candidates, criteria.select());
                        bestProbe.ifPresent(tiledProbes::add);
                    }

                    return tiledProbes.stream();
                })
                .toList();

        // Compute rejected regions based on what has been covered by the probes.
        String rejectionReason = "No probe covering region, producing valid tiling, and meeting criteria " + criteria;
        List<RejectedRegion> rejectedRegions =
                computeUncoveredRegions(region.baseRegion(), probes.stream().map(probe -> probe.candidate().probeRegion().baseRegion()))
                        .stream()
                        .map(r -> new RejectedRegion(ChrBaseRegion.from(chromosome, r), context.targetRegion(), rejectionReason))
                        .toList();

        // Target is not added here because it would be added multiple times if there are multiple calls to this function for 1 target.
        // Target must be added by the caller.
        return new ProbeGenerationResult(emptyList(), probes, rejectedRegions);
    }

    // Generates the 1 best acceptable probe that is contained within the specified region.
    public ProbeGenerationResult coverOneSubregion(final ChrBaseRegion region, final CandidateProbeContext context,
            final ProbeSelectCriteria criteria)
    {
        Stream<CandidateProbe> candidates = mCandidateGenerator.coverOneSubregion(region, context);
        Stream<EvaluatedProbe> evaluatedCandidates = mProbeEvaluator.evaluateCandidates(candidates, criteria.eval());
        return selectBestProbe(evaluatedCandidates, criteria.select())
                .map(probe ->
                        new ProbeGenerationResult(List.of(context.targetRegion()), List.of(probe), emptyList()))
                .orElseGet(() ->
                {
                    String rejectionReason = "No probe in region meeting criteria " + criteria.eval();
                    return new ProbeGenerationResult(
                            List.of(context.targetRegion()),
                            emptyList(),
                            List.of(RejectedRegion.fromTargetRegion(context.targetRegion(), rejectionReason)));
                });
    }

    // Generates the 1 best acceptable probe which covers a position.
    public ProbeGenerationResult coverPosition(final BasePosition position, final CandidateProbeContext context,
            final ProbeSelectCriteria criteria)
    {
        Stream<CandidateProbe> candidates = mCandidateGenerator.coverPosition(position, context);
        Stream<EvaluatedProbe> evaluatedCandidates = mProbeEvaluator.evaluateCandidates(candidates, criteria.eval());
        return selectBestProbe(evaluatedCandidates, criteria.select())
                .map(probe -> new ProbeGenerationResult(List.of(context.targetRegion()), List.of(probe), emptyList()))
                .orElseGet(() ->
                {
                    String rejectionReason = "No probe covering position meeting criteria " + criteria.eval();
                    return new ProbeGenerationResult(
                            List.of(context.targetRegion()),
                            emptyList(),
                            List.of(RejectedRegion.fromTargetRegion(context.targetRegion(), rejectionReason)));
                });
    }
}
