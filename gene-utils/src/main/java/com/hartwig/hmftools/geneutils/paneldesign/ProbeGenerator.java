package com.hartwig.hmftools.geneutils.paneldesign;

import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.calculateOptimalProbeTiling;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.maxProbeStartOverlapping;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.minProbeStartOverlapping;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.probeRegionStartingAt;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.computeUncoveredRegions;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.regionIntersection;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
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
    public final ProbeSelector mProbeSelector;

    private static final Logger LOGGER = LogManager.getLogger(ProbeGenerator.class);

    public ProbeGenerator(final CandidateProbeGenerator candidateGenerator, final ProbeSelector probeSelector)
    {
        mCandidateGenerator = candidateGenerator;
        mProbeSelector = probeSelector;
    }

    // Generates the best acceptable probes to cover an entire region. The probes may overlap and extend outside the target region.
    // If `coverage` is not null, avoid placing probes in already covered regions.
    public ProbeGenerationResult coverRegion(final ChrBaseRegion region, final CandidateProbeContext context,
            final ProbeEvaluator.Criteria criteria, @Nullable final PanelCoverage coverage)
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
            final ProbeEvaluator.Criteria criteria)
    {
        // Methodology:
        // Generate all possible probes and check which are acceptable, which gives us subregions which are acceptable and unacceptable.
        // Then tile probes within each of the acceptable regions.

        String chromosome = region.chromosome();
        BaseRegion baseRegion = region.baseRegion();

        int minPlausibleProbeStart = minProbeStartOverlapping(baseRegion);
        int maxPlausibleProbeStart = maxProbeStartOverlapping(baseRegion);
        Stream<EvaluatedProbe> allPlausibleProbes = IntStream.range(minPlausibleProbeStart, maxPlausibleProbeStart + 1)
                .mapToObj(start -> context.createProbe(probeRegionStartingAt(start)))
                .map(candidate -> mProbeSelector.mProbeEvaluator.evaluateCandidate(candidate, criteria))
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

                    return calculateOptimalProbeTiling(tilingTarget, probeBounds).stream()
                            // TODO: can try shifting probe a little and pick best probe
                            .map(acceptableProbes::get)
                            // TODO: is this valid? is it possible to have no acceptable probe here?
                            .map(Objects::requireNonNull);
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
            final ProbeSelector.Criteria criteria)
    {
        Stream<CandidateProbe> candidates = mCandidateGenerator.coverOneSubregion(region, context);
        return mProbeSelector.selectBestCandidate(candidates, criteria)
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
            final ProbeSelector.Criteria criteria)
    {
        Stream<CandidateProbe> candidates = mCandidateGenerator.coverPosition(position, context);
        return mProbeSelector.selectBestCandidate(candidates, criteria)
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
