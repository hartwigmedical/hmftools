package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.maxProbeStartCovering;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.minProbeStartContaining;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.minProbeStartCovering;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.nextProbeStartPosition;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.computeUncoveredRegions;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

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
            final ProbeSelector.Criteria criteria, @Nullable final PanelCoverage coverage)
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
                .map(subregion -> coverSubregion(subregion, context, criteria))
                .reduce(new ProbeGenerationResult(), ProbeGenerationResult::add);
        // Add in the target region that's not added by coverSubregion().
        result = result.add(new ProbeGenerationResult(List.of(context.targetRegion()), emptyList(), emptyList()));
        return result;
    }

    private ProbeGenerationResult coverSubregion(final ChrBaseRegion region, final CandidateProbeContext context,
            final ProbeSelector.Criteria criteria)
    {
        // TODO: this is not good because
        //  - can make probe that only covers tiny bit of region
        //  - can make probe that extends far past region - should centre probes instead
        //  - shouldn't try to find best probe?

        // Methodology:
        //   - For each position in the target region, try to find the 1 best acceptable probe that covers it.
        //   - If an acceptable probe is found, advance the position to the next position after the probe.
        //   - If a position can't be covered by a probe, move to the next position.

        String chromosome = region.chromosome();
        BaseRegion targetBaseRegion = region.baseRegion();

        List<EvaluatedProbe> probes = new ArrayList<>();
        for(int position = targetBaseRegion.start(); position <= targetBaseRegion.end(); )
        {
            // Ensure:
            //  - Probe covers this position;
            //  - Limited overlap with the previous probe;
            //  - Minimum overlap with the target region
            int initialProbeStart = min(position, maxProbeStartCovering(targetBaseRegion));
            int minProbeStart = max(minProbeStartContaining(position), minProbeStartCovering(targetBaseRegion));
            if(!probes.isEmpty())
            {
                minProbeStart = max(minProbeStart, nextProbeStartPosition(probes.get(probes.size() - 1).candidate().probeRegion().end()));
            }

            if(initialProbeStart < minProbeStart)
            {
                // Not allowed to place any more probes due to the target region coverage constraint.
                break;
            }

            Stream<CandidateProbe> candidates = mCandidateGenerator.leftMovingLeftAlignedProbes(
                    new BasePosition(chromosome, initialProbeStart), minProbeStart, context);
            Optional<EvaluatedProbe> bestCandidate = mProbeSelector.selectBestCandidate(candidates, criteria);

            if(bestCandidate.isPresent())
            {
                probes.add(bestCandidate.get());
                position = bestCandidate.get().candidate().probeRegion().end() + 1;
            }
            else
            {
                ++position;
            }
        }

        // Compute rejected regions based on what has been covered by the probes.
        String rejectionReason = "No probe covering region meeting criteria " + criteria.eval();
        List<RejectedRegion> rejectedRegions = computeUncoveredRegions(
                targetBaseRegion, probes.stream().map(probe -> probe.candidate().probeRegion().baseRegion()))
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
