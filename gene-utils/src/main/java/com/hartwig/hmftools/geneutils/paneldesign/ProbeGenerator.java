package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.maxProbeStartCovering;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.minProbeStartContaining;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.minProbeStartCovering;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeUtils.nextProbeStartPosition;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.computeUncoveredRegions;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class ProbeGenerator
{
    public final CandidateProbeGenerator mCandidateGenerator;
    public final ProbeSelector mProbeSelector;

    public ProbeGenerator(final CandidateProbeGenerator candidateGenerator, final ProbeSelector probeSelector)
    {
        mCandidateGenerator = candidateGenerator;
        mProbeSelector = probeSelector;
    }

    // Generates the best acceptable probes to cover an entire region. The probes may overlap and extend outside the target region.
    public ProbeGenerationResult coverRegion(final TargetRegion target, final ProbeSelectCriteria criteria)
    {
        // TODO: this is not good because
        //  - can make probe that only covers tiny bit of region
        //  - can make probe that extends far past region - should centre probes instead
        //  - shouldn't try to find best probe?

        // Methodology:
        //   - For each position in the target, try to find the 1 best acceptable probe that covers it.
        //   - If an acceptable probe is found, advance the position to the next position after the probe.
        //   - If a position can't be covered by a probe, move to the next position.

        ChrBaseRegion targetChrBaseRegion = target.region();
        String chromosome = targetChrBaseRegion.chromosome();
        BaseRegion targetBaseRegion = targetChrBaseRegion.baseRegion();

        CandidateProbeFactory candidateFactory = new CandidateProbeFactory(target);
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
                    new BasePosition(chromosome, initialProbeStart), minProbeStart, candidateFactory);
            Optional<EvaluatedProbe> bestCandidate = mProbeSelector.selectBestProbe(candidates, criteria);

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
                .map(region ->
                        new RejectedRegion(ChrBaseRegion.from(chromosome, region), target, rejectionReason))
                .toList();

        return new ProbeGenerationResult(List.of(target), probes, rejectedRegions);
    }

    // Generates the 1 best acceptable probe that is contained within the specified region.
    public ProbeGenerationResult coverOneSubregion(final TargetRegion target, final ProbeSelectCriteria criteria)
    {
        return mProbeSelector.selectBestProbe(mCandidateGenerator.coverOneSubregion(target), criteria)
                .map(probe -> new ProbeGenerationResult(List.of(target), List.of(probe), Collections.emptyList()))
                .orElseGet(() ->
                {
                    String rejectionReason = "No probe in region meeting criteria " + criteria.eval();
                    return new ProbeGenerationResult(
                            List.of(target),
                            Collections.emptyList(),
                            List.of(RejectedRegion.fromTargetRegion(target, rejectionReason)));
                });
    }

    // Generates the 1 best acceptable probe which covers a position.
    public ProbeGenerationResult coverPosition(final BasePosition position, final TargetMetadata metadata,
            final ProbeSelectCriteria criteria)
    {
        TargetRegion target = new TargetRegion(ChrBaseRegion.from(position), metadata);
        return mProbeSelector.selectBestProbe(mCandidateGenerator.coverPosition(position, metadata), criteria)
                .map(probe -> new ProbeGenerationResult(List.of(target), List.of(probe), Collections.emptyList()))
                .orElseGet(() ->
                {
                    String rejectionReason = "No probe covering position meeting criteria " + criteria.eval();
                    return new ProbeGenerationResult(
                            List.of(target),
                            Collections.emptyList(),
                            List.of(RejectedRegion.fromTargetRegion(target, rejectionReason)));
                });
    }
}
