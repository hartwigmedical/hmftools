package com.hartwig.hmftools.panelbuilder;

import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CDR3_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CDR3_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CDR3_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionEndingAt;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionStartingAt;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.cider.IgTcrGene;
import com.hartwig.hmftools.common.cider.IgTcrGeneFile;
import com.hartwig.hmftools.common.cider.IgTcrRegion;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Fixed set of probes based on CDR3 regions.
// Methodology:
//   - 1 probe at the end of V regions
//   - 1 probe at the start of J regions
public class Cdr3Regions
{
    private static final ProbeEvaluator.Criteria PROBE_CRITERIA = new ProbeEvaluator.Criteria(
            CDR3_QUALITY_MIN, CDR3_GC_TARGET, CDR3_GC_TOLERANCE);

    private static final Logger LOGGER = LogManager.getLogger(Cdr3Regions.class);

    public static void generateProbes(final RefGenomeVersion refGenomeVersion, final ProbeGenerator probeGenerator, PanelData panelData)
    {
        LOGGER.info("Generating CDR3 probes");

        List<IgTcrGene> genes = IgTcrGeneFile.read(refGenomeVersion).stream()
                .filter(gene -> gene.region() == IgTcrRegion.V_REGION || gene.region() == IgTcrRegion.J_REGION)
                .filter(IgTcrGene::inPrimaryAssembly)
                .filter(gene -> gene.anchorLocation() != null)
                .toList();
        LOGGER.info("Loaded {} V/J genes", genes.size());

        ProbeGenerationResult result = new ProbeGenerationResult();
        ArrayList<ChrBaseRegion> coveredRegions = new ArrayList<>();
        for(IgTcrGene gene : genes)
        {
            result = result.add(generateProbe(gene, probeGenerator, panelData, coveredRegions));
        }

        panelData.addResult(result);

        LOGGER.info("Done generating CDR3 probes");
    }

    private static ProbeGenerationResult generateProbe(final IgTcrGene gene, final ProbeGenerator probeGenerator,
            final PanelCoverage coverage, ArrayList<ChrBaseRegion> coveredRegions)
    {
        ChrBaseRegion region = calculateProbeRegion(gene);

        if(coveredRegions.stream().anyMatch(region::overlaps))
        {
            // It's possible regions overlap other regions, in which case just take the first and discard the rest.
            LOGGER.trace("CDR3 region overlaps with another; discarding");
            return new ProbeGenerationResult();
        }
        coveredRegions.add(region);

        TargetMetadata metadata = createTargetMetadata(gene);
        TargetRegion target = new TargetRegion(region, metadata);
        Probe probe = probeGenerator.mProbeFactory.createProbeFromRegion(region, metadata).orElseThrow();
        probe = probeGenerator.mProbeEvaluator.evaluateProbe(probe, PROBE_CRITERIA);

        if(probe.accepted())
        {
            if(coverage.isCovered(region))
            {
                LOGGER.trace("CDR3 target already covered by panel: {}", target);
                return ProbeGenerationResult.alreadyCoveredTarget(target);
            }
            else
            {
                return ProbeGenerationResult.coveredTarget(target, probe);
            }
        }
        else
        {
            LOGGER.debug("No acceptable probe for CDR3 target: {}", target);
            String rejectionReason = "Probe does not meet criteria " + PROBE_CRITERIA;
            return ProbeGenerationResult.rejectTarget(target, rejectionReason);
        }
    }

    private static ChrBaseRegion calculateProbeRegion(final IgTcrGene gene)
    {
        ChrBaseRegion anchor = requireNonNull(gene.anchorLocation());
        boolean vForward = gene.region() == IgTcrRegion.V_REGION && gene.geneStrand() == Strand.FORWARD;
        boolean jReverse = gene.region() == IgTcrRegion.J_REGION && gene.geneStrand() == Strand.REVERSE;
        if(vForward || jReverse)
        {
            return probeRegionEndingAt(anchor.chromosome(), anchor.end());
        }
        else
        {
            return probeRegionStartingAt(anchor.chromosome(), anchor.start());
        }
    }

    private static TargetMetadata createTargetMetadata(final IgTcrGene gene)
    {
        String extraInfo = gene.geneName();
        return new TargetMetadata(TargetMetadata.Type.CDR3, extraInfo);
    }
}
