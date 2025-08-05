package com.hartwig.hmftools.panelbuilder;

import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CDR3_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENERAL_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENERAL_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionEndingAt;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionStartingAt;

import java.util.stream.Stream;

import com.hartwig.hmftools.common.cider.IgTcrGene;
import com.hartwig.hmftools.common.cider.IgTcrGeneFile;
import com.hartwig.hmftools.common.cider.IgTcrRegion;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Fixed set of probes based on CDR3 regions.
// TODO: describe the methodology
public class Cdr3Regions
{
    // TODO: confirm constraints
    private static final ProbeEvaluator.Criteria PROBE_CRITERIA = new ProbeEvaluator.Criteria(
            CDR3_QUALITY_MIN, GENERAL_GC_TARGET, GENERAL_GC_TOLERANCE);

    private static final Logger LOGGER = LogManager.getLogger(Cdr3Regions.class);

    public static void generateProbes(final RefGenomeVersion refGenomeVersion, final ProbeGenerator probeGenerator, PanelData panelData)
    {
        LOGGER.info("Generating CDR3 probes");

        Stream<IgTcrGene> genes = IgTcrGeneFile.read(refGenomeVersion).stream()
                .filter(gene -> gene.region() == IgTcrRegion.V_REGION || gene.region() == IgTcrRegion.J_REGION)
                .filter(IgTcrGene::inPrimaryAssembly)
                .filter(gene -> gene.anchorLocation() != null);

        ProbeGenerationResult result = genes
                .map(gene -> generateProbe(gene, probeGenerator, panelData))
                .reduce(new ProbeGenerationResult(), ProbeGenerationResult::add);

        panelData.addResult(result);

        LOGGER.info("Done generating CDR3 probes");
    }

    private static ProbeGenerationResult generateProbe(final IgTcrGene gene, final ProbeGenerator probeGenerator,
            final PanelCoverage coverage)
    {
        ChrBaseRegion region = calculateProbeRegion(gene);
        TargetMetadata metadata = createTargetMetadata(gene);
        TargetRegion target = new TargetRegion(region, metadata);
        Probe probe = probeGenerator.mProbeFactory.createProbeFromRegion(region, metadata).orElseThrow();
        probe = probeGenerator.mProbeEvaluator.evaluateProbe(probe, PROBE_CRITERIA);

        if(probe.accepted())
        {
            // Note this coverage check also checks previously added CDR3 regions, since overlaps may result.
            if(coverage.isCovered(region))
            {
                LOGGER.debug("CDR3 target already covered by panel: {}", target);
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
            return ChrBaseRegion.from(anchor.chromosome(), probeRegionEndingAt(anchor.end()));
        }
        else
        {
            return ChrBaseRegion.from(anchor.chromosome(), probeRegionStartingAt(anchor.start()));
        }
    }

    private static TargetMetadata createTargetMetadata(final IgTcrGene gene)
    {
        String extraInfo = gene.geneName();
        return new TargetMetadata(TargetMetadata.Type.CDR3, extraInfo);
    }
}
