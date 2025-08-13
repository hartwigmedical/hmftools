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
    private static final TargetMetadata.Type TARGET_TYPE = TargetMetadata.Type.CDR3;

    private static final ProbeEvaluator.Criteria PROBE_CRITERIA = new ProbeEvaluator.Criteria(
            CDR3_QUALITY_MIN, CDR3_GC_TARGET, CDR3_GC_TOLERANCE);
    private static final ProbeSelector.Strategy PROBE_SELECT = new ProbeSelector.Strategy.FirstAcceptable();

    private static final Logger LOGGER = LogManager.getLogger(Cdr3Regions.class);

    public static void generateProbes(final RefGenomeVersion refGenomeVersion, final ProbeGenerator probeGenerator, PanelData panelData)
    {
        LOGGER.info("Generating CDR3 probes");

        List<IgTcrGene> genes = loadIgTcrGenes(refGenomeVersion);

        ProbeGenerationResult result = generateProbes(genes, probeGenerator, panelData);
        // Overlaps between CDR3 probes are checked as they are generated, so it's safe to add them to the result all at once at the end.
        panelData.addResult(result);

        LOGGER.info("Done generating CDR3 probes");
    }

    private static List<IgTcrGene> loadIgTcrGenes(final RefGenomeVersion refGenomeVersion)
    {
        List<IgTcrGene> genes = IgTcrGeneFile.read(refGenomeVersion).stream()
                .filter(gene -> gene.region() == IgTcrRegion.V_REGION || gene.region() == IgTcrRegion.J_REGION)
                .filter(IgTcrGene::inPrimaryAssembly)
                .filter(gene -> gene.anchorLocation() != null)
                .toList();
        LOGGER.info("Loaded {} V/J genes", genes.size());
        return genes;
    }

    private static ProbeGenerationResult generateProbes(final List<IgTcrGene> genes, final ProbeGenerator probeGenerator,
            final PanelCoverage coverage)
    {
        ProbeGenerationResult result = new ProbeGenerationResult();
        List<ChrBaseRegion> coveredRegions = new ArrayList<>();
        for(IgTcrGene gene : genes)
        {
            result = result.add(generateProbe(gene, probeGenerator, coverage, coveredRegions));
        }
        return result;
    }

    private static ProbeGenerationResult generateProbe(final IgTcrGene gene, final ProbeGenerator probeGenerator,
            final PanelCoverage coverage, List<ChrBaseRegion> coveredRegions)
    {
        ChrBaseRegion targetRegion = calculateTargetRegion(gene);

        if(coveredRegions.stream().anyMatch(targetRegion::overlaps))
        {
            // It's possible regions overlap other regions, in which case just take the first and discard the rest.
            LOGGER.trace("CDR3 region overlaps with another; discarding");
            return new ProbeGenerationResult();
        }
        else
        {
            TargetMetadata metadata = createTargetMetadata(gene);
            // TODO: need to prevent the probe from going past the end of the gene?
            // TODO: this can produce more than 1 probe due to acceptable subregion splitting. ok?
            ProbeGenerationResult result = probeGenerator.coverRegion(targetRegion, metadata, PROBE_CRITERIA, PROBE_SELECT, coverage);
            if(!result.probes().isEmpty())
            {
                coveredRegions.add(targetRegion);
            }
            return result;
        }
    }

    private static ChrBaseRegion calculateTargetRegion(final IgTcrGene gene)
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
        return new TargetMetadata(TARGET_TYPE, extraInfo);
    }
}
