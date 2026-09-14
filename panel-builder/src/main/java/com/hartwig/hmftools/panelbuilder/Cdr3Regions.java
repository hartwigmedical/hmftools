package com.hartwig.hmftools.panelbuilder;

import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CDR3_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CDR3_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CDR3_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CDR3_V_UPSTREAM_PROBE_OFFSET;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionEndingAt;
import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeRegionStartingAt;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
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
// Methodology:
//   - 1 probe at the end of V regions
//   - 1 probe at the start of J regions
//   - for IG V genes, an additional upstream probe (over FR2) to capture somatically hypermutated sequences
public class Cdr3Regions
{
    private static final TargetMetadata.Type TARGET_TYPE = TargetMetadata.Type.CDR3;

    private static final ProbeEvaluator.Criteria PROBE_CRITERIA = new ProbeEvaluator.Criteria(
            CDR3_QUALITY_MIN, CDR3_GC_TARGET, CDR3_GC_TOLERANCE);

    private static final Logger LOGGER = LogManager.getLogger(Cdr3Regions.class);

    public static void generateProbes(final RefGenomeVersion refGenomeVersion, final ProbeGenerator probeGenerator, PanelData panelData)
    {
        LOGGER.info("Generating CDR3 probes");

        List<IgTcrGene> genes = loadIgTcrGenes(refGenomeVersion);
        genes = filterAlleles(genes);

        generateProbes(genes, probeGenerator, panelData);

        LOGGER.info("Done generating CDR3 probes");
    }

    private static List<IgTcrGene> loadIgTcrGenes(final RefGenomeVersion refGenomeVersion)
    {
        List<IgTcrGene> genes = IgTcrGeneFile.read(refGenomeVersion).stream()
                .filter(gene -> gene.region() == IgTcrRegion.V_REGION || gene.region() == IgTcrRegion.J_REGION)
                .filter(IgTcrGene::inPrimaryAssembly)
                .filter(gene -> gene.anchorLocation() != null)
                .toList();
        LOGGER.debug("Loaded {} V/J gene alleles", genes.size());
        return genes;
    }

    private static List<IgTcrGene> filterAlleles(final List<IgTcrGene> genes)
    {
        // No need to compute produce probes for duplicate allele locations.

        Map<String, List<IgTcrGene>> allelesbyGene = new HashMap<>();
        for(IgTcrGene gene : genes)
        {
            allelesbyGene.computeIfAbsent(gene.geneName(), k -> new ArrayList<>()).add(gene);
        }

        List<IgTcrGene> filteredGenes = new ArrayList<>();
        allelesbyGene.forEach((geneName, alleles) ->
        {
            Set<ChrBaseRegion> geneAnchorLocations = new HashSet<>();
            for(IgTcrGene allele : alleles)
            {
                if(geneAnchorLocations.add(allele.anchorLocation()))
                {
                    filteredGenes.add(allele);
                }
            }
        });

        LOGGER.debug("Filtered to {} gene allele locations", filteredGenes.size());

        return filteredGenes;
    }

    private static void generateProbes(final List<IgTcrGene> genes, final ProbeGenerator probeGenerator,
            PanelData panelData)
    {
        Stream<ProbeGenerationSpec> probeGenerationSpecs = genes.stream().flatMap(gene -> createProbeGenerationSpecs(gene).stream());
        probeGenerator.generateBatch(probeGenerationSpecs, panelData);
    }

    private static List<ProbeGenerationSpec> createProbeGenerationSpecs(final IgTcrGene gene)
    {
        List<ProbeGenerationSpec> specs = new ArrayList<>();
        specs.add(createProbeSpec(calculateAnchorProbeRegion(gene), gene.geneName()));
        if(hasUpstreamProbe(gene))
        {
            specs.add(createProbeSpec(calculateUpstreamProbeRegion(gene), gene.geneName() + " upstream"));
        }
        return specs;
    }

    private static ProbeGenerationSpec createProbeSpec(final ChrBaseRegion targetRegion, final String extraInfo)
    {
        // Produce a probe exactly at the determined region or not at all. Shifting probes is not acceptable here.
        // Need to produce a probe which is aligned with the edge of the gene. And want to match the previous CDR3 panel design closely.
        SequenceDefinition sequenceDefinition = SequenceDefinition.singleRegion(targetRegion);
        TargetedRange targetedRange = TargetedRange.wholeRegion(sequenceDefinition.baseLength());
        TargetMetadata metadata = new TargetMetadata(TARGET_TYPE, extraInfo);
        return new ProbeGenerationSpec.SingleProbe(sequenceDefinition, targetedRange, metadata, PROBE_CRITERIA);
    }

    // Probe aligned to the anchor, at the CDR3 edge of the V or J region.
    private static ChrBaseRegion calculateAnchorProbeRegion(final IgTcrGene gene)
    {
        ChrBaseRegion anchor = requireNonNull(gene.anchorLocation());
        if(cdr3SideIsRegionEnd(gene))
        {
            return probeRegionEndingAt(anchor.chromosome(), anchor.end());
        }
        else
        {
            return probeRegionStartingAt(anchor.chromosome(), anchor.start());
        }
    }

    // Probe shifted upstream (toward FR1) from the V anchor, to improve capture for hypermutated V regions.
    // This region is known to experience fewer mutations than the region near the anchor.
    private static ChrBaseRegion calculateUpstreamProbeRegion(final IgTcrGene gene)
    {
        ChrBaseRegion anchor = requireNonNull(gene.anchorLocation());
        if(cdr3SideIsRegionEnd(gene))
        {
            return probeRegionEndingAt(anchor.chromosome(), anchor.end() - CDR3_V_UPSTREAM_PROBE_OFFSET);
        }
        else
        {
            return probeRegionStartingAt(anchor.chromosome(), anchor.start() + CDR3_V_UPSTREAM_PROBE_OFFSET);
        }
    }

    // The anchor sits at the CDR3-facing edge of the region: the high-coordinate end for a forward-strand V or reverse-strand J.
    private static boolean cdr3SideIsRegionEnd(final IgTcrGene gene)
    {
        boolean vForward = gene.region() == IgTcrRegion.V_REGION && gene.geneStrand() == Strand.FORWARD;
        boolean jReverse = gene.region() == IgTcrRegion.J_REGION && gene.geneStrand() == Strand.REVERSE;
        return vForward || jReverse;
    }

    private static boolean hasUpstreamProbe(final IgTcrGene gene)
    {
        return gene.region() == IgTcrRegion.V_REGION
                && gene.geneName().startsWith("IG")
                && !UPSTREAM_PROBE_EXCLUDED_GENES.contains(gene.geneName());
    }

    // IGK-KDE recombination elements, not framework-bearing V/J genes, so an upstream framework probe is meaningless.
    private static final Set<String> UPSTREAM_PROBE_EXCLUDED_GENES = Set.of("IGKINTR", "IGKDEL");
}
