package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates.refGenomeCoordinates;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_BACKBONE_CENTROMERE_MARGIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_BACKBONE_QUALITY_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_BACKBONE_GNOMAD_FREQ_MAX;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_BACKBONE_GNOMAD_FREQ_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_BACKBONE_PARTITION_SIZE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_GC_TARGET;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_GC_TOLERANCE;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Probes based on Amber heterozygous sites, used to deduce copy number.
// Methodology:
//   - Divide chromosomes into large partitions;
//   - In each partition, generate candidate probes on each Amber site;
//   - In each partition, select the acceptable probe with the best GC content.
public class CopyNumberBackbone
{
    private static final ProbeSourceType PROBE_SOURCE = ProbeSourceType.CN_BACKBONE;

    private static final ProbeSelectCriteria PROBE_SELECT_CRITERIA = new ProbeSelectCriteria(
            new ProbeEvalCriteria(CN_BACKBONE_QUALITY_MIN, CN_GC_TARGET, CN_GC_TOLERANCE),
            ProbeSelectStrategy.BEST_GC);

    private static final Logger LOGGER = LogManager.getLogger(CopyNumberBackbone.class);

    public static ProbeGenerationResult generateProbes(final String amberSitesFile, final RefGenomeVersion refGenVersion,
            final ProbeGenerator probeGenerator)
    {
        LOGGER.info("Generating copy number backbone probes");

        Map<String, List<Partition>> partitions = createPartitions(refGenVersion);

        populateAmberSites(partitions, amberSitesFile);

        ProbeGenerationResult result = generateProbes(partitions, probeGenerator);

        LOGGER.info("Done generating copy number backbone probes");
        return result;
    }

    private static class Partition
    {
        public final ChrBaseRegion Region;
        public final List<AmberSite> Sites;

        public Partition(final ChrBaseRegion region)
        {
            Region = region;
            Sites = Lists.newArrayList();
        }
    }

    private static Map<String, List<Partition>> createPartitions(final RefGenomeVersion refGenomeVersion)
    {
        RefGenomeCoordinates refGenomeCoordinates = refGenomeCoordinates(refGenomeVersion);

        return Arrays.stream(HumanChromosome.values()).map(chromosome ->
        {
            // Partitions spaced across the chromosome, avoiding the centromeres.

            String chrStr = refGenomeVersion.versionedChromosome(chromosome.toString());
            int centromere = refGenomeCoordinates.centromeres().get(chromosome);
            int centromereMin = centromere - CN_BACKBONE_CENTROMERE_MARGIN;
            int centromereMax = centromere + CN_BACKBONE_CENTROMERE_MARGIN;

            List<Partition> chrPartitions = partitionChromosome(chrStr, refGenomeVersion, CN_BACKBONE_PARTITION_SIZE).stream()
                    .filter(region -> !region.overlaps(chrStr, centromereMin, centromereMax))
                    .map(Partition::new)
                    .toList();

            chrPartitions.forEach(partition -> LOGGER.trace("Copy number backbone partition: {}", partition.Region));

            return Map.entry(chrStr, chrPartitions);
        }).collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
    }

    private static void populateAmberSites(Map<String, List<Partition>> partitions, final String sitesFilePath)
    {
        enum SiteFilter
        {
            CANDIDATE,
            EXCLUDED_REGION,
            GNOMAD_FREQ
        }

        List<AmberSite> sites = AmberSites.loadAmberSitesFile(sitesFilePath);
        int[] siteFilterCounts = new int[SiteFilter.values().length];
        for(AmberSite site : sites)
        {
            if(!(site.gnomadFreq() >= CN_BACKBONE_GNOMAD_FREQ_MIN && site.gnomadFreq() <= CN_BACKBONE_GNOMAD_FREQ_MAX))
            {
                ++siteFilterCounts[SiteFilter.GNOMAD_FREQ.ordinal()];
                continue;
            }
            else if(!HumanChromosome.contains(site.position().Chromosome))
            {
                ++siteFilterCounts[SiteFilter.EXCLUDED_REGION.ordinal()];
                continue;
            }
            List<Partition> chrPartitions = partitions.get(site.position().Chromosome);
            if(chrPartitions == null)
            {
                ++siteFilterCounts[SiteFilter.EXCLUDED_REGION.ordinal()];
                continue;
            }
            Partition partition = chrPartitions.stream()
                    .filter(p -> p.Region.containsPosition(site.position().Position)).findFirst().orElse(null);
            if(partition == null)
            {
                ++siteFilterCounts[SiteFilter.EXCLUDED_REGION.ordinal()];
                continue;
            }
            ++siteFilterCounts[SiteFilter.CANDIDATE.ordinal()];
            partition.Sites.add(site);
        }

        LOGGER.info("Loaded {} Amber sites from {}", sites.size(), sitesFilePath);

        String filterCountStr = Arrays.stream(SiteFilter.values())
                .map(x -> format("%s=%d", x.toString(), siteFilterCounts[x.ordinal()])).collect(Collectors.joining(", "));
        LOGGER.debug("Amber site filters: {}", filterCountStr);
    }

    private static ProbeGenerationResult generateProbes(final Map<String, List<Partition>> partitions, final ProbeGenerator probeGenerator)
    {
        return partitions.values().stream()
                .flatMap(chrPartitions ->
                        chrPartitions.stream().map(partition -> generateProbe(partition, probeGenerator)))
                .reduce(new ProbeGenerationResult(), ProbeGenerationResult::add);
    }

    private static ProbeGenerationResult generateProbe(final Partition partition, final ProbeGenerator probeGenerator)
    {
        LOGGER.trace("Generating probes for {} with {} Amber sites", partition.Region, partition.Sites.size());

        Optional<EvaluatedProbe> bestCandidate = probeGenerator.selectBestProbe(generateCandidateProbes(partition), PROBE_SELECT_CRITERIA);
        LOGGER.trace("{}: Best probe: {}", partition.Region, bestCandidate);

        ProbeSourceInfo source = new ProbeSourceInfo(PROBE_SOURCE, partition.Region.toString());
        TargetRegion target = new TargetRegion(source, partition.Region);

        ProbeGenerationResult result = bestCandidate
                .map(bestProbe -> new ProbeGenerationResult(List.of(target), List.of(bestProbe), Collections.emptyList()))
                .orElseGet(() ->
                {
                    // Given the Amber sites are predetermined and there's a lot of them, in typical use we should find an acceptable probe.
                    // If no probe is found then the input data is probably wrong (or the code is being tested).
                    LOGGER.warn("No acceptable probe for copy number backbone partition: {}", partition.Region);

                    String rejectionReason;
                    if(partition.Sites.isEmpty())
                    {
                        rejectionReason = "No Amber sites in partition";
                    }
                    else
                    {
                        rejectionReason = "No probe covering Amber sites meets criteria " + PROBE_SELECT_CRITERIA.eval();
                    }

                    RejectedRegion rejectedRegion = new RejectedRegion(partition.Region, source, rejectionReason);
                    return new ProbeGenerationResult(List.of(target), Collections.emptyList(), List.of(rejectedRegion));
                });
        return result;
    }

    private static Stream<CandidateProbe> generateCandidateProbes(final Partition partition)
    {
        return partition.Sites.stream()
                .flatMap(CopyNumberBackbone::generateCandidateProbes)
                // Plausible that a site is near the edge of a partition such that the probe goes outside the partition and/or chromosome.
                .filter(probe -> partition.Region.containsRegion(probe.probeRegion()));
    }

    private static Stream<CandidateProbe> generateCandidateProbes(final AmberSite site)
    {
        return ProbeGenerator.coverPositionCandidates(site.position(), createProbeSourceInfo(site));
    }

    private static ProbeSourceInfo createProbeSourceInfo(final AmberSite site)
    {
        String extraInfo = site.position().toString();
        return new ProbeSourceInfo(PROBE_SOURCE, extraInfo);
    }
}
