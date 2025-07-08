package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates.refGenomeCoordinates;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_BACKBONE_CENTROMERE_MARGIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_BACKBONE_GC_RATIO_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_BACKBONE_GNMOD_FREQ_MAX;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_BACKBONE_GNMOD_FREQ_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_BACKBONE_MAPPABILITY_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_BACKBONE_PARTITION_SIZE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;

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
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Probes based on Amber heterozygous sites, used to deduce copy number.
// Methodology:
//   - Divide chromosomes into large partitions;
//   - In each partition, generate candidate probes on each Amber site;
//   - Select the best probe within each partition.
public class CopyNumberBackbone
{
    private static final ProbeSource PROBE_SOURCE = ProbeSource.CN_BACKBONE;

    private static final Logger LOGGER = LogManager.getLogger(CopyNumberBackbone.class);

    public static ProbeGenerationResult generateProbes(final String amberSitesFile, final RefGenomeVersion refGenVersion,
            final ProbeEvaluator probeEvaluator)
    {
        Map<String, List<Partition>> partitions = createPartitions(refGenVersion);
        populateAmberSites(partitions, amberSitesFile);
        ProbeGenerationResult result = generateProbes(partitions, probeEvaluator);
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
            return Map.entry(chrStr, chrPartitions);
        }).collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
    }

    private static void populateAmberSites(Map<String, List<Partition>> partitions, final String sitesFilePath)
    {
        enum SiteFilter
        {
            SELECTED,
            EXCLUDED,
            GNOMAD,
            GC_RATIO,
            MAPPABILITY
        }

        List<AmberSite> sites = AmberSites.loadAmberSitesFile(sitesFilePath);
        int[] siteFilterCounts = new int[SiteFilter.values().length];
        for(AmberSite site : sites)
        {
            if(site.mappability() < CN_BACKBONE_MAPPABILITY_MIN)
            {
                ++siteFilterCounts[SiteFilter.MAPPABILITY.ordinal()];
                continue;
            }
            else if(site.gnomadFreq() < CN_BACKBONE_GNMOD_FREQ_MIN || site.gnomadFreq() > CN_BACKBONE_GNMOD_FREQ_MAX)
            {
                ++siteFilterCounts[SiteFilter.GNOMAD.ordinal()];
                continue;
            }
            else if(site.gcRatio() < CN_BACKBONE_GC_RATIO_MIN)
            {
                ++siteFilterCounts[SiteFilter.GC_RATIO.ordinal()];
                continue;
            }
            else if(!HumanChromosome.contains(site.position().Chromosome))
            {
                ++siteFilterCounts[SiteFilter.EXCLUDED.ordinal()];
                continue;
            }
            List<Partition> chrPartitions = partitions.get(site.position().Chromosome);
            if(chrPartitions == null)
            {
                ++siteFilterCounts[SiteFilter.EXCLUDED.ordinal()];
                continue;
            }
            Partition partition = chrPartitions.stream()
                    .filter(p -> p.Region.containsPosition(site.position().Position)).findFirst().orElse(null);
            if(partition == null)
            {
                ++siteFilterCounts[SiteFilter.EXCLUDED.ordinal()];
                continue;
            }
            ++siteFilterCounts[SiteFilter.SELECTED.ordinal()];
            partition.Sites.add(site);
        }

        LOGGER.info("Loaded {} Amber sites from {}", sites.size(), sitesFilePath);

        String filterCountStr = Arrays.stream(SiteFilter.values())
                .map(x -> format("%s=%d", x.toString(), siteFilterCounts[x.ordinal()])).collect(Collectors.joining(", "));
        LOGGER.debug("Amber site filters: {}", filterCountStr);
    }

    private static ProbeGenerationResult generateProbes(final Map<String, List<Partition>> partitions, final ProbeEvaluator probeEvaluator)
    {
        List<ProbeCandidate> probes = partitions.values().stream()
                .flatMap(chrPartitions ->
                        chrPartitions.stream().map(partition -> generateProbe(partition, probeEvaluator)))
                .toList();
        // Note there are never any rejected regions because it's assumed each partition has at least one acceptable probe.
        return new ProbeGenerationResult(probes, Collections.emptyList());
    }

    private static ProbeCandidate generateProbe(final Partition partition, final ProbeEvaluator probeEvaluator)
    {
        Stream<ProbeCandidate> candidates = generateProbeCandidates(partition);
        Optional<ProbeCandidate> bestCandidate = probeEvaluator.selectBestProbe(candidates);
        // Assume there is at least one acceptable probe in each partition to make things simpler.
        // It's likely a valid assumption because the set of Amber sites is fixed and each partition includes many sites.
        return bestCandidate.orElseThrow(
                () -> new RuntimeException(
                        format("Expected at least one Amber site to be an acceptable probe in CN backbone partition %s", partition.Region)));
    }

    private static Stream<ProbeCandidate> generateProbeCandidates(final Partition partition)
    {
        return partition.Sites.stream().map(CopyNumberBackbone::amberSiteToProbe);
    }

    private static ProbeCandidate amberSiteToProbe(final AmberSite site)
    {
        BasePosition position = site.position();
        // Probe centered on the site.
        int probeStart = max(1, position.Position - PROBE_LENGTH / 2);
        int probeEnd = probeStart + PROBE_LENGTH - 1;
        ChrBaseRegion probeRegion = new ChrBaseRegion(position.Chromosome, probeStart, probeEnd);
        ChrBaseRegion targetRegion = new ChrBaseRegion(position.Chromosome, position.Position, position.Position);
        String extraInfo = position.toString();
        ProbeSourceInfo source = new ProbeSourceInfo(PROBE_SOURCE, extraInfo);
        return new ProbeCandidate(source, probeRegion, targetRegion);
    }
}
