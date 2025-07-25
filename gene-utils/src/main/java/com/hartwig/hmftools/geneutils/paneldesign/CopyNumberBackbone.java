package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.String.format;
import static java.util.Collections.emptyList;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates.refGenomeCoordinates;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_BACKBONE_CENTROMERE_MARGIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_BACKBONE_QUALITY_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_BACKBONE_GNOMAD_FREQ_MAX;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_BACKBONE_GNOMAD_FREQ_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_BACKBONE_PARTITION_SIZE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_GC_OPTIMAL_TOLERANCE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_GC_TARGET;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_GC_TOLERANCE;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeSelector.selectBestProbe;

import java.util.Arrays;
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
import org.jetbrains.annotations.Nullable;

// Probes based on Amber heterozygous sites, used to deduce copy number.
// Methodology:
//   - Divide chromosomes into large partitions;
//   - In each partition, generate candidate probes on each Amber site;
//   - In each partition, select the acceptable probe with the best GC content.
public class CopyNumberBackbone
{
    private final String mAmberSitesFile;
    private final RefGenomeVersion mRefGenomeVersion;
    private final ProbeGenerator mProbeGenerator;
    private final PanelData mPanelData;

    private static final TargetMetadata.Type TARGET_REGION_TYPE = TargetMetadata.Type.CN_BACKBONE;

    private static final ProbeSelectCriteria PROBE_CRITERIA = new ProbeSelectCriteria(
            new ProbeEvaluator.Criteria(CN_BACKBONE_QUALITY_MIN, CN_GC_TARGET, CN_GC_TOLERANCE),
            new ProbeSelector.Strategy.BestGc(CN_GC_OPTIMAL_TOLERANCE));

    private static final Logger LOGGER = LogManager.getLogger(CopyNumberBackbone.class);

    public CopyNumberBackbone(final String amberSitesFile, final RefGenomeVersion refGenomeVersion, final ProbeGenerator probeGenerator,
            PanelData panelData)
    {
        mAmberSitesFile = amberSitesFile;
        mRefGenomeVersion = refGenomeVersion;
        mProbeGenerator = probeGenerator;
        mPanelData = panelData;
    }

    public void generateProbes()
    {
        LOGGER.info("Generating copy number backbone probes");

        Map<String, List<Partition>> partitions = createPartitions();

        populateAmberSites(partitions);

        ProbeGenerationResult result = generateProbes(partitions);
        // Probes generated here cannot overlap with themselves since there is 1 probe per partition.
        // So it's safe to generate all the probes together and then add them to the result at the end.
        // No need to check overlap of generated probes with themselves.
        mPanelData.addResult(result);

        LOGGER.info("Done generating copy number backbone probes");
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

    private Map<String, List<Partition>> createPartitions()
    {
        RefGenomeCoordinates refGenomeCoordinates = refGenomeCoordinates(mRefGenomeVersion);

        return Arrays.stream(HumanChromosome.values()).map(chromosome ->
        {
            // Partitions spaced across the chromosome, avoiding the centromeres.

            String chrStr = mRefGenomeVersion.versionedChromosome(chromosome.toString());
            int centromere = refGenomeCoordinates.centromeres().get(chromosome);
            int centromereMin = centromere - CN_BACKBONE_CENTROMERE_MARGIN;
            int centromereMax = centromere + CN_BACKBONE_CENTROMERE_MARGIN;
            LOGGER.debug("Excluded centromere region {}:{}-{}", chromosome, centromereMin, centromereMax);

            List<Partition> chrPartitions = partitionChromosome(chrStr, mRefGenomeVersion, CN_BACKBONE_PARTITION_SIZE).stream()
                    .filter(region -> !region.overlaps(chrStr, centromereMin, centromereMax))
                    .map(Partition::new)
                    .toList();

            chrPartitions.forEach(partition -> LOGGER.trace("Copy number backbone partition: {}", partition.Region));

            return Map.entry(chrStr, chrPartitions);
        }).collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
    }

    private void populateAmberSites(Map<String, List<Partition>> partitions)
    {
        enum SiteFilter
        {
            CANDIDATE,
            EXCLUDED_REGION,
            GNOMAD_FREQ
        }

        List<AmberSite> sites = AmberSites.loadAmberSitesFile(mAmberSitesFile);
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

        String filterCountStr = Arrays.stream(SiteFilter.values())
                .map(x -> format("%s=%d", x.toString(), siteFilterCounts[x.ordinal()])).collect(Collectors.joining(", "));
        LOGGER.debug("Amber site filters: {}", filterCountStr);
    }

    private ProbeGenerationResult generateProbes(final Map<String, List<Partition>> partitions)
    {
        return partitions.values().stream()
                .flatMap(chrPartitions ->
                        chrPartitions.stream().map(this::generateProbe))
                .reduce(new ProbeGenerationResult(), ProbeGenerationResult::add);
    }

    private ProbeGenerationResult generateProbe(final Partition partition)
    {
        LOGGER.trace("Generating probes for {} with {} Amber sites", partition.Region, partition.Sites.size());

        Stream<Probe> candidates = generateCandidateProbes(partition);
        Stream<Probe> evaluatedCandidates = mProbeGenerator.mProbeEvaluator.evaluateProbes(candidates, PROBE_CRITERIA.eval());
        Optional<Probe> bestCandidate = selectBestProbe(evaluatedCandidates, PROBE_CRITERIA.select());
        LOGGER.trace("{}: Best probe: {}", partition.Region, bestCandidate);

        ProbeGenerationResult result = bestCandidate
                .map(bestProbe ->
                {
                    TargetMetadataExtra metadataExtra = (TargetMetadataExtra) requireNonNull(bestProbe.metadata().extraData());
                    AmberSite site = requireNonNull(metadataExtra.site());
                    TargetRegion target = new TargetRegion(ChrBaseRegion.from(site.position()), bestProbe.metadata());
                    if(mPanelData.isCovered(target.region()))
                    {
                        LOGGER.debug("Copy number backbone target already covered by panel: {}", target);
                        return new ProbeGenerationResult(emptyList(), List.of(target), emptyList(), emptyList());
                    }
                    else
                    {
                        return new ProbeGenerationResult(List.of(bestProbe), List.of(target), List.of(target), emptyList());
                    }
                })
                .orElseGet(() ->
                {
                    LOGGER.debug("No acceptable probe for copy number backbone partition: {}", partition.Region);

                    TargetMetadata metadata = createTargetMetadata(partition.Region, null);
                    TargetRegion target = new TargetRegion(partition.Region, metadata);

                    String rejectionReason;
                    if(partition.Sites.isEmpty())
                    {
                        rejectionReason = "No Amber sites in partition";
                    }
                    else
                    {
                        rejectionReason = "No probe covering Amber sites meets criteria " + PROBE_CRITERIA.eval();
                    }

                    return ProbeGenerationResult.rejectTarget(target, rejectionReason);
                });
        return result;
    }

    private Stream<Probe> generateCandidateProbes(final Partition partition)
    {
        return partition.Sites.stream()
                .flatMap(site -> generateCandidateProbes(partition, site))
                // Plausible that a site is near the edge of a partition such that the probe goes outside the partition and/or chromosome.
                .filter(probe -> partition.Region.containsRegion(probe.region()));
    }

    private Stream<Probe> generateCandidateProbes(final Partition partition, final AmberSite site)
    {
        TargetMetadata metadata = createTargetMetadata(partition.Region, site);
        ProbeContext context = new ProbeContext(metadata);
        return mProbeGenerator.mCandidateGenerator.coverPosition(site.position(), context);
    }

    public record TargetMetadataExtra(
            @Nullable AmberSite site
    ) implements TargetMetadata.ExtraData
    {
    }

    private static TargetMetadata createTargetMetadata(final ChrBaseRegion partitionRegion, @Nullable final AmberSite site)
    {
        String extraInfo = partitionRegion.toString();
        if(site != null)
        {
            extraInfo += format(":%d", site.position().Position);
        }
        TargetMetadataExtra extraData = new TargetMetadataExtra(site);
        return new TargetMetadata(TARGET_REGION_TYPE, extraInfo, extraData);
    }
}
