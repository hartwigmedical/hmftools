package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.String.format;
import static java.util.Collections.emptyList;

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

// TODO: unit test

// Probes based on Amber heterozygous sites, used to deduce copy number.
// Methodology:
//   - Divide chromosomes into large partitions;
//   - In each partition, generate candidate probes on each Amber site;
//   - In each partition, select the acceptable probe with the best GC content.
// TODO: this should probably be a predefined list rather than generate every time
public class CopyNumberBackbone
{
    private final String mAmberSitesFile;
    private final RefGenomeVersion mRefGenomeVersion;
    private final ProbeGenerator mProbeGenerator;
    private final PanelData mPanelData;

    private static final TargetMetadata.Type TARGET_REGION_TYPE = TargetMetadata.Type.CN_BACKBONE;

    private static final ProbeSelector.Criteria PROBE_CRITERIA = new ProbeSelector.Criteria(
            new ProbeEvaluator.Criteria(CN_BACKBONE_QUALITY_MIN, CN_GC_TARGET, CN_GC_TOLERANCE),
            new ProbeSelector.Strategy.BestGc(CN_GC_OPTIMAL_TOLERANCE));

    private static final Logger LOGGER = LogManager.getLogger(CopyNumberBackbone.class);

    public CopyNumberBackbone(final String amberSitesFile, final RefGenomeVersion refGenomeVersion, final ProbeGenerator probeGenerator,
            final PanelData panelData)
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

        Stream<CandidateProbe> candidates = generateCandidateProbes(partition);
        Optional<EvaluatedProbe> bestCandidate = mProbeGenerator.mProbeSelector.selectBestCandidate(candidates, PROBE_CRITERIA);
        LOGGER.trace("{}: Best probe: {}", partition.Region, bestCandidate);

        ProbeGenerationResult result = bestCandidate
                .map(bestProbe ->
                {
                    TargetRegion target = bestProbe.candidate().target();
                    if(mPanelData.isCovered(target.region()))
                    {
                        LOGGER.debug("Copy number backbone target already covered by panel: {}", target);
                        return new ProbeGenerationResult(List.of(target), emptyList(), emptyList());
                    }
                    else
                    {
                        // TODO: is this the best target region to use here?
                        return new ProbeGenerationResult(List.of(target), List.of(bestProbe), emptyList());
                    }
                })
                .orElseGet(() ->
                {
                    LOGGER.debug("No acceptable probe for copy number backbone partition: {}", partition.Region);

                    TargetMetadata metadata = createTargetMetadata(partition.Region, null);
                    // TODO: is this the best target region to use here?
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
                    RejectedRegion rejectedRegion = new RejectedRegion(partition.Region, target, rejectionReason);

                    return new ProbeGenerationResult(List.of(target), emptyList(), List.of(rejectedRegion));
                });
        return result;
    }

    private Stream<CandidateProbe> generateCandidateProbes(final Partition partition)
    {
        return partition.Sites.stream()
                .flatMap(site -> generateCandidateProbes(partition, site))
                // Plausible that a site is near the edge of a partition such that the probe goes outside the partition and/or chromosome.
                .filter(probe -> partition.Region.containsRegion(probe.probeRegion()));
    }

    private Stream<CandidateProbe> generateCandidateProbes(final Partition partition, final AmberSite site)
    {
        TargetRegion target = new TargetRegion(ChrBaseRegion.from(site.position()), createTargetMetadata(partition.Region, site));
        CandidateProbeContext candidateContext = new CandidateProbeContext(target);
        return mProbeGenerator.mCandidateGenerator.coverPosition(site.position(), candidateContext);
    }

    private static TargetMetadata createTargetMetadata(final ChrBaseRegion partitionRegion, @Nullable final AmberSite site)
    {
        String extraInfo = partitionRegion.toString();
        if(site != null)
        {
            extraInfo += format(":%d", site.position().Position);
        }
        return new TargetMetadata(TARGET_REGION_TYPE, extraInfo);
    }
}
