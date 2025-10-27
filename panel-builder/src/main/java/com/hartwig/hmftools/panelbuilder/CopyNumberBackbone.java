package com.hartwig.hmftools.panelbuilder;

import static java.lang.String.format;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates.refGenomeCoordinates;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CN_BACKBONE_ALTERNATE_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CN_BACKBONE_CENTROMERE_MARGIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CN_BACKBONE_GNOMAD_FREQ_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CN_BACKBONE_GNOMAD_FREQ_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CN_BACKBONE_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CN_GC_OPTIMAL_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CN_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CN_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Probes based on Amber heterozygous sites, used to calculate accurate copy number.
// Methodology:
//   - Divide chromosomes into large partitions;
//   - In each partition, generate candidate probes on each Amber site.
//      Or if on chrY and there are no Amber sites, consider all probe positions in the partition.
//   - In each partition, select the acceptable probe with the best GC content.
public class CopyNumberBackbone
{
    private final String mAmberSitesFile;
    private final int mResolution;
    private final RefGenomeVersion mRefGenomeVersion;
    private final ProbeGenerator mProbeGenerator;
    private final PanelData mPanelData;

    private static final TargetMetadata.Type TARGET_TYPE = TargetMetadata.Type.CN_BACKBONE;

    private static final ProbeEvaluator.Criteria PROBE_CRITERIA = new ProbeEvaluator.Criteria(
            CN_BACKBONE_QUALITY_MIN, CN_GC_TARGET, CN_GC_TOLERANCE);
    private static final ProbeEvaluator.Criteria PROBE_CRITERIA_ALTERNATIVE = new ProbeEvaluator.Criteria(
            CN_BACKBONE_ALTERNATE_QUALITY_MIN, CN_GC_TARGET, CN_GC_TOLERANCE);
    private static final ProbeSelector.Strategy PROBE_SELECT = new ProbeSelector.Strategy.BestGc(CN_GC_OPTIMAL_TOLERANCE);

    private static final String FLD_GNOMAD_FREQ = "GnomadFreq";

    private static final Logger LOGGER = LogManager.getLogger(CopyNumberBackbone.class);

    public CopyNumberBackbone(final String amberSitesFile, int resolution, final RefGenomeVersion refGenomeVersion,
            final ProbeGenerator probeGenerator, PanelData panelData)
    {
        mAmberSitesFile = amberSitesFile;
        mResolution = resolution;
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
        // Probes generated here cannot overlap with themselves since there is one probe per partition.
        // So it's safe to generate all the probes together and then add them to the result at the end.
        // No need to check overlap of generated probes with themselves.
        mPanelData.addResult(result);

        LOGGER.info("Done generating copy number backbone probes");
    }

    private record AmberSite(
            BasePosition position,
            double gnomadFreq
    )
    {
    }

    private static List<AmberSite> loadAmberSitesFile(final String path)
    {
        LOGGER.debug("Loading Amber sites file: {}", path);

        try(DelimFileReader reader = new DelimFileReader(path))
        {
            int chrIndex = requireNonNull(reader.getColumnIndex(FLD_CHROMOSOME));
            int posIndex = requireNonNull(reader.getColumnIndex(FLD_POSITION));
            int gnomadIndex = requireNonNull(reader.getColumnIndex(FLD_GNOMAD_FREQ));

            List<AmberSite> sites = reader.stream().map(row ->
            {
                String chromosome = row.get(chrIndex);
                int position = row.getInt(posIndex);
                double gnomadFreq = row.getDouble(gnomadIndex);
                return new AmberSite(new BasePosition(chromosome, position), gnomadFreq);
            }).toList();

            LOGGER.debug("Loaded {} Amber sites from {}", sites.size(), path);
            return sites;
        }
    }

    private static class Partition
    {
        public final ChrBaseRegion Region;
        public final List<AmberSite> Sites;

        public Partition(final ChrBaseRegion region)
        {
            Region = region;
            Sites = new ArrayList<>();
        }
    }

    private Map<String, List<Partition>> createPartitions()
    {
        LOGGER.debug("Creating copy number backbone partitions with resolution: {}b", mResolution);
        if(mResolution < 1000)
        {
            throw new IllegalArgumentException("Copy number backbone resolution too small");
        }

        RefGenomeCoordinates refGenomeCoordinates = refGenomeCoordinates(mRefGenomeVersion);

        return Arrays.stream(HumanChromosome.values()).map(chromosome ->
        {
            // Partitions spaced across the chromosome, avoiding the centromeres.

            String chrStr = mRefGenomeVersion.versionedChromosome(chromosome);
            int centromere = refGenomeCoordinates.centromeres().get(chromosome);
            int centromereMin = centromere - CN_BACKBONE_CENTROMERE_MARGIN;
            int centromereMax = centromere + CN_BACKBONE_CENTROMERE_MARGIN;
            LOGGER.debug("Excluded centromere region {}:{}-{}", chromosome, centromereMin, centromereMax);

            List<Partition> chrPartitions = partitionChromosome(chrStr, mRefGenomeVersion, mResolution).stream()
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

        List<AmberSite> sites = loadAmberSitesFile(mAmberSitesFile);
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
        return partitions.entrySet().stream()
                // Sort to ensure deterministic ordering.
                .sorted(Map.Entry.comparingByKey())
                .flatMap(entry ->
                        entry.getValue().stream().map(this::generateProbe))
                .reduce(new ProbeGenerationResult(), ProbeGenerationResult::add);
    }

    private ProbeGenerationResult generateProbe(final Partition partition)
    {
        LOGGER.trace("Generating probes for {} with {} Amber sites", partition.Region, partition.Sites.size());

        TargetMetadata metadata = createTargetMetadata(partition.Region);

        // On the Y chromosome there are no heterozygous sites, but we still want copy number determination.
        // The methodology is the pick the first position that meets stricter criteria.
        boolean alternativeMethod = useAlternativeCandidateMethod(partition);
        if(alternativeMethod)
        {
            LOGGER.trace("Using alternative methodology");
            return mProbeGenerator.coverOneSubregion(partition.Region, metadata, PROBE_CRITERIA_ALTERNATIVE, PROBE_SELECT, mPanelData);
        }
        else
        {
            BaseRegion positionBounds = new BaseRegion(partition.Region.start() + PROBE_LENGTH, partition.Region.end() - PROBE_LENGTH);
            Stream<BasePosition> candidatePositions = partition.Sites.stream()
                    .map(AmberSite::position)
                    // Plausible that a site is near the edge of a partition such that the probe goes outside the partition.
                    // Filter these out to avoid probe overlap, since there are many probes to choose from.
                    .filter(position -> positionBounds.containsPosition(position.Position));
            // TODO: should only centre a probe on the site. simply code and more accurate coverage
            return mProbeGenerator.coverOnePosition(candidatePositions, metadata, PROBE_CRITERIA, PROBE_SELECT, mPanelData);
        }
    }

    private static boolean useAlternativeCandidateMethod(final Partition partition)
    {
        return partition.Region.humanChromosome() == HumanChromosome._Y && partition.Sites.isEmpty();
    }

    private static TargetMetadata createTargetMetadata(final ChrBaseRegion partitionRegion)
    {
        String extraInfo = partitionRegion.toString();
        return new TargetMetadata(TARGET_TYPE, extraInfo);
    }
}
