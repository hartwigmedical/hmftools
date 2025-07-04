package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates.refGenomeCoordinates;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.CN_BACKBONE_CENTROMERE_MARGIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.CN_BACKBONE_GC_RATIO_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.CN_BACKBONE_GNMOD_FREQ_MAX;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.CN_BACKBONE_GNMOD_FREQ_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.CN_BACKBONE_MAPPABILITY;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.CN_BACKBONE_PARTITION_SIZE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.MIN_PROBE_QUALITY_SCORE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeCandidate.createProbeCandidate;

import java.io.BufferedReader;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class CopyNumberBackbone
{
    private final PanelConfig mConfig;
    private final PanelCache mPanelCache;
    private final ProbeQualityProfile mProbeQualityProfile;

    private final Map<String, List<Partition>> mChrPartitionsMap;

    public CopyNumberBackbone(final PanelConfig config, final PanelCache panelCache, final ProbeQualityProfile probeQualityProfile)
    {
        mConfig = config;
        mPanelCache = panelCache;
        mProbeQualityProfile = probeQualityProfile;

        mChrPartitionsMap = Maps.newHashMap();

        createPartitions();
    }

    public void run()
    {
        // evaluate sites for any partition which doesn't already have a region or probe
        markCoveredPartitions();

        loadSites();

        selectSiteProbes();
    }

    private void markCoveredPartitions()
    {
        for(Map.Entry<String, List<Partition>> entry : mChrPartitionsMap.entrySet())
        {
            String chromosome = entry.getKey();
            List<Partition> partitions = entry.getValue();

            List<PanelRegion> panelRegions = mPanelCache.chromosomeRegions(chromosome);

            for(Partition partition : partitions)
            {
                if(panelRegions.stream().anyMatch(partition.Region::overlaps))
                {
                    partition.HasExistingProbes = true;
                }
            }
        }
    }

    private void selectSiteProbes()
    {
        List<AmberSite> amberSites = Lists.newArrayList();
        List<ChrBaseRegion> regions = Lists.newArrayList();

        for(List<Partition> partitions : mChrPartitionsMap.values())
        {
            for(Partition partition : partitions)
            {
                if(partition.HasExistingProbes)
                {
                    continue;
                }

                for(AmberSite amberSite : partition.Sites)
                {
                    int probeStart = max(1, amberSite.Position - PROBE_LENGTH / 2);
                    int probeEnd = probeStart + PROBE_LENGTH - 1;

                    ProbeCandidate probe = createProbeCandidate(
                            new ChrBaseRegion(amberSite.Chromosome, probeStart, probeEnd), mConfig.RefGenome);

                    amberSite.setProbe(probe);

                    amberSites.add(amberSite);
                    regions.add(probe.region());
                }
            }
        }

        if(amberSites.isEmpty())
        {
            GU_LOGGER.info("no valid Amber sites found");
            return;
        }

        GU_LOGGER.info("Computing quality scores for {} Amber site probes", regions.size());

        List<Optional<Double>> qualityScores = regions.stream().map(mProbeQualityProfile::computeQualityScore).toList();

        for(int i = 0; i < amberSites.size(); ++i)
        {
            AmberSite amberSite = amberSites.get(i);
            ProbeCandidate probe = amberSite.probe();
            qualityScores.get(i).ifPresent(probe::setQualityScore);
        }

        // take the lowest scoring site for each partition
        for(List<Partition> partitions : mChrPartitionsMap.values())
        {
            for(Partition partition : partitions)
            {
                if(partition.HasExistingProbes)
                {
                    continue;
                }

                double bestScore = MIN_PROBE_QUALITY_SCORE;
                AmberSite topAmberSite = null;

                for(AmberSite amberSite : partition.Sites)
                {
                    ProbeCandidate probe = amberSite.probe();

                    if(probe.getQualityScore().isEmpty())
                    {
                        continue;
                    }

                    if(topAmberSite == null || probe.getQualityScore().get() > bestScore)
                    {
                        bestScore = probe.getQualityScore().get();
                        topAmberSite = amberSite;
                    }
                }

                if(topAmberSite != null)
                {
                    ProbeCandidate probe = topAmberSite.probe();

                    PanelRegion amberProbe = new PanelRegion(
                            probe.region(), RegionType.CN_BACKBONE, format("%s:%d", topAmberSite.Chromosome, topAmberSite.Position),
                            probe.getSequence(), probe.getGcContent(), probe.getQualityScore().get());

                    mPanelCache.addRegion(amberProbe);
                }
            }
        }
    }

    private static class Partition
    {
        public final BaseRegion Region;
        public final List<AmberSite> Sites;
        public boolean HasExistingProbes;

        public Partition(final BaseRegion region)
        {
            Region = region;
            Sites = Lists.newArrayList();
            HasExistingProbes = false;
        }

        public String toString()
        {
            return format("region(%s) sites(%d) hasExistingProbes(%s)", Region, Sites.size(), HasExistingProbes);
        }
    }

    private void createPartitions()
    {
        RefGenomeCoordinates refGenomeCoordinates = refGenomeCoordinates(mConfig.RefGenVersion);

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            int centromere = refGenomeCoordinates.centromeres().get(chromosome);
            int centromereMin = centromere - CN_BACKBONE_CENTROMERE_MARGIN;
            int centromereMax = centromere + CN_BACKBONE_CENTROMERE_MARGIN;

            List<Partition> chrPartitions = Lists.newArrayList();
            mChrPartitionsMap.put(chrStr, chrPartitions);

            List<ChrBaseRegion> regions = partitionChromosome(chrStr, mConfig.RefGenVersion, CN_BACKBONE_PARTITION_SIZE);

            for(ChrBaseRegion region : regions)
            {
                if(region.overlaps(chrStr, centromereMin, centromereMax))
                {
                    continue;
                }

                chrPartitions.add(new Partition(new BaseRegion(region.start(), region.end())));
            }
        }
    }

    private static final String FLD_GNOMAD_FREQ = "GnomadFreq";
    private static final String FLD_MAPPABILITY = "Mappability";
    private static final String FLD_GC_RATIO = "GcRatio";

    private enum SiteFilter
    {
        REQUIRED,
        EXCLUDED,
        GNOMAD,
        GC_RATIO,
        MAPPABILITY,
        POPULATED;
    }

    public void loadSites()
    {
        try(BufferedReader reader = createBufferedReader(mConfig.AmberSitesFile))
        {
            String header = reader.readLine();

            Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posIndex = fieldsIndexMap.get(FLD_POSITION);
            int mapIndex = fieldsIndexMap.get(FLD_MAPPABILITY);
            int gnomadIndex = fieldsIndexMap.get(FLD_GNOMAD_FREQ);
            int gcRatioIndex = fieldsIndexMap.get(FLD_GC_RATIO);

            int[] siteFilterCounts = new int[SiteFilter.values().length];

            String line = null;
            List<Partition> chrPartitions = null;
            Partition currentPartition = null;
            String currentChromosome = "";
            int siteCount = 0;

            while((line = reader.readLine()) != null)
            {
                String[] values = line.split(TSV_DELIM, -1);

                ++siteCount;

                double mappability = Double.parseDouble(values[mapIndex]);

                if(mappability < CN_BACKBONE_MAPPABILITY)
                {
                    ++siteFilterCounts[SiteFilter.MAPPABILITY.ordinal()];
                    continue;
                }

                double gnomadFreq = Double.parseDouble(values[gnomadIndex]);

                if(gnomadFreq < CN_BACKBONE_GNMOD_FREQ_MIN || gnomadFreq > CN_BACKBONE_GNMOD_FREQ_MAX)
                {
                    ++siteFilterCounts[SiteFilter.GNOMAD.ordinal()];
                    continue;
                }

                double gcRatio = Double.parseDouble(values[gcRatioIndex]);

                if(gcRatio < CN_BACKBONE_GC_RATIO_MIN)
                {
                    ++siteFilterCounts[SiteFilter.GC_RATIO.ordinal()];
                    continue;
                }

                String chrStr = values[chrIndex];

                if(!HumanChromosome.contains(chrStr))
                {
                    continue;
                }

                if(!currentChromosome.equals(chrStr))
                {
                    currentChromosome = chrStr;
                    chrPartitions = mChrPartitionsMap.get(chrStr);
                    currentPartition = null;
                }

                int position = Integer.parseInt(values[posIndex]);

                if(currentPartition == null || !currentPartition.Region.containsPosition(position))
                {
                    currentPartition = chrPartitions.stream().filter(x -> x.Region.containsPosition(position)).findFirst().orElse(null);

                    if(currentPartition == null)
                    {
                        ++siteFilterCounts[SiteFilter.EXCLUDED.ordinal()];
                        continue;
                    }

                    if(currentPartition.HasExistingProbes)
                    {
                        ++siteFilterCounts[SiteFilter.POPULATED.ordinal()];
                        continue;
                    }
                }

                currentPartition.Sites.add(new AmberSite(chrStr, position));
                ++siteFilterCounts[SiteFilter.REQUIRED.ordinal()];
            }

            GU_LOGGER.info("loaded {} Amber sites from {}", siteCount, mConfig.AmberSitesFile);

            String filterCountStr = Arrays.stream(SiteFilter.values())
                    .map(x -> format("%s=%d", x.toString(), siteFilterCounts[x.ordinal()])).collect(Collectors.joining(", "));
            GU_LOGGER.debug("Amber site filters: {}", filterCountStr);
        }
        catch(Exception e)
        {
            GU_LOGGER.error("failed to load Amber sites file: {}", e.toString());
        }
    }
}
