package com.hartwig.hmftools.orange.algo.linx;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.driver.DriverType.DRIVERS_LINX_GERMLINE;
import static com.hartwig.hmftools.common.driver.DriverType.DRIVERS_LINX_SOMATIC;
import static com.hartwig.hmftools.common.linx.LinxCommonTypes.SV_VIS_CLUSTER_PREFIX;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.PNG_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getIntValue;
import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCatalogFile;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxCluster;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.linx.LinxGermlineDisruption;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.orange.OrangeConfig;

import org.jetbrains.annotations.Nullable;

public final class LinxDataLoader
{
    public static LinxData load(final OrangeConfig config)throws IOException
    {
        String linxSomaticDir = config.LinxSomaticDataDirectory;
        String linxGermlineDir = config.LinxGermlineDataDirectory;
        String tumorSample = config.TumorId;

        String somaticSvAnnotationFile = LinxSvAnnotation.generateFilename(linxSomaticDir, tumorSample, false);
        String somaticBreakendFile = LinxBreakend.generateFilename(linxSomaticDir, tumorSample);
        String somaticClusterFile = LinxCluster.generateFilename(linxSomaticDir, tumorSample);
        String somaticFusionFile = LinxFusion.generateFilename(linxSomaticDir, tumorSample);
        String somaticDriversFile = LinxDriver.generateFilename(linxSomaticDir, tumorSample);
        String somaticDriverCatalogFile = LinxDriver.generateCatalogFilename(linxSomaticDir, tumorSample, true);

        String germlineBreakendFile = null;
        String germlineDisruptionFile = null;
        String germlineDriverCatalogFile = null;
        String germlineSvAnnotationFile = null;

        if(linxGermlineDir != null)
        {
            germlineSvAnnotationFile = LinxSvAnnotation.generateFilename(linxGermlineDir, tumorSample, true);
            germlineBreakendFile = LinxBreakend.generateFilename(linxGermlineDir, tumorSample, true);
            germlineDisruptionFile = LinxGermlineDisruption.generateFilename(linxGermlineDir, tumorSample);
            germlineDriverCatalogFile = LinxDriver.generateCatalogFilename(linxGermlineDir, tumorSample, false);
        }

        return load(
                somaticSvAnnotationFile, somaticFusionFile, somaticBreakendFile, somaticClusterFile, somaticDriverCatalogFile,
                somaticDriversFile, config.LinxPlotDirectory, germlineSvAnnotationFile, germlineBreakendFile, germlineDisruptionFile,
                germlineDriverCatalogFile);
    }

    private static LinxData load(
            final String somaticSvAnnotationFile, final String fusionFile, final String somaticBreakendFile, final String somaticClusterFile,
            final String somaticDriverCatalogTsv, final String somaticDriverTsv, final String linxPlotDirectory,
            @Nullable String germlineSvAnnotationFile, @Nullable String germlineBreakendFile, @Nullable String germlineDisruptionFile,
            @Nullable String germlineDriverCatalogFile)
            throws IOException
    {
        List<LinxDriver> somaticDriverData = LinxDriver.read(somaticDriverTsv);
        List<LinxFusion> fusions = loadFusions(fusionFile);
        List<LinxBreakend> somaticBreakends = loadBreakends(somaticBreakendFile, fusions);

        List<LinxCluster> somaticClusters = LinxCluster.read(somaticClusterFile);

        List<DriverCatalog> somaticDrivers = DriverCatalogFile.read(somaticDriverCatalogTsv).stream()
                .filter(x -> DRIVERS_LINX_SOMATIC.contains(x.driver())).collect(Collectors.toList());

        // limit SV annotations to reportable breakends
        List<LinxSvAnnotation> somaticSvAnnotations = LinxSvAnnotation.read(somaticSvAnnotationFile);

        Map<LinxBreakend,Integer> somaticBreakendClusterIds = findBreakendClusterId(somaticBreakends, somaticSvAnnotations);

        Map<LinxFusion,Integer> fusionClusterIds = findFusionClusterId(fusions, somaticBreakends, somaticSvAnnotations);

        // likely now redundant since the list will not be persisted
        List<LinxSvAnnotation> somaticBreakendSvAnnotations = Lists.newArrayList(somaticSvAnnotations);
        restrictToMatchingBreakends(somaticBreakendSvAnnotations, somaticBreakends);

        List<DriverCatalog> somaticHomozygousDisruptions = somaticDrivers.stream()
                .filter(x -> x.driver() == DriverType.HOM_DUP_DISRUPTION || x.driver() == DriverType.HOM_DEL_DISRUPTION)
                .collect(Collectors.toList());

        List<LinxSvAnnotation> germlineSvAnnotations = null;
        List<LinxBreakend> reportableGermlineBreakends = null;
        List<DriverCatalog> germlineDrivers = null;
        List<LinxGermlineDisruption> reportableGermlineDisruptions = null;

        if(germlineSvAnnotationFile != null && germlineBreakendFile != null && germlineDisruptionFile != null)
        {
            germlineDrivers = DriverCatalogFile.read(germlineDriverCatalogFile).stream()
                    .filter(x -> DRIVERS_LINX_GERMLINE.contains(x.driver())).collect(Collectors.toList());

            germlineSvAnnotations = LinxSvAnnotation.read(germlineSvAnnotationFile);
            reportableGermlineBreakends = loadBreakends(germlineBreakendFile, Collections.emptyList());

            restrictToMatchingBreakends(germlineSvAnnotations, reportableGermlineBreakends);

            List<LinxGermlineDisruption> allGermlineDisruptions = LinxGermlineDisruption.read(germlineDisruptionFile);
            reportableGermlineDisruptions = selectGermlineSvs(allGermlineDisruptions, reportableGermlineBreakends);
        }

        List<String> reportableEventPlots = Lists.newArrayList();
        List<String> reportableEventPlotTitles = Lists.newArrayList();

        if(linxPlotDirectory != null)
        {
            for(String file : new File(linxPlotDirectory).list())
            {
                if(file.contains(SV_VIS_CLUSTER_PREFIX) && file.endsWith(PNG_EXTENSION))
                {
                    reportableEventPlots.add(linxPlotDirectory + file);
                    reportableEventPlotTitles.add(formPlotTitle(file, somaticSvAnnotations, somaticBreakends, somaticClusters));
                }
            }

            LOGGER.debug(" loaded {} Linx plots from {}", reportableEventPlots.size(), linxPlotDirectory);
        }

        return ImmutableLinxData.builder()
                .somaticSvAnnotations(somaticBreakendSvAnnotations)
                .somaticDrivers(somaticDrivers)
                .somaticDriverData(somaticDriverData)
                .fusions(fusions)
                .fusionClusterIds(fusionClusterIds)
                .somaticBreakends(somaticBreakends)
                .somaticBreakendClusterIds(somaticBreakendClusterIds)
                .somaticHomozygousDisruptions(somaticHomozygousDisruptions)
                .reportableEventPlots(reportableEventPlots)
                .reportableEventPlotTitles(reportableEventPlotTitles)
                .germlineDrivers(germlineDrivers)
                .germlineSvAnnotations(germlineSvAnnotations)
                .germlineBreakends(reportableGermlineBreakends)
                .germlineDisruptions(reportableGermlineDisruptions)
                .build();
    }

    private static List<LinxFusion> loadFusions(final String fusionTsv) throws IOException
    {
        List<LinxFusion> fusions = LinxFusion.read(fusionTsv);

        List<LinxFusion> reportableFusions = new ArrayList<>();
        for(LinxFusion fusion : fusions)
        {
            if(fusion.reported())
            {
                reportableFusions.add(fusion);
            }
        }
        return reportableFusions;
    }

    private static List<LinxBreakend> loadBreakends(final String somaticBreakendTsv, final List<LinxFusion> fusions) throws IOException
    {
        List<LinxBreakend> breakends = LinxBreakend.read(somaticBreakendTsv);

        List<LinxBreakend> reportableBreakends = new ArrayList<>();
        for(LinxBreakend breakend : breakends)
        {
            if(breakend.reportedStatus() == ReportedStatus.REPORTED)
            {
                reportableBreakends.add(breakend);
            }
            else if(fusions.stream().anyMatch(x -> x.fivePrimeBreakendId() == breakend.id() || x.threePrimeBreakendId() == breakend.id()))
            {
                reportableBreakends.add(breakend);
            }
        }

        return reportableBreakends;
    }

    private static Map<LinxBreakend,Integer> findBreakendClusterId(
            final List<LinxBreakend> breakends, final List<LinxSvAnnotation> svAnnotations)
    {
        Map<LinxBreakend,Integer> breakendClusterIds = Maps.newHashMap();

        for(LinxBreakend breakend : breakends)
        {
            Integer clusterId = findBreakendClusterId(breakend, svAnnotations);

            if(clusterId != null)
                breakendClusterIds.put(breakend, clusterId);
        }

        return breakendClusterIds;
    }

    private static Map<LinxFusion,Integer> findFusionClusterId(
            final List<LinxFusion> fusions, final List<LinxBreakend> breakends, final List<LinxSvAnnotation> svAnnotations)
    {
        Map<LinxFusion,Integer> fusionClusterIds = Maps.newHashMap();

        for(LinxFusion fusion : fusions)
        {
            LinxBreakend breakend = breakends.stream()
                    .filter(x -> x.id() == fusion.fivePrimeBreakendId() || x.id() == fusion.threePrimeBreakendId())
                    .findFirst().orElse(null);

            if(breakend != null)
            {
                Integer clusterId = findBreakendClusterId(breakend, svAnnotations);

                if(clusterId != null)
                    fusionClusterIds.put(fusion, clusterId);

            }
        }

        return fusionClusterIds;
    }

    private static Integer findBreakendClusterId(final LinxBreakend breakend, final List<LinxSvAnnotation> svAnnotations)
    {
        LinxSvAnnotation svAnnotation = svAnnotations.stream().filter(x -> x.svId() == breakend.svId()).findFirst().orElse(null);
        return svAnnotation != null ? svAnnotation.clusterId() : null;
    }

    private static List<LinxGermlineDisruption> selectGermlineSvs(
            final List<LinxGermlineDisruption> germlineSvs, final List<LinxBreakend> reportableGermlineBreakends)
    {
        List<LinxGermlineDisruption> reportableGermlineSvs = new ArrayList<>();

        if(reportableGermlineBreakends == null)
            return reportableGermlineSvs;

        for(LinxGermlineDisruption germlineSv : germlineSvs)
        {
            if(reportableGermlineBreakends.stream().anyMatch(x -> x.svId() == germlineSv.SvId))
            {
                reportableGermlineSvs.add(germlineSv);
            }
        }
        return reportableGermlineSvs;
    }

    private static void restrictToMatchingBreakends(final List<LinxSvAnnotation> svAnnotations, final List<LinxBreakend> breakends)
    {
        int index = 0;

        while(index < svAnnotations.size())
        {
            LinxSvAnnotation svAnnotation = svAnnotations.get(index);

            if(breakends.stream().anyMatch(x -> x.svId() == svAnnotation.svId()))
            {
                ++index;
            }
            else
            {
                svAnnotations.remove(index);
            }
        }
    }

    private static String formPlotTitle(
            final String plotFilename, final List<LinxSvAnnotation> svAnnotations, final List<LinxBreakend> breakends,
            final List<LinxCluster> clusters)
    {
        // use cluster ID to look-up SV count, resolved type and gene list
        // expected file name: SAMPLE_ID.cluster-102.007.png

        try
        {
            int clusterIdIndex = plotFilename.indexOf(SV_VIS_CLUSTER_PREFIX);
            String fileRemainder = plotFilename.substring(clusterIdIndex + SV_VIS_CLUSTER_PREFIX.length());
            int endIndex = fileRemainder.indexOf(".");
            String clusterIdStr = fileRemainder.substring(1, endIndex);
            int clusterId = Integer.parseInt(clusterIdStr);

            LinxCluster cluster = clusters.stream().filter(x -> x.clusterId() == clusterId).findFirst().orElse(null);

            if(cluster == null)
            {
                LOGGER.warn("Linx plot filename({}) cluster({}) not found", plotFilename, clusterId);
                return "";
            }

            Set<String> distinctGenes = Sets.newHashSet();
            int svCount = 0;

            for(LinxSvAnnotation annotation : svAnnotations)
            {
                if(annotation.clusterId() == clusterId)
                {
                    ++svCount;

                    String[] genes = annotation.geneStart().split(ITEM_DELIM, -1);
                    Arrays.stream(genes).forEach(x -> distinctGenes.add(x));

                    genes = annotation.geneEnd().split(ITEM_DELIM, -1);
                    Arrays.stream(genes).forEach(x -> distinctGenes.add(x));
                }
            }

            List<String> geneList = distinctGenes.stream().filter(x -> !x.isEmpty()).collect(Collectors.toList());
            Collections.sort(geneList);
            String geneStr = geneList.stream().collect(Collectors.joining(", "));

            return format("Cluster(%d: %s) SVs(%d) Genes(%s)", clusterId, cluster.resolvedType(), svCount, geneStr);

        }
        catch(Exception e)
        {
            LOGGER.warn("Linx plot filename({}) failed to extract cluster ID", plotFilename);
            return "";
        }
    }
}

