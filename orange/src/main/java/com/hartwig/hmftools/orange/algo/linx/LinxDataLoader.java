package com.hartwig.hmftools.orange.algo.linx;

import static com.hartwig.hmftools.common.driver.DriverType.DRIVERS_LINX_GERMLINE;
import static com.hartwig.hmftools.common.driver.DriverType.DRIVERS_LINX_SOMATIC;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getIntValue;
import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCatalogFile;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.linx.LinxBreakend;
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
                somaticSvAnnotationFile, somaticFusionFile, somaticBreakendFile, somaticDriverCatalogFile, somaticDriversFile,
                config.LinxPlotDirectory, germlineSvAnnotationFile, germlineBreakendFile, germlineDisruptionFile, germlineDriverCatalogFile);
    }

    private static LinxData load(
            final String somaticSvAnnotationFile, final String fusionTsv, final String somaticBreakendTsv,
            final String somaticDriverCatalogTsv, final String somaticDriverTsv, final String linxPlotDirectory,
            @Nullable String germlineSvAnnotationFile, @Nullable String germlineBreakendTsv, @Nullable String germlineDisruptionTsv,
            @Nullable String germlineDriverCatalogTsv)
            throws IOException
    {
        List<LinxDriver> somaticDriverData = LinxDriver.read(somaticDriverTsv);
        List<LinxFusion> fusions = loadFusions(fusionTsv);
        List<LinxBreakend> somaticBreakends = loadBreakends(somaticBreakendTsv, fusions);

        List<DriverCatalog> somaticDrivers = DriverCatalogFile.read(somaticDriverCatalogTsv).stream()
                .filter(x -> DRIVERS_LINX_SOMATIC.contains(x.driver())).collect(Collectors.toList());

        // limit SV annotations to reportable breakends
        List<LinxSvAnnotation> somaticSvAnnotations = LinxSvAnnotation.read(somaticSvAnnotationFile);

        Map<LinxBreakend,Integer> somaticBreakendClusterIds = findBreakendClusterId(somaticBreakends, somaticSvAnnotations);

        Map<LinxFusion,Integer> fusionClusterIds = findFusionClusterId(fusions, somaticBreakends, somaticSvAnnotations);

        // likely now redundant since the list will not be persisted
        restrictToMatchingBreakends(somaticSvAnnotations, somaticBreakends);

        List<DriverCatalog> somaticHomozygousDisruptions = somaticDrivers.stream()
                .filter(x -> x.driver() == DriverType.HOM_DUP_DISRUPTION || x.driver() == DriverType.HOM_DEL_DISRUPTION)
                .collect(Collectors.toList());

        List<LinxSvAnnotation> germlineSvAnnotations = null;
        List<LinxBreakend> reportableGermlineBreakends = null;
        List<DriverCatalog> germlineDrivers = null;
        List<LinxGermlineDisruption> reportableGermlineDisruptions = null;

        if(germlineSvAnnotationFile != null && germlineBreakendTsv != null && germlineDisruptionTsv != null)
        {
            germlineDrivers = DriverCatalogFile.read(germlineDriverCatalogTsv).stream()
                    .filter(x -> DRIVERS_LINX_GERMLINE.contains(x.driver())).collect(Collectors.toList());

            germlineSvAnnotations = LinxSvAnnotation.read(germlineSvAnnotationFile);
            reportableGermlineBreakends = loadBreakends(germlineBreakendTsv, Collections.emptyList());

            restrictToMatchingBreakends(germlineSvAnnotations, reportableGermlineBreakends);

            List<LinxGermlineDisruption> allGermlineDisruptions = LinxGermlineDisruption.read(germlineDisruptionTsv);
            reportableGermlineDisruptions = selectGermlineSvs(allGermlineDisruptions, reportableGermlineBreakends);
        }

        List<String> reportableEventPlots = Lists.newArrayList();

        if(linxPlotDirectory != null)
        {
            for(String file : new File(linxPlotDirectory).list())
            {
                reportableEventPlots.add(linxPlotDirectory + File.separator + file);
            }

            LOGGER.debug(" loaded {} Linx plots from {}", reportableEventPlots.size(), linxPlotDirectory);
        }

        return ImmutableLinxData.builder()
                .somaticSvAnnotations(somaticSvAnnotations)
                .somaticDrivers(somaticDrivers)
                .somaticDriverData(somaticDriverData)
                .fusions(fusions)
                .fusionClusterIds(fusionClusterIds)
                .somaticBreakends(somaticBreakends)
                .somaticBreakendClusterIds(somaticBreakendClusterIds)
                .somaticHomozygousDisruptions(somaticHomozygousDisruptions)
                .reportableEventPlots(reportableEventPlots)
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

    private static Set<Integer> loadFusionClusters(final String filename) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(filename).toPath());
        if(lines.isEmpty())
        {
            throw new IllegalStateException(String.format("File lacks header: %s", filename));
        }

        Set<Integer> clusterIds = new HashSet<>();
        String header = lines.get(0);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        lines.remove(0);

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM);
            clusterIds.add(getIntValue(fieldsIndexMap, "ClusterId", 0, values));
        }

        return clusterIds;
    }

    private static void loadSvToCluster(
            final String somaticVisSvDataTsv, final Map<Integer,Integer> svIdToClusterId, final Map<Integer,Integer> clusterIdToLinkCount) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(somaticVisSvDataTsv).toPath());
        if(lines.isEmpty())
        {
            throw new IllegalStateException(String.format("File lacks header: %s", somaticVisSvDataTsv));
        }

        String header = lines.get(0);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        lines.remove(0);

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM);
            int clusterId = getIntValue(fieldsIndexMap, "ClusterId", 0, values);
            int svId = getIntValue(fieldsIndexMap, "SvId", 0, values);

            if(!svIdToClusterId.containsKey(svId))
            {
                svIdToClusterId.put(svId, clusterId);
            }

            if(clusterIdToLinkCount.containsKey(clusterId))
            {
                clusterIdToLinkCount.put(clusterId, clusterIdToLinkCount.get(clusterId) + 1);
            }
            else
            {
                clusterIdToLinkCount.put(clusterId, 1);
            }
        }
    }

    private static Map<Integer, Integer> loadClusterExonCounts(final String filename) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(filename).toPath());
        if(lines.isEmpty())
        {
            throw new IllegalStateException(String.format("File lacks header: %s", filename));
        }

        Map<Integer, Integer> counts = new HashMap<>();
        String header = lines.get(0);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        lines.remove(0);

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM);
            int clusterId = getIntValue(fieldsIndexMap, "ClusterId", 0, values);

            if(counts.containsKey(clusterId))
            {
                counts.put(clusterId, counts.get(clusterId) + 1);
            }
            else
            {
                counts.put(clusterId, 1);
            }
        }

        return counts;
    }
}

