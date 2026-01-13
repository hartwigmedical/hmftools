package com.hartwig.hmftools.orange.algo.linx;

import static com.hartwig.hmftools.common.linx.LinxCommonTypes.generateVisExonFilename;
import static com.hartwig.hmftools.common.linx.LinxCommonTypes.generateVisFusionFilename;
import static com.hartwig.hmftools.common.linx.LinxCommonTypes.generateVisSvFilename;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getIntValue;

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
        String linxSomaticDir = config.linxSomaticDataDirectory();
        String linxGermlineDir = config.wgsRefConfig() != null ? config.wgsRefConfig().linxGermlineDataDirectory() : null;
        String tumorSample = config.tumorSampleId();

        String somaticSvAnnotationFile = LinxSvAnnotation.generateFilename(linxSomaticDir, tumorSample, false);
        String somaticBreakendFile = LinxBreakend.generateFilename(linxSomaticDir, tumorSample);
        String somaticFusionFile = LinxFusion.generateFilename(linxSomaticDir, tumorSample);
        String somaticDriversFile = LinxDriver.generateFilename(linxSomaticDir, tumorSample);
        String somaticDriverCatalogFile = LinxDriver.generateCatalogFilename(linxSomaticDir, tumorSample, true);
        String somaticVisFusionFile = generateVisFusionFilename(linxSomaticDir, tumorSample, false);
        String somaticVisSvDataFile = generateVisSvFilename(linxSomaticDir, tumorSample, false);
        String somaticVisGeneExonFile = generateVisExonFilename(linxSomaticDir, tumorSample, false);

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
                somaticVisFusionFile, somaticVisSvDataFile, somaticVisGeneExonFile,
                germlineSvAnnotationFile, germlineBreakendFile, germlineDisruptionFile, germlineDriverCatalogFile);
    }

    private static LinxData load(
            final String somaticSvAnnotationFile, final String fusionTsv, final String somaticBreakendTsv,
            final String somaticDriverCatalogTsv, final String somaticDriverTsv,
            final String somaticVisFusionTsv, final String somaticVisSvDataTsv, final String somaticVisGeneExonTsv,
            @Nullable String germlineSvAnnotationFile, @Nullable String germlineBreakendTsv, @Nullable String germlineDisruptionTsv,
            @Nullable String germlineDriverCatalogTsv)
            throws IOException
    {
        List<LinxDriver> somaticDrivers = LinxDriver.read(somaticDriverTsv);
        List<LinxFusion> reportableSomaticFusions = loadReportableFusions(fusionTsv);
        List<LinxBreakend> reportableSomaticBreakends = loadReportableBreakends(somaticBreakendTsv);

        // limit SV annotations to reportable breakends
        List<LinxSvAnnotation> somaticSvAnnotations = LinxSvAnnotation.read(somaticSvAnnotationFile);
        restrictToMatchingBreakends(somaticSvAnnotations, reportableSomaticBreakends);

        List<HomozygousDisruption> somaticHomozygousDisruptions = extractHomozygousDisruptions(somaticDriverCatalogTsv);

        Set<Integer> fusionClusterIds = loadFusionClusters(somaticVisFusionTsv);

        Map<Integer,Integer> svIdToClusterId = Maps.newHashMap();
        Map<Integer,Integer> clusterIdToLinkCount = Maps.newHashMap();

        loadSvToCluster(somaticVisSvDataTsv, svIdToClusterId, clusterIdToLinkCount);
        Map<Integer, Integer> clusterIdToExonCount = loadClusterExonCounts(somaticVisGeneExonTsv);

        List<LinxSvAnnotation> germlineSvAnnotations = null;
        List<LinxBreakend> reportableGermlineBreakends = null;
        List<LinxGermlineDisruption> reportableGermlineDisruptions = null;
        List<HomozygousDisruption> germlineHomozygousDisruptions = Collections.emptyList();

        if(germlineSvAnnotationFile != null && germlineBreakendTsv != null && germlineDisruptionTsv != null)
        {
            germlineSvAnnotations = LinxSvAnnotation.read(germlineSvAnnotationFile);
            reportableGermlineBreakends = loadReportableBreakends(germlineBreakendTsv);

            restrictToMatchingBreakends(germlineSvAnnotations, reportableGermlineBreakends);

            List<LinxGermlineDisruption> allGermlineDisruptions = LinxGermlineDisruption.read(germlineDisruptionTsv);
            reportableGermlineDisruptions = selectReportableGermlineSvs(allGermlineDisruptions, reportableGermlineBreakends);
        }

        return ImmutableLinxData.builder()
                .somaticSvAnnotations(somaticSvAnnotations)
                .somaticDrivers(somaticDrivers)
                .fusions(reportableSomaticFusions)
                .somaticBreakends(reportableSomaticBreakends)
                .somaticHomozygousDisruptions(somaticHomozygousDisruptions)
                .fusionClusterIds(fusionClusterIds)
                .svIdToClusterId(svIdToClusterId)
                .clusterIdToLinkCount(clusterIdToLinkCount)
                .clusterIdToExonCount(clusterIdToExonCount)
                .germlineSvAnnotations(germlineSvAnnotations)
                .germlineBreakends(reportableGermlineBreakends)
                .germlineDisruptions(reportableGermlineDisruptions)
                .germlineHomozygousDisruptions(germlineHomozygousDisruptions)
                .build();
    }

    private static List<LinxFusion> loadReportableFusions(final String fusionTsv) throws IOException
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

    private static List<LinxBreakend> loadReportableBreakends(final String somaticBreakendTsv) throws IOException
    {
        List<LinxBreakend> breakends = LinxBreakend.read(somaticBreakendTsv);

        List<LinxBreakend> reportableBreakends = new ArrayList<>();
        for(LinxBreakend breakend : breakends)
        {
            if(breakend.reportedStatus() == ReportedStatus.REPORTED)
            {
                reportableBreakends.add(breakend);
            }
        }
        return reportableBreakends;
    }

    private static List<LinxGermlineDisruption> selectReportableGermlineSvs(
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

    public static List<HomozygousDisruption> extractHomozygousDisruptions(final String driverCatalogTsv)
            throws IOException
    {
        List<DriverCatalog> linxDriversCatalog = DriverCatalogFile.read(driverCatalogTsv);

        List<HomozygousDisruption> homozygousDisruptions = extractSomaticHomozygousDisruptions(linxDriversCatalog);
        return homozygousDisruptions;
    }

    private static List<HomozygousDisruption> extractSomaticHomozygousDisruptions(final List<DriverCatalog> driverCatalog)
    {
        List<HomozygousDisruption> homozygousDisruptions = Lists.newArrayList();

        for(DriverCatalog driver : driverCatalog)
        {
            if(driver.driver() == DriverType.HOM_DUP_DISRUPTION || driver.driver() == DriverType.HOM_DEL_DISRUPTION)
            {
                homozygousDisruptions.add(create(driver));
            }
        }

        return homozygousDisruptions;
    }

    private static HomozygousDisruption create(final DriverCatalog driverCatalog)
    {
        return ImmutableHomozygousDisruption.builder()
                .chromosome(driverCatalog.chromosome())
                .chromosomeBand(driverCatalog.chromosomeBand())
                .gene(driverCatalog.gene())
                .transcript(driverCatalog.transcript())
                .isCanonical(driverCatalog.isCanonical())
                .build();
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

