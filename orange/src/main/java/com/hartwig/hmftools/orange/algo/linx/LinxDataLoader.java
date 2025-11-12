package com.hartwig.hmftools.orange.algo.linx;

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

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.linx.LinxGermlineDisruption;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.orange.OrangeConfig;

import org.jetbrains.annotations.Nullable;

public final class LinxDataLoader
{
    private static final String VIS_FUSION_FILE_EXTENSION = ".linx.vis_fusion.tsv";
    private static final String VIS_SV_DATA_FILE_EXTENSION = ".linx.vis_sv_data.tsv";
    private static final String VIS_GENE_EXON_FILE_EXTENSION = ".linx.vis_gene_exon.tsv";

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
        String somaticVisFusionFile = generateVisFusionFilename(linxSomaticDir, tumorSample);
        String somaticVisSvDataFile = generateVisSvDataFilename(linxSomaticDir, tumorSample);
        String somaticVisGeneExonFile = generateVisGeneExonFilename(linxSomaticDir, tumorSample);

        String germlineSvAnnotationFile = null;
        String germlineBreakendFile = null;
        String germlineDisruptionFile = null;
        String germlineDriverCatalogFile = null;

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
                germlineSvAnnotationFile, germlineBreakendFile, germlineDisruptionFile, germlineDriverCatalogFile, config.includeNonGenePanelEvents());
    }

    private static LinxData load(final String somaticStructuralVariantTsv, final String somaticFusionTsv,
            final String somaticBreakendTsv,
            final String somaticDriverCatalogTsv, final String somaticDriverTsv,
            final String somaticVisFusionTsv, final String somaticVisSvDataTsv, final String somaticVisGeneExonTsv,
            @Nullable String germlineStructuralVariantTsv,
            @Nullable String germlineBreakendTsv, @Nullable String germlineDisruptionTsv, @Nullable String germlineDriverCatalogTsv,
            boolean includeNonGenePanelEvents)
            throws IOException
    {
        List<LinxSvAnnotation> allSomaticStructuralVariants = includeNonGenePanelEvents ?
                LinxSvAnnotation.read(somaticStructuralVariantTsv) : Collections.emptyList();

        List<LinxDriver> allSomaticDrivers = LinxDriver.read(somaticDriverTsv);

        List<LinxFusion> allSomaticFusions = LinxFusion.read(somaticFusionTsv);
        List<LinxFusion> reportableSomaticFusions = selectReportableFusions(allSomaticFusions);

        if(!includeNonGenePanelEvents)
            allSomaticFusions.clear();

        List<LinxBreakend> allSomaticBreakends = LinxBreakend.read(somaticBreakendTsv);
        List<LinxBreakend> reportableSomaticBreakends = selectReportableBreakends(allSomaticBreakends);

        if(!includeNonGenePanelEvents)
            allSomaticBreakends.clear();

        List<HomozygousDisruption> somaticHomozygousDisruptions = HomozygousDisruptionFactory.extractSomaticFromLinxDriverCatalogTsv(
                somaticDriverCatalogTsv);

        Set<Integer> fusionClusterIds = loadFusionClusters(somaticVisFusionTsv);

        Map<Integer,Integer> svIdToClusterId = Maps.newHashMap();
        Map<Integer,Integer> clusterIdToLinkCount = Maps.newHashMap();

        loadSvToCluster(somaticVisSvDataTsv, svIdToClusterId, clusterIdToLinkCount);
        Map<Integer, Integer> clusterIdToExonCount = loadClusterExonCounts(somaticVisGeneExonTsv);

        List<LinxSvAnnotation> allGermlineStructuralVariants = null;
        if(germlineStructuralVariantTsv != null)
        {
            allGermlineStructuralVariants = LinxSvAnnotation.read(germlineStructuralVariantTsv);
        }

        List<LinxBreakend> allGermlineBreakends = null;
        List<LinxBreakend> reportableGermlineBreakends = null;
        if(germlineBreakendTsv != null)
        {
            allGermlineBreakends = LinxBreakend.read(germlineBreakendTsv);
            reportableGermlineBreakends = selectReportableBreakends(allGermlineBreakends);
        }

        List<LinxGermlineDisruption> allGermlineDisruptions = null;
        List<LinxGermlineDisruption> reportableGermlineDisruptions = null;
        if(germlineDisruptionTsv != null)
        {
            allGermlineDisruptions = LinxGermlineDisruption.read(germlineDisruptionTsv);
            reportableGermlineDisruptions = selectReportableGermlineSvs(allGermlineDisruptions, reportableGermlineBreakends);
        }

        List<HomozygousDisruption> germlineHomozygousDisruptions = null;
        if(germlineDriverCatalogTsv != null)
        {
            germlineHomozygousDisruptions = Collections.emptyList();
        }

        return ImmutableLinxData.builder()
                .allSomaticSvAnnotations(allSomaticStructuralVariants)
                .somaticDrivers(allSomaticDrivers)
                .allSomaticFusions(allSomaticFusions)
                .reportableSomaticFusions(reportableSomaticFusions)
                .allSomaticBreakends(allSomaticBreakends)
                .reportableSomaticBreakends(reportableSomaticBreakends)
                .somaticHomozygousDisruptions(somaticHomozygousDisruptions)
                .fusionClusterIds(fusionClusterIds)
                .svIdToClusterId(svIdToClusterId)
                .clusterIdToLinkCount(clusterIdToLinkCount)
                .clusterIdToExonCount(clusterIdToExonCount)
                .allGermlineSvAnnotations(allGermlineStructuralVariants)
                .allGermlineBreakends(allGermlineBreakends)
                .reportableGermlineBreakends(reportableGermlineBreakends)
                .allGermlineDisruptions(allGermlineDisruptions)
                .reportableGermlineDisruptions(reportableGermlineDisruptions)
                .germlineHomozygousDisruptions(germlineHomozygousDisruptions)
                .build();
    }

    private static List<LinxFusion> selectReportableFusions(final List<LinxFusion> fusions)
    {
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

    private static List<LinxBreakend> selectReportableBreakends(final List<LinxBreakend> breakends)
    {
        List<LinxBreakend> reportableBreakends = new ArrayList<>();
        for(LinxBreakend breakend : breakends)
        {
            if(breakend.reportedDisruption())
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
        {
            return reportableGermlineSvs;
        }

        for(LinxGermlineDisruption germlineSv : germlineSvs)
        {

            if(reportableGermlineBreakends.stream().anyMatch(x -> x.svId() == germlineSv.SvId))
            {
                reportableGermlineSvs.add(germlineSv);
            }
        }
        return reportableGermlineSvs;
    }

    private static String generateVisFusionFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + VIS_FUSION_FILE_EXTENSION;
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

    private static String generateVisSvDataFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + VIS_SV_DATA_FILE_EXTENSION;
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

    private static String generateVisGeneExonFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + VIS_GENE_EXON_FILE_EXTENSION;
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

