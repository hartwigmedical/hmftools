package com.hartwig.hmftools.common.linx;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getIntValue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class LinxDataLoader
{
    private static final String VIS_FUSION_FILE_EXTENSION = ".linx.vis_fusion.tsv";
    private static final String VIS_SV_DATA_FILE_EXTENSION = ".linx.vis_sv_data.tsv";

    @NotNull
    public static LinxData load(final String tumorSample, final String linxSomaticDir, @Nullable String linxGermlineDir)
            throws IOException
    {
        final String somaticSvAnnotationFile = LinxSvAnnotation.generateFilename(linxSomaticDir, tumorSample, false);
        final String somaticBreakendFile = LinxBreakend.generateFilename(linxSomaticDir, tumorSample);
        final String somaticFusionFile = LinxFusion.generateFilename(linxSomaticDir, tumorSample);
        final String somaticDriversFile = LinxDriver.generateFilename(linxSomaticDir, tumorSample);
        final String somaticDriverCatalogFile = LinxDriver.generateCatalogFilename(linxSomaticDir, tumorSample, true);
        final String somaticVisFusionFile = generateVisFusionFilename(linxSomaticDir, tumorSample);
        final String somaticVisSvDataFile = generateVisSvDataFilename(linxSomaticDir, tumorSample);

        String germlineSvAnnotationFile = null;
        String germlineBreakendFile = null;
        String germlineDisruptionFile = null;
        String germlineDriverCatalogFile = null;

        if(linxGermlineDir != null)
        {
            germlineSvAnnotationFile = LinxSvAnnotation.generateFilename(linxGermlineDir, tumorSample, true);
            germlineBreakendFile = LinxBreakend.generateFilename(linxGermlineDir, tumorSample, true);
            germlineDisruptionFile = LinxGermlineSv.generateFilename(linxGermlineDir, tumorSample);
            germlineDriverCatalogFile = LinxDriver.generateCatalogFilename(linxGermlineDir, tumorSample, false);
        }

        return load(somaticSvAnnotationFile,
                somaticFusionFile,
                somaticBreakendFile,
                somaticDriverCatalogFile,
                somaticDriversFile,
                somaticVisFusionFile,
                somaticVisSvDataFile,
                germlineSvAnnotationFile,
                germlineBreakendFile,
                germlineDisruptionFile,
                germlineDriverCatalogFile);
    }

    @NotNull
    private static LinxData load(@NotNull String somaticStructuralVariantTsv, @NotNull String somaticFusionTsv,
            @NotNull String somaticBreakendTsv,
            @NotNull String somaticDriverCatalogTsv, @NotNull String somaticDriverTsv,
            @NotNull String somaticVisFusionTsv, @NotNull String somaticVisSvDataTsv,
            @Nullable String germlineStructuralVariantTsv,
            @Nullable String germlineBreakendTsv, @Nullable String germlineDisruptionTsv, @Nullable String germlineDriverCatalogTsv)
            throws IOException
    {
        List<LinxSvAnnotation> allSomaticStructuralVariants = LinxSvAnnotation.read(somaticStructuralVariantTsv);
        List<LinxDriver> allSomaticDrivers = LinxDriver.read(somaticDriverTsv);

        List<LinxFusion> allSomaticFusions = LinxFusion.read(somaticFusionTsv);
        List<LinxFusion> reportableSomaticFusions = selectReportableFusions(allSomaticFusions);

        List<LinxBreakend> allSomaticBreakends = LinxBreakend.read(somaticBreakendTsv);
        List<LinxBreakend> reportableSomaticBreakends = selectReportableBreakends(allSomaticBreakends);

        List<HomozygousDisruption> somaticHomozygousDisruptions =
                HomozygousDisruptionFactory.extractSomaticFromLinxDriverCatalogTsv(somaticDriverCatalogTsv);

        Set<Integer> fusionClusterIds = loadFusionClusters(somaticVisFusionTsv);
        Map<Integer, Integer> svIdToClusterId = loadSvToCluster(somaticVisSvDataTsv);

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

        List<LinxGermlineSv> allGermlineDisruptions = null;
        List<LinxGermlineSv> reportableGermlineDisruptions = null;
        if(germlineDisruptionTsv != null)
        {
            allGermlineDisruptions = LinxGermlineSv.read(germlineDisruptionTsv);
            reportableGermlineDisruptions = selectReportableGermlineSvs(allGermlineDisruptions, reportableGermlineBreakends);
        }

        List<HomozygousDisruption> germlineHomozygousDisruptions = null;
        if(germlineDriverCatalogTsv != null)
        {
            germlineHomozygousDisruptions = HomozygousDisruptionFactory.extractGermlineFromLinxDriverCatalogTsv(germlineDriverCatalogTsv);
        }

        return ImmutableLinxData.builder()
                .allSomaticStructuralVariants(allSomaticStructuralVariants)
                .somaticDrivers(allSomaticDrivers)
                .allSomaticFusions(allSomaticFusions)
                .reportableSomaticFusions(reportableSomaticFusions)
                .allSomaticBreakends(allSomaticBreakends)
                .reportableSomaticBreakends(reportableSomaticBreakends)
                .somaticHomozygousDisruptions(somaticHomozygousDisruptions)
                .fusionClusterIds(fusionClusterIds)
                .putAllSvIdToClusterId(svIdToClusterId)
                .allGermlineStructuralVariants(allGermlineStructuralVariants)
                .allGermlineBreakends(allGermlineBreakends)
                .reportableGermlineBreakends(reportableGermlineBreakends)
                .allGermlineDisruptions(allGermlineDisruptions)
                .reportableGermlineDisruptions(reportableGermlineDisruptions)
                .germlineHomozygousDisruptions(germlineHomozygousDisruptions)
                .build();
    }

    @NotNull
    private static List<LinxFusion> selectReportableFusions(@NotNull List<LinxFusion> fusions)
    {
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        for(LinxFusion fusion : fusions)
        {
            if(fusion.reported())
            {
                reportableFusions.add(fusion);
            }
        }
        return reportableFusions;
    }

    @NotNull
    private static List<LinxBreakend> selectReportableBreakends(@NotNull List<LinxBreakend> breakends)
    {
        List<LinxBreakend> reportableBreakends = Lists.newArrayList();
        for(LinxBreakend breakend : breakends)
        {
            if(breakend.reportedDisruption())
            {
                reportableBreakends.add(breakend);
            }
        }
        return reportableBreakends;
    }

    @NotNull
    private static List<LinxGermlineSv> selectReportableGermlineSvs(
            @NotNull List<LinxGermlineSv> germlineSvs, @NotNull List<LinxBreakend> reportableGermlineBreakends)
    {
        List<LinxGermlineSv> reportableGermlineSvs = Lists.newArrayList();

        if(reportableGermlineBreakends == null)
        {
            return reportableGermlineSvs;
        }

        for(LinxGermlineSv germlineSv : germlineSvs)
        {

            if(reportableGermlineBreakends.stream().anyMatch(x -> x.svId() == germlineSv.SvId))
            {
                reportableGermlineSvs.add(germlineSv);
            }
        }
        return reportableGermlineSvs;
    }

    @NotNull
    private static String generateVisFusionFilename(@NotNull String basePath, @NotNull String sample)
    {
        return basePath + File.separator + sample + VIS_FUSION_FILE_EXTENSION;
    }

    @NotNull
    private static Set<Integer> loadFusionClusters(@NotNull String filename) throws IOException
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

    @NotNull
    private static String generateVisSvDataFilename(@NotNull String basePath, @NotNull String sample)
    {
        return basePath + File.separator + sample + VIS_SV_DATA_FILE_EXTENSION;
    }

    @NotNull
    private static Map<Integer, Integer> loadSvToCluster(@NotNull String filename) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(filename).toPath());
        if(lines.isEmpty())
        {
            throw new IllegalStateException(String.format("File lacks header: %s", filename));
        }

        Map<Integer, Integer> svToCluster = new HashMap<>();
        String header = lines.get(0);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        lines.remove(0);

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM);
            int clusterId = getIntValue(fieldsIndexMap, "ClusterId", 0, values);
            int svId = getIntValue(fieldsIndexMap, "SvId", 0, values);

            if(!svToCluster.containsKey(svId))
            {
                svToCluster.put(svId, clusterId);
            }
        }

        return svToCluster;
    }
}
