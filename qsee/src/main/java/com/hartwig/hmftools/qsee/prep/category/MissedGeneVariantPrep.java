package com.hartwig.hmftools.qsee.prep.category;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.metrics.GeneDepth;
import com.hartwig.hmftools.common.metrics.GeneDepthFile;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;

import org.jetbrains.annotations.NotNull;

public class MissedGeneVariantPrep implements CategoryPrep
{
    private final CommonPrepConfig mConfig;

    private static final SourceTool SOURCE_TOOL = SourceTool.BAM_METRICS;

    private static final String FIELD_GENE = "Gene";

    public MissedGeneVariantPrep(CommonPrepConfig config)
    {
        mConfig = config;
    }

    public SourceTool sourceTool()
    {
        return SOURCE_TOOL;
    }

    private List<GeneDepth> loadGeneCoverage(String sampleId, SampleType sampleType) throws IOException
    {
        String baseDir = mConfig.getBamMetricsDir(sampleId, sampleType);
        String filePath = GeneDepthFile.generateGeneCoverageFilename(baseDir, sampleId);
        return GeneDepthFile.read(filePath);
    }

    private static List<Feature> getMissedVariantLikelihoods(List<GeneDepth> geneDepths, List<DriverGene> driverGenes)
    {
        List<String> selectedGenes = driverGenes.stream().filter(DriverGene::reportSomatic).map(DriverGene::gene).toList();
        List<GeneDepth> selectedGeneDepths = geneDepths.stream().filter(x -> selectedGenes.contains(x.Gene)).toList();

        return selectedGeneDepths.stream()
                .map(x ->
                {
                    String featureName = FeatureKey.formSingleFieldName(FIELD_GENE, x.Gene);
                    FeatureKey key = new FeatureKey(featureName, FeatureType.MISSED_VARIANT_LIKELIHOOD, SOURCE_TOOL);
                    return new Feature(key, x.MissedVariantLikelihood);
                })
                .toList();
    }

    @Override
    public List<Feature> extractSampleData(String sampleId, @NotNull SampleType sampleType) throws IOException
    {
        if(sampleType != SampleType.TUMOR)
        {
            return List.of();
        }

        List<GeneDepth> geneCoverage = loadGeneCoverage(sampleId, sampleType);
        List<Feature> features = getMissedVariantLikelihoods(geneCoverage, mConfig.DriverGenes);
        return features;
    }
}
