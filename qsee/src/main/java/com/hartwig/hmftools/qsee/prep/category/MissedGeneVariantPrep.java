package com.hartwig.hmftools.qsee.prep.category;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.metrics.GeneDepth;
import com.hartwig.hmftools.common.metrics.GeneDepthFile;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;

import org.jetbrains.annotations.NotNull;

public class MissedGeneVariantPrep implements CategoryPrep
{
    private final CommonPrepConfig mConfig;

    public MissedGeneVariantPrep(CommonPrepConfig config)
    {
        mConfig = config;
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
                .map(x -> new Feature(FeatureType.MISSED_VARIANT_LIKELIHOOD, x.Gene, x.MissedVariantLikelihood, SourceTool.BAM_METRICS))
                .toList();
    }

    @Override
    public List<Feature> extractSampleData(String sampleId, @NotNull SampleType sampleType) throws IOException
    {
        List<GeneDepth> geneCoverage = loadGeneCoverage(sampleId, sampleType);
        List<Feature> features = getMissedVariantLikelihoods(geneCoverage, mConfig.DriverGenes);
        return features;
    }
}
