package com.hartwig.hmftools.qsee.prep.category;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.metrics.GeneDepth;
import com.hartwig.hmftools.common.metrics.GeneDepthFile;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.common.MultiFieldStringBuilder;
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

    private static List<String> getReportableGenes(List<DriverGene> driverGenes, SampleType sampleType)
    {
        List<String> reportableGenes = new ArrayList<>();

        for(DriverGene driverGene : driverGenes)
        {
            boolean isRelevantReportableGene = sampleType == SampleType.TUMOR
                    ? driverGene.reportSomatic()
                    : driverGene.reportGermline();

            if(isRelevantReportableGene)
            {
                reportableGenes.add(driverGene.gene());
            }
        }

        return reportableGenes;
    }

    private static List<Feature> getMissedVariantLikelihoods(List<GeneDepth> geneDepths, List<String> reportableGenes)
    {
        List<GeneDepth> selectedGeneDepths = geneDepths.stream().filter(x -> reportableGenes.contains(x.Gene)).toList();

        List<Feature> features = new ArrayList<>();
        for(GeneDepth geneDepth : selectedGeneDepths)
        {
            String featureName = MultiFieldStringBuilder.formSingleField(FIELD_GENE, geneDepth.Gene);
            FeatureKey key = new FeatureKey(featureName, FeatureType.MISSED_VARIANT_LIKELIHOOD, SOURCE_TOOL);
            Feature feature = new Feature(key, geneDepth.MissedVariantLikelihood);
            features.add(feature);
        }

        return features;
    }

    @Override
    public List<Feature> extractSampleData(String sampleId, @NotNull SampleType sampleType) throws IOException
    {
        List<GeneDepth> geneCoverage = loadGeneCoverage(sampleId, sampleType);
        List<String> reportableGenes = getReportableGenes(mConfig.DriverGenes, sampleType);
        List<Feature> features = getMissedVariantLikelihoods(geneCoverage, reportableGenes);
        return features;
    }
}
