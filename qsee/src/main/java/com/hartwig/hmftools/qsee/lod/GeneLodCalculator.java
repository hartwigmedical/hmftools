package com.hartwig.hmftools.qsee.lod;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.metrics.ExonDepth;
import com.hartwig.hmftools.common.metrics.GeneDepthFile;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.qsee.QseeConfig;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.category.MissedGeneVariantPrep;

public class GeneLodCalculator
{
    private final QseeConfig mConfig;

    // Assuming diploid genome
    private static final double TUMOR_VAF_CLONAL_HOM_VARIANT_DETECTABLE = 1;
    private static final double TUMOR_VAF_CLONAL_HET_VARIANT_DETECTABLE = 0.5;

    public GeneLodCalculator(QseeConfig config)
    {
        mConfig = config;
    }

    private List<ExonDepth> loadExonCoverage(String sampleId, SampleType sampleType) throws IOException
    {
        String baseDir = mConfig.getBamMetricsDir(sampleId, sampleType);
        String filePath = GeneDepthFile.generateExonCoverageFilename(baseDir, sampleId);
        return GeneDepthFile.readExonCoverage(filePath);
    }

    private double loadTumorPurity(String sampleId, SampleType sampleType) throws IOException
    {
        if(sampleType != SampleType.TUMOR)
            return LimitOfDetectionCalcs.REF_SAMPLE_PURITY;

        String baseDir = mConfig.getPurpleDir(sampleId);
        String purityFile = PurplePurity.generateFilename(baseDir, sampleId);
        return PurplePurity.read(purityFile).Purity;
    }

    private static VariantDetectabilityStatus determineVariantDetectabilityStatus(double tumorVafLod)
    {
        if(tumorVafLod >= TUMOR_VAF_CLONAL_HOM_VARIANT_DETECTABLE)
            return VariantDetectabilityStatus.FAIL;

        if(tumorVafLod > TUMOR_VAF_CLONAL_HET_VARIANT_DETECTABLE)
            return VariantDetectabilityStatus.WARN_HET_UNRELIABLE;

        return VariantDetectabilityStatus.PASS;
    }

    @VisibleForTesting
    static List<GeneLod> calcGeneLods(
            String sampleId, SampleType sampleType, List<ExonDepth> exonDepths, List<String> reportableGenes, double tumorPurity)
    {
        List<GeneLod> geneLods = new ArrayList<>();

        for(ExonDepth exonDepth : exonDepths)
        {
            if(!reportableGenes.contains(exonDepth.Gene))
                continue;

            double tumorVafLod = LimitOfDetectionCalcs.calcMinTumorVafToDetectVariant(exonDepth.DepthCounts, tumorPurity);
            VariantDetectabilityStatus variantDetectabilityStatus = determineVariantDetectabilityStatus(tumorVafLod);

            GeneLod geneLod = new GeneLod(
                    sampleId, sampleType,
                    exonDepth.Gene, exonDepth.ExonRank, exonDepth.MedianDepth,
                    tumorVafLod, variantDetectabilityStatus
            );

            geneLods.add(geneLod);
        }

        return geneLods;
    }

    public List<GeneLod> calcSampleGeneLods(String sampleId, SampleType sampleType) throws IOException
    {
        List<ExonDepth> exonDepths = loadExonCoverage(sampleId, sampleType);
        List<String> reportableGenes = MissedGeneVariantPrep.getReportableGenes(mConfig.DriverGenes, sampleType);
        double tumorPurity = loadTumorPurity(sampleId, sampleType);

        return calcGeneLods(sampleId, sampleType, exonDepths, reportableGenes, tumorPurity);
    }
}
