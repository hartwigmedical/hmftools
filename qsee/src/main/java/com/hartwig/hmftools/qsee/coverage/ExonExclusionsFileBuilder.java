package com.hartwig.hmftools.qsee.coverage;

import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.qsee.QseeConfig;
import com.hartwig.hmftools.qsee.common.SampleType;

public class ExonExclusionsFileBuilder
{
    private static final String MIN_LOW_COVERAGE_SAMPLE_FRACTION = "min_low_coverage_sample_fraction";
    private static final String MIN_LOW_COVERAGE_SAMPLE_FRACTION_DESC = "Min fraction of tumor samples with exon status WARN for an exon to be excluded";
    private static final double DEFAULT_MIN_LOW_COVERAGE_SAMPLE_FRACTION = 0.95;

    private final QseeConfig mConfig;
    private final double mMinLowCoverageSampleFraction;

    public ExonExclusionsFileBuilder(QseeConfig config, double minLowCoverageSampleFraction)
    {
        mConfig = config;
        mMinLowCoverageSampleFraction = minLowCoverageSampleFraction;
    }

    private record Exon(String geneName, int exonRank) {}

    @VisibleForTesting
    static ExonExclusions calcExonExclusions(List<ExonCoverageQc> coverageQcPerExon, double minLowCoverageSampleFraction)
    {
        Map<Exon,List<ExonCoverageQc>> coverageQcByExon = new HashMap<>();
        for(ExonCoverageQc exonCoverageQc : coverageQcPerExon)
        {
            if(exonCoverageQc.sampleType() != SampleType.TUMOR)
                continue;

            Exon exon = new Exon(exonCoverageQc.geneName(), exonCoverageQc.exonRank());
            coverageQcByExon.computeIfAbsent(exon, k -> new ArrayList<>()).add(exonCoverageQc);
        }

        ExonExclusions exonExclusions = new ExonExclusions();
        for(Exon exon : coverageQcByExon.keySet())
        {
            List<ExonCoverageQc> exonCoverageQcPerSample = coverageQcByExon.get(exon);

            long lowCoverageSampleCount = exonCoverageQcPerSample.stream().filter(x -> x.variantDetectabilityStatus() == VariantDetectabilityStatus.WARN).count();
            int minLowCoverageSampleCount = (int) Math.ceil(exonCoverageQcPerSample.size() * minLowCoverageSampleFraction);

            if(lowCoverageSampleCount >= minLowCoverageSampleCount)
                exonExclusions.add(exon.geneName(), exon.exonRank());
        }

        return exonExclusions;
    }

    public void run()
    {
        if(mConfig.getSampleIds(SampleType.TUMOR).size() <= 1)
        {
            QC_LOGGER.error("Building exon exclusions requires a cohort of tumor samples");
            System.exit(1);
        }

        List<ExonCoverageQc> coverageQcPerExon = new ExonCoverageQcFileBuilder(mConfig).run();
        ExonExclusions exonExclusions = calcExonExclusions(coverageQcPerExon, mMinLowCoverageSampleFraction);

        String filename = ExonExclusions.generateFilename(mConfig.OutputDir, mConfig.OutputId);
        exonExclusions.write(filename);
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        QseeConfig.registerConfig(configBuilder);
        configBuilder.addDecimal(MIN_LOW_COVERAGE_SAMPLE_FRACTION, MIN_LOW_COVERAGE_SAMPLE_FRACTION_DESC, DEFAULT_MIN_LOW_COVERAGE_SAMPLE_FRACTION);
        configBuilder.checkAndParseCommandLine(args);

        QseeConfig config = new QseeConfig(configBuilder);
        double minLowCoverageSampleFraction = configBuilder.getDecimal(MIN_LOW_COVERAGE_SAMPLE_FRACTION);
        ExonExclusionsFileBuilder exonExclusionsFileBuilder = new ExonExclusionsFileBuilder(config, minLowCoverageSampleFraction);
        exonExclusionsFileBuilder.run();
    }
}
