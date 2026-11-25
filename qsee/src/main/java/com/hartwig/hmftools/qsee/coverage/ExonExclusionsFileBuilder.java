package com.hartwig.hmftools.qsee.coverage;

import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.qsee.common.SampleType;

public class ExonExclusionsFileBuilder
{
    private final ExonExclusionsConfig mConfig;

    public ExonExclusionsFileBuilder(ExonExclusionsConfig config)
    {
        mConfig = config;
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
        if(mConfig.Qsee.getSampleIds(SampleType.TUMOR).size() <= 1)
        {
            QC_LOGGER.error("Building exon exclusions requires a cohort of tumor samples");
            System.exit(1);
        }

        List<ExonCoverageQc> coverageQcPerExon = new ExonCoverageQcFileBuilder(mConfig.Qsee).run();
        ExonExclusions exonExclusions = calcExonExclusions(coverageQcPerExon, mConfig.MinLowCoverageSampleFraction);

        String filename = ExonExclusions.generateFilename(mConfig.Qsee.OutputDir, mConfig.Qsee.OutputId);
        exonExclusions.write(filename);
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        ExonExclusionsConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        ExonExclusionsConfig config = new ExonExclusionsConfig(configBuilder);
        ExonExclusionsFileBuilder exonExclusionsFileBuilder = new ExonExclusionsFileBuilder(config);
        exonExclusionsFileBuilder.run();
    }
}
