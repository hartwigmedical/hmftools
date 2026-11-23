package com.hartwig.hmftools.qsee.coverage;

import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;
import static com.hartwig.hmftools.qsee.coverage.LimitOfDetectionCalcs.LOD_VAF_WARN_THRESHOLD;

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
    private static final String MIN_HIGH_LOD_VAF_SAMPLE_FRACTION = "min_high_lod_vaf_sample_fraction";
    private static final String MIN_HIGH_LOD_VAF_SAMPLE_FRACTION_DESC = "Min fraction of tumor samples with LodVaf > " + LOD_VAF_WARN_THRESHOLD + " for an exon to be excluded";
    private static final double DEFAULT_MIN_HIGH_LOD_VAF_SAMPLE_FRACTION = 0.95;

    private final QseeConfig mConfig;
    private final double mMinHighLodVafSampleFraction;

    public ExonExclusionsFileBuilder(QseeConfig config, double minHighLodVafSampleFraction)
    {
        mConfig = config;
        mMinHighLodVafSampleFraction = minHighLodVafSampleFraction;
    }

    private record Exon(String geneName, int exonRank) {}

    @VisibleForTesting
    static ExonExclusions calcExonExclusions(List<ExonCoverageQc> coverageQcPerExon, double minHighLodVafSampleFraction)
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

            long highLodVafSampleCount = exonCoverageQcPerSample.stream().filter(x -> x.lodVaf() > LOD_VAF_WARN_THRESHOLD).count();
            int minHighLodVafSampleCount = (int) Math.ceil(exonCoverageQcPerSample.size() * minHighLodVafSampleFraction);

            if(highLodVafSampleCount >= minHighLodVafSampleCount)
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
        ExonExclusions exonExclusions = calcExonExclusions(coverageQcPerExon, mMinHighLodVafSampleFraction);

        String filename = ExonExclusions.generateFilename(mConfig.OutputDir, mConfig.OutputId);
        exonExclusions.write(filename);
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        QseeConfig.registerConfig(configBuilder);
        configBuilder.addDecimal(MIN_HIGH_LOD_VAF_SAMPLE_FRACTION, MIN_HIGH_LOD_VAF_SAMPLE_FRACTION_DESC, DEFAULT_MIN_HIGH_LOD_VAF_SAMPLE_FRACTION);
        configBuilder.checkAndParseCommandLine(args);

        QseeConfig config = new QseeConfig(configBuilder);
        double minHighLodVafSampleFraction = configBuilder.getDecimal(MIN_HIGH_LOD_VAF_SAMPLE_FRACTION);
        ExonExclusionsFileBuilder exonExclusionsFileBuilder = new ExonExclusionsFileBuilder(config, minHighLodVafSampleFraction);
        exonExclusionsFileBuilder.run();
    }
}
