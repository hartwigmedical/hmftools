package com.hartwig.hmftools.qsee.coverage;

import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.common.metrics.ExonDepth;
import com.hartwig.hmftools.common.metrics.GeneDepthFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.qsee.QseeConfig;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.category.MissedGeneVariantPrep;

import org.apache.logging.log4j.Level;

public class ExonCoverageQcFileBuilder
{
    private final QseeConfig mConfig;

    public ExonCoverageQcFileBuilder(QseeConfig config)
    {
        mConfig = config;
    }

    private List<ExonDepth> loadExonCoverage(String sampleId, SampleType sampleType) throws IOException
    {
        String baseDir = mConfig.getBamMetricsDir(sampleId, sampleType);
        String filePath = GeneDepthFile.generateExonCoverageFilename(baseDir, sampleId);
        return GeneDepthFile.readExonCoverage(filePath);
    }

    public List<ExonCoverageQc> run()
    {
        List<ExonCoverageQc> coverageQcPerExon = new ArrayList<>();

        for(SampleType sampleType : List.of(SampleType.TUMOR, SampleType.NORMAL))
        {
            for(String sampleId : mConfig.getSampleIds(sampleType))
            {
                try
                {
                    List<ExonDepth> exonDepths = loadExonCoverage(sampleId, sampleType);
                    List<String> reportableGenes = MissedGeneVariantPrep.getReportableGenes(mConfig.DriverGenes, sampleType);

                    List<ExonCoverageQc> coverageQcPerExonForSample = ExonCoverageQc.fromExonDepths(
                            sampleId,
                            sampleType,
                            exonDepths,
                            mConfig.ExcludedExons,
                            reportableGenes,
                            mConfig.QcThresholds
                    );

                    coverageQcPerExon.addAll(coverageQcPerExonForSample);
                }
                catch(IOException e)
                {
                    Level logLevel = mConfig.AllowMissingInput ? Level.WARN : Level.ERROR;
                    QC_LOGGER.log(logLevel, "sampleType({}) exon variant LoD - missing input file(s): {}", sampleType, e.getMessage());

                    if(!mConfig.AllowMissingInput)
                        System.exit(1);
                }
            }
        }

        if(!mConfig.isSinglePatient())
            coverageQcPerExon.sort(Comparator.comparing(ExonCoverageQc::sampleId));

        String exonCoverageQcFile = ExonCoverageQcFile.generateFilename(mConfig);
        ExonCoverageQcFile.write(exonCoverageQcFile, coverageQcPerExon, mConfig.WriteLodDetails);

        return coverageQcPerExon;
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        QseeConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        QseeConfig config = new QseeConfig(configBuilder);
        ExonCoverageQcFileBuilder exonCoverageQcFileBuilder = new ExonCoverageQcFileBuilder(config);
        exonCoverageQcFileBuilder.run();
    }
}
