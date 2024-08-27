package com.hartwig.hmftools.virusinterpreter;

import static com.hartwig.hmftools.virusinterpreter.VirusInterpreterConfig.registerConfig;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.AnnotatedVirusFile;
import com.hartwig.hmftools.common.virus.VirusBreakend;
import com.hartwig.hmftools.common.virus.VirusBreakendFile;
import com.hartwig.hmftools.virusinterpreter.algo.VirusBlacklistingDbFile;
import com.hartwig.hmftools.virusinterpreter.algo.VirusReportingDbFile;
import com.hartwig.hmftools.virusinterpreter.algo.VirusReportingDbModel;
import com.hartwig.hmftools.virusinterpreter.coverages.CoveragesAnalysis;
import com.hartwig.hmftools.virusinterpreter.coverages.CoveragesAnalyzer;
import com.hartwig.hmftools.virusinterpreter.taxonomy.TaxonomyDb;
import com.hartwig.hmftools.virusinterpreter.taxonomy.TaxonomyDbFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class VirusInterpreterApplication
{
    public static final Logger VI_LOGGER = LogManager.getLogger(VirusInterpreterApplication.class);

    public static void main(String[] args) throws IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder("VirusInterpreter");
        registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        VirusInterpreterConfig config = new VirusInterpreterConfig(configBuilder);

        VI_LOGGER.info("Loading taxonomy db from {}", config.TaxonomyDbTsv);
        TaxonomyDb taxonomyDb = TaxonomyDbFile.loadFromTsv(config.TaxonomyDbTsv);

        VI_LOGGER.info("Loading virus blacklisting db from {}", config.VirusBlacklistedDbTsv);
        List<Integer> blacklistedTaxids = VirusBlacklistingDbFile.loadFromTsv(config.VirusBlacklistedDbTsv);

        VI_LOGGER.info("Building virus reporting db model from {}", config.VirusReportedDbTsv);
        VirusReportingDbModel virusReportingDbModel = VirusReportingDbFile.buildFromTsv(config.VirusReportedDbTsv);

        VI_LOGGER.info("Loading virus breakends from {}", new File(config.VirusBreakendTsv).getParent());
        List<VirusBreakend> virusBreakends = VirusBreakendFile.read(config.VirusBreakendTsv);
        VI_LOGGER.info("Loaded {} virus breakends from {}", virusBreakends.size(), config.VirusBreakendTsv);

        VI_LOGGER.info("Loading purity context from Purple dir {}", config.PurpleDir);
        PurityContext purityContext = PurityContextFile.read(config.PurpleDir, config.SampleId);

        VI_LOGGER.info("Running coverage analysis based on tumor WGS metrics file {}", config.TumorSampleWGSMetricsFile);
        CoveragesAnalysis coveragesAnalysis =
                CoveragesAnalyzer.run(purityContext, config.TumorSampleWGSMetricsFile);
        VI_LOGGER.info("Determined the expected clonal coverage to be {}", coveragesAnalysis.ExpectedClonalCoverage);

        VirusInterpreterAlgo algo = new VirusInterpreterAlgo(taxonomyDb, blacklistedTaxids, virusReportingDbModel, coveragesAnalysis);
        List<AnnotatedVirus> annotatedViruses = algo.analyze(virusBreakends, purityContext);
        VI_LOGGER.info("Interpreter classified {} viruses as reportable", annotatedViruses.stream().filter(x -> x.reported()).count());

        String annotatedVirusTsv = AnnotatedVirusFile.generateFileName(config.OutputDir, config.SampleId);
        AnnotatedVirusFile.write(annotatedVirusTsv, annotatedViruses);
        VI_LOGGER.info("Written {} annotated viruses to {}", annotatedViruses.size(), annotatedVirusTsv);
    }

}
