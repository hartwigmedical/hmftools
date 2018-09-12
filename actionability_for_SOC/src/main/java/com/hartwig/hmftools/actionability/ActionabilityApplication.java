package com.hartwig.hmftools.actionability;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.sql.SQLException;
import java.util.List;
import java.io.File;

import com.hartwig.hmftools.actionability.variants.ActionabilityVariantsAnalyzer;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.io.path.PathExtensionFinder;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionabilityApplication {
    private static final org.apache.logging.log4j.Logger LOGGER = LogManager.getLogger(ActionabilityApplication.class);
    //add tumor locations of database
    private static final String RUN_DIRECTORY = "run_dir";
    private static final String SOMATIC_VCF_EXTENSION_V3 = "_post_processed_v2.2.vcf.gz";
    private static final String SOMATIC_VCF_EXTENSION_V4 = "_post_processed.vcf.gz";

    public static void main(final String... args) throws ParseException, IOException, SQLException {
        LOGGER.info("Determining actionability variants.");
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);
        final String runDir = cmd.getOptionValue(RUN_DIRECTORY);
        final RunContext run = ProductionRunContextFactory.fromRunDirectory(runDir);
        final List<SomaticVariant> variants = loadPassedSomaticVariants(run.tumorSample(), runDir);
        LOGGER.info("Tumor sample: " + run.tumorSample());

        LOGGER.info("");
        LOGGER.info("Start processing actionability variants");

        String fileActionabilityVariants = "/data/common/dbs/knowledgebases/output/actionableVariants.tsv";
     //   String fileActionabilityRanges = "/data/common/dbs/knowledgebases/output/actionableRanges.tsv";

        if (Files.exists(new File(fileActionabilityVariants).toPath())) {
            ActionabilityVariantsAnalyzer analyzer = ActionabilityVariantsAnalyzer.loadFromFile(fileActionabilityVariants);
        } else {
            LOGGER.warn("File does not exist: " + fileActionabilityVariants);
        }

        // look for evy variant in file from fileActionabilityVariants
      //  LOGGER.info("variants van tumorSample" + variants.get(1).alt());

        LOGGER.info("");
        LOGGER.info("Start processing actionability fusions");
     //   String fileActionabilityFusionPairs = "/data/common/dbs/knowledgebases/output/actionableFusionPairs.tsv";
     //   String fileActionabilityPromiscuousFive = "/data/common/dbs/knowledgebases/output/actionablePromiscuousFive.tsv";
      //  String fileActionabilityPromiscuousThree = "/data/common/dbs/knowledgebases/output/actionablePromiscuousThree.tsv";


        LOGGER.info("");
        LOGGER.info("Start processing actionability cnvs");
       // String fileActionabilityCNVs = "/data/common/dbs/knowledgebases/output/actionableCNVs.tsv";


        LOGGER.info("");
        LOGGER.info("Writing output data to file");

        LOGGER.info("");
        LOGGER.info("Finish orocessing actionability variants");
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(RUN_DIRECTORY, true, "Complete path towards a single run dir where patient reporter will run on.");
        //add tumor locations of database
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    public static List<SomaticVariant> loadPassedSomaticVariants(@NotNull final String sample, @NotNull final String path) throws IOException {
        // TODO (KODU): Clean up once pipeline v3 no longer exists
        Path vcfPath;
        try {
            vcfPath = PathExtensionFinder.build().findPath(path, SOMATIC_VCF_EXTENSION_V3);
        } catch (FileNotFoundException exception) {
            vcfPath = PathExtensionFinder.build().findPath(path, SOMATIC_VCF_EXTENSION_V4);
        }
        return SomaticVariantFactory.passOnlyInstance().fromVCFFile(sample, vcfPath.toString());
    }
}
