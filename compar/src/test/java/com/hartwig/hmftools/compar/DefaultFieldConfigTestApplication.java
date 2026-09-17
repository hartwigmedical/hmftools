package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.compar.ComparConfig.CATEGORIES;
import static com.hartwig.hmftools.compar.ComparConfig.MATCH_LEVEL;
import static com.hartwig.hmftools.compar.common.CategoryType.ALL_CATEGORIES;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.List;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.compar.common.MatchLevel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// generates the default field config file (thresholds, comparison type etc. per field) for every category, for reference in the README
public class DefaultFieldConfigTestApplication
{
    private static final Logger LOGGER = LogManager.getLogger(DefaultFieldConfigTestApplication.class);

    // resolved from this class's own compiled location (compar/target/test-classes) rather than the JVM's working
    // directory, which varies depending on how the application is launched (IDE, Maven, module vs project root)
    private static final String OUTPUT_DIR = new File(moduleTargetDir(), "default_field_config").getPath();

    private static File moduleTargetDir()
    {
        try
        {
            File codeSourceLocation =
                    new File(DefaultFieldConfigTestApplication.class.getProtectionDomain().getCodeSource().getLocation().toURI());

            // codeSourceLocation is compar/target/test-classes (or the test jar); its parent is compar/target
            return codeSourceLocation.getParentFile();
        }
        catch(URISyntaxException e)
        {
            throw new RuntimeException(e);
        }
    }

    public static void main(final String[] args)
    {
        for(MatchLevel matchLevel : MatchLevel.values())
        {
            writeFieldConfigFile(matchLevel);
        }
    }

    private static void writeFieldConfigFile(final MatchLevel matchLevel)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        ComparConfig.addConfig(configBuilder);

        configBuilder.setValue(CATEGORIES, ALL_CATEGORIES);
        configBuilder.setValue(MATCH_LEVEL, matchLevel.toString());

        configBuilder.setValue(SAMPLE, "FAKE");
        ComparConfig config = new ComparConfig(configBuilder);

        String outputDir = OUTPUT_DIR + "/" + matchLevel.toString().toLowerCase();
        new File(outputDir).mkdirs();

        String filePrefix = outputDir + "/compar_cohort";
        String filename = FieldConfigFile.generateFileName(filePrefix);

        FieldCheckCache fieldConfig = new FieldCheckCache(); // initialiseFieldConfig(config);
        List<ItemComparer> comparers = ComparerUtils.buildComparers(config, fieldConfig);
        FieldConfigFile.write(filename, comparers);

        LOGGER.info("wrote default field config file: {}", filename);
    }
}
