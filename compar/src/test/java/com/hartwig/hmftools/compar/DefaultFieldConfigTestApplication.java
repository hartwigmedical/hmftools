package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.common.CommonUtils.buildComparers;
import static com.hartwig.hmftools.compar.common.CommonUtils.initialiseFieldConfig;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.FieldConfig;
import com.hartwig.hmftools.compar.common.FieldConfigFile;
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

    public static void main(final String[] args) throws IOException
    {
        for(MatchLevel matchLevel : MatchLevel.values())
        {
            writeFieldConfigFile(matchLevel);
        }
    }

    private static void writeFieldConfigFile(final MatchLevel matchLevel) throws IOException
    {
        ComparConfig config = new ComparConfig();

        for(CategoryType category : CategoryType.values())
        {
            config.Categories.put(category, matchLevel);
        }

        FieldConfig fieldConfig = initialiseFieldConfig(config);

        String outputDir = OUTPUT_DIR + "/" + matchLevel.toString().toLowerCase();
        new File(outputDir).mkdirs();

        String filename = FieldConfigFile.generateFileName(outputDir);

        Set<CategoryType> categories = buildComparers(config).stream()
                .map(c -> c.category())
                .collect(Collectors.toSet());
        FieldConfigFile.write(filename, fieldConfig, categories);

        LOGGER.info("wrote default field config file: {}", filename);
    }
}