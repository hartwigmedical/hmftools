package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.common.CommonUtils.buildComparers;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.FieldConfigFile;
import com.hartwig.hmftools.compar.common.MatchLevel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// generates the default field config file (thresholds, comparison type etc. per field) for every category, for reference in the README
public class DefaultFieldConfigTestApplication
{
    private static final Logger LOGGER = LogManager.getLogger(DefaultFieldConfigTestApplication.class);

    private static final String OUTPUT_DIR = "target/field_config";

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

        buildComparers(config);

        String outputDir = OUTPUT_DIR + "/" + matchLevel.toString().toLowerCase();
        new File(outputDir).mkdirs();

        String filename = FieldConfigFile.generateFileName(outputDir);

        FieldConfigFile.write(filename, config.FieldConfig, config.Categories.keySet());

        LOGGER.info("wrote default field config file: {}", filename);
    }
}