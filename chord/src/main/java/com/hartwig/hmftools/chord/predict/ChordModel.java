package com.hartwig.hmftools.chord.predict;

import static com.hartwig.hmftools.chord.ChordConstants.APP_NAME;
import static com.hartwig.hmftools.chord.ChordConstants.CHORD_LOGGER;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.RExecutor;

import org.jetbrains.annotations.NotNull;

public class ChordModel
{
    private final String mModelPath;

    private static final String SCRIPT_RESOURCE_PATH = "chord_predict.R";
    private static final String MODEL_RESOURCE_PATH = "CHORD.rds";
    private static final String MODEL_TMP_PATH = System.getProperty("java.io.tmpdir") + "/CHORD.rds";

    private ChordModel(String modelPath)
    {
        mModelPath = modelPath;
    }

    public static ChordModel fromResources()
    {
        return new ChordModel(MODEL_TMP_PATH);
    }

    public static void copyModelToTmpPath() throws IOException
    {
        InputStream source = ChordModel.class.getClassLoader().getResourceAsStream(MODEL_RESOURCE_PATH);
        Path target = new File(MODEL_TMP_PATH).toPath();

        Files.copy(source, target, StandardCopyOption.REPLACE_EXISTING);
    }

    public void predict(String mutContextsFile, String outputFile)
    {
        CHORD_LOGGER.info("Starting CHORD predict");

        try
        {
            copyModelToTmpPath();

            List<String> scriptArgs = new ArrayList<>(List.of(
                    mModelPath,
                    mutContextsFile,
                    outputFile,
                    CHORD_LOGGER.getLevel().toString()
            ));

            int result = RExecutor.executeFromClasspath(
                    SCRIPT_RESOURCE_PATH,
                    true,
                    scriptArgs.toArray(new String[0])
            );

            if(result != 0)
            {
                throw new IOException("R execution failed for script: " + SCRIPT_RESOURCE_PATH);
            }
        }
        catch(Exception e)
        {
            CHORD_LOGGER.error("Failed to get predictions from CHORD model: " + e);
            e.printStackTrace();
            System.exit(1);
        }
        finally
        {
            File modelTmpPath = new File(MODEL_TMP_PATH);
            if(modelTmpPath.isFile())
                modelTmpPath.delete();
        }

        CHORD_LOGGER.info("Completed CHORD predict");
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        PredictConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        ChordModel model = ChordModel.fromResources();

        PredictConfig config = new PredictConfig(configBuilder);
        model.predict(config.MutContextsFile, config.OutputFile);
    }
}
