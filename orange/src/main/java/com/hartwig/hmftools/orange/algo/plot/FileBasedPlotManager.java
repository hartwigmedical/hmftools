package com.hartwig.hmftools.orange.algo.plot;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import com.google.common.annotations.VisibleForTesting;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FileBasedPlotManager implements PlotManager
{
    private static final String PLOT_DIRECTORY = "plot";

    @NotNull
    private final String outputDir;

    public FileBasedPlotManager(@NotNull final String outputDir)
    {
        this.outputDir = outputDir;
    }

    @Override
    public void createPlotDirectory() throws IOException
    {
        File plotDir = new File(plotDirectoryPath());
        if(plotDir.exists())
        {
            File[] files = plotDir.listFiles();
            if(files == null)
            {
                throw new IllegalStateException(String.format(
                        "Plot directory of [%s] is not a directory. Please check configured plot directory.",
                        plotDirectoryPath()));
            }
            else if(files.length > 0)
            {
                LOGGER.warn("Plot directory already existed at path [{}], continuing, but output may be mixed with older files. "
                        + "It is recommended to start ORANGE with a clean output directory", plotDirectoryPath());
            }
        }
        else
        {
            if(!plotDir.mkdirs())
            {
                throw new IOException("Unable to create plot directory: " + plotDir);
            }

            LOGGER.debug("Created plot directory '{}'", plotDir.getPath());
        }
    }

    @Nullable
    @Override
    public String processPlotFile(@Nullable String sourcePlotPath) throws IOException
    {
        if(sourcePlotPath == null)
        {
            return null;
        }

        if(!Files.exists(Paths.get(sourcePlotPath)))
        {
            LOGGER.warn("Missing source plot path: " + sourcePlotPath);
            return null;
        }

        String targetPath = checkAddDirSeparator(plotDirectoryPath()) + extractFileName(sourcePlotPath);

        if(!Files.exists(Paths.get(targetPath)))
        {
            LOGGER.debug("Copying '{}' to '{}'", sourcePlotPath, targetPath);
            Files.copy(new File(sourcePlotPath).toPath(), new File(targetPath).toPath());
        }

        return relativePath(targetPath, outputDir);
    }

    @NotNull
    private String plotDirectoryPath()
    {
        return outputDir + File.separator + PLOT_DIRECTORY;
    }

    @NotNull
    @VisibleForTesting
    static String extractFileName(@NotNull String sourcePlotPath)
    {
        return sourcePlotPath.substring(sourcePlotPath.lastIndexOf(File.separator) + 1);
    }

    @NotNull
    @VisibleForTesting
    static String relativePath(@NotNull String target, @NotNull String rootDir)
    {
        String pathToRemove = rootDir.endsWith(File.separator) ? rootDir : rootDir + File.separator;

        if(!target.contains(pathToRemove))
        {
            throw new IllegalStateException("Cannot make relative path of '" + target + "' based on '" + rootDir + "'");
        }

        return target.substring(target.indexOf(pathToRemove) + pathToRemove.length());
    }
}
