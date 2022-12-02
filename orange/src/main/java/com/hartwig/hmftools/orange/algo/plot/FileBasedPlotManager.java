package com.hartwig.hmftools.orange.algo.plot;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FileBasedPlotManager implements PlotManager {

    private static final Logger LOGGER = LogManager.getLogger(FileBasedPlotManager.class);

    private static final String PLOT_DIRECTORY = "plot";

    @NotNull
    private final String outputDir;

    public FileBasedPlotManager(@NotNull final String outputDir) {
        this.outputDir = outputDir;
    }

    @Override
    public void createPlotDirectory() throws IOException {
        File plotDir = new File(plotDirectoryPath());
        if (plotDir.exists()) {
            throw new IOException("Plot directory exists already: " + plotDir);
        }

        if (!plotDir.mkdirs()) {
            throw new IOException("Unable to create plot directory: " + plotDir);
        }

        LOGGER.debug("Created plot directory '{}'", plotDir.getPath());
    }

    @Nullable
    @Override
    public String processPlotFile(@Nullable String sourcePlotPath) throws IOException {
        if (sourcePlotPath == null) {
            return null;
        }

        String targetPath = plotDirectoryPath() + File.separator + extractFileName(sourcePlotPath);

        LOGGER.debug("Copying '{}' to '{}'", sourcePlotPath, targetPath);
        Files.copy(new File(sourcePlotPath).toPath(), new File(targetPath).toPath());

        return relativePath(targetPath, outputDir);
    }

    @NotNull
    private String plotDirectoryPath() {
        return outputDir + File.separator + PLOT_DIRECTORY;
    }

    @NotNull
    @VisibleForTesting
    static String extractFileName(@NotNull String sourcePlotPath) {
        return sourcePlotPath.substring(sourcePlotPath.lastIndexOf(File.separator) + 1);
    }

    @NotNull
    @VisibleForTesting
    static String relativePath(@NotNull String target, @NotNull String rootDir) {
        String pathToRemove = rootDir.endsWith(File.separator) ? rootDir : rootDir + File.separator;

        if (!target.contains(pathToRemove)) {
            throw new IllegalStateException("Cannot make relative path of '" + target + "' based on '" + rootDir + "'");
        }

        return target.substring(target.indexOf(pathToRemove) + pathToRemove.length());
    }
}
