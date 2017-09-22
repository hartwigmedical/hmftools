package com.hartwig.hmftools.metadata;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.io.FolderChecker;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class MetadataCreator {

    private static final Logger LOGGER = LogManager.getLogger(MetadataCreator.class);

    private static final String RUNS_DIR = "runs_dir";

    public static void main(String[] args) throws ParseException, IOException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        final String runsDir = cmd.getOptionValue(RUNS_DIR);
        if (runsDir == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Metadata Creator", options);
            System.exit(1);
        }
        new MetadataCreator(runsDir).run();
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(RUNS_DIR, true, "Path towards runs");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private final String runsDir;

    MetadataCreator(@NotNull final String runsDir) {
        this.runsDir = runsDir;
    }

    @VisibleForTesting
    void run() throws IOException {
        final File runDirectories = new File(runsDir);
        if (runDirectories.isDirectory()) {
            final File[] folders = runDirectories.listFiles(File::isDirectory);
            if (folders == null) {
                throw new IOException("List files in " + runDirectories.getName() + " returned null.");
            }
            for (final File folder : folders) {
                String runDirectory = null;
                try {
                    runDirectory = FolderChecker.build().checkFolder(folder.getPath());
                } catch (HartwigException e) {
                    LOGGER.error("Could not generate run directory!");
                }

                final String metaDataPath = runDirectory + File.separator + "metadata";

                if (runDirectory != null) {
                    try {
                        ProductionRunContextFactory.fromRunDirectory(runDirectory);
                        LOGGER.info("Skipping " + runDirectory + ": can already generate run context.");
                    } catch (HartwigException e) {
                        LOGGER.info("Generating meta data for " + runDirectory);
                        final List<String> metaData = Lists.newArrayList(makeMetaDataStringForRun(runDirectory));
                        LOGGER.info("Writing metadata to " + metaDataPath);
                        Files.write(new File(metaDataPath).toPath(), metaData);
                        try {
                            ProductionRunContextFactory.fromRunDirectory(runDirectory);
                        } catch (HartwigException e1) {
                            LOGGER.error(" -> Still cannot resolve run context for " + runDirectory + "!");
                        }
                    }
                }
            }
        }
    }

    @NotNull
    private static String makeMetaDataStringForRun(@NotNull final String runDirPath) throws IOException {
        final File runDir = new File(runDirPath);
        String metadata = "{\"set_name\": %set%, \"ref_sample\": \"%ref%\", \"tumor_sample\": \"%tum%\"}";
        metadata = metadata.replace("%set%", runDir.getName());
        final File[] folders = runDir.listFiles(File::isDirectory);
        assert folders != null;
        for (File folder : folders) {
            String name = folder.getName();
            if (name.contains("CPCT") || name.contains("DRUP")) {
                name = name.substring(5, name.length());
                if (name.contains("R")) {
                    metadata = metadata.replace("%ref%", folder.getName());
                } else if (name.contains("T")) {
                    metadata = metadata.replace("%tum%", folder.getName());
                }
            }
        }

        LOGGER.info(metadata);
        return metadata;
    }
}
