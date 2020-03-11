package com.hartwig.hmftools.knowledgebasegenerator.transvar;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.knowledgebasegenerator.RefGenomeVersion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class TransvarProcess {

    private static final Logger LOGGER = LogManager.getLogger(TransvarProcess.class);

    private static final int TRANSVAR_TIMEOUT_SEC = 60;

    @NotNull
    private final RefGenomeVersion refGenomeVersion;
    @NotNull
    private final String refGenomeFastaFile;

    TransvarProcess(@NotNull RefGenomeVersion refGenomeVersion, @NotNull String refGenomeFastaFile) {
        this.refGenomeVersion = refGenomeVersion;
        this.refGenomeFastaFile = refGenomeFastaFile;
    }

    @NotNull
    List<TransvarRecord> runTransvarPanno(@NotNull String gene, @NotNull String proteinAnnotation)
            throws InterruptedException, IOException {
        ProcessBuilder processBuilder = new ProcessBuilder("transvar",
                "panno",
                "--reference",
                refGenomeFastaFile,
                "--refversion",
                refGenomeVersion.refVersionString(),
                "--noheader",
                "--ensembl",
                "-i",
                gene + ":p." + proteinAnnotation);

        // Below is required on environments where LC_CTYPE is not properly configured (usually on apple).
        processBuilder.environment().put("LC_CTYPE", "UTF-8");

        LOGGER.debug("Running '{}'", command(processBuilder));
        Process process = processBuilder.start();
        if (!process.waitFor(TRANSVAR_TIMEOUT_SEC, TimeUnit.SECONDS)) {
            throw new RuntimeException(String.format("Timeout. '%s' took more than '%s %s' to execute",
                    command(processBuilder),
                    TRANSVAR_TIMEOUT_SEC,
                    TimeUnit.SECONDS));
        }

        if (process.exitValue() != 0) {
            throw new RuntimeException(String.format("'%s' failed with non-zero exit code '%s'",
                    command(processBuilder),
                    process.exitValue()));
        }

        List<String> stderr = captureStderr(process);
        if (!stderr.isEmpty()) {
            LOGGER.warn("Non-empty stderr when running '{}'!", command(processBuilder));
            for (String errLine : stderr) {
                if (!errLine.trim().isEmpty()) {
                    LOGGER.warn(" {}", errLine);
                }
            }
        }

        List<TransvarRecord> records = Lists.newArrayList();
        for (String stdoutLine : captureStdout(process)) {
            LOGGER.debug("Converting transvar output line to TransvarRecord: '{}'", stdoutLine);
            TransvarRecord record = TransvarConverter.toTransvarRecord(stdoutLine);
            if (record != null) {
                records.add(record);
            }
        }

        return records;
    }

    @NotNull
    private static List<String> captureStdout(@NotNull Process process) throws IOException {
        return captureStream(process.getInputStream());
    }

    @NotNull
    private static List<String> captureStderr(@NotNull Process process) throws IOException {
        return captureStream(process.getErrorStream());
    }

    @NotNull
    private static List<String> captureStream(@NotNull InputStream stream) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(stream));

        List<String> strings = toStringList(reader);

        reader.close();

        return strings;

    }

    @NotNull
    private static List<String> toStringList(@NotNull BufferedReader reader) throws IOException {
        List<String> output = Lists.newArrayList();
        String line;
        while ((line = reader.readLine()) != null) {
            output.add(line);
        }
        return output;
    }

    @NotNull
    private static String command(@NotNull ProcessBuilder processBuilder) {
        return processBuilder.command().stream().collect(Collectors.joining(" "));
    }
}
