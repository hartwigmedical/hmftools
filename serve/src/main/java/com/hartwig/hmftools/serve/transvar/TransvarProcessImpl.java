package com.hartwig.hmftools.serve.transvar;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarRecord;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class TransvarProcessImpl {

    private static final Logger LOGGER = LogManager.getLogger(TransvarProcessImpl.class);

    private static final int TRANSVAR_TIMEOUT_SEC = 90;

    @NotNull
    private final RefGenomeVersion refGenomeVersion;
    @NotNull
    private final String refGenomeFastaFile;

    TransvarProcessImpl(@NotNull RefGenomeVersion refGenomeVersion, @NotNull String refGenomeFastaFile) {
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

        String command = command(processBuilder);

        LOGGER.debug("Running '{}'", command);
        Process process = processBuilder.start();
        if (!process.waitFor(TRANSVAR_TIMEOUT_SEC, TimeUnit.SECONDS)) {
            String warning =
                    String.format("Timeout. '%s' took more than '%s %s' to execute", command, TRANSVAR_TIMEOUT_SEC, TimeUnit.SECONDS);
            LOGGER.warn(warning);
            // We still continue to wait for ever.
            process.waitFor();
        }

        if (process.exitValue() != 0) {
            throw new RuntimeException(String.format("'%s' failed with non-zero exit code '%s'", command, process.exitValue()));
        }

        List<String> stderr = captureStderr(process);
        if (!stderr.isEmpty()) {
            LOGGER.warn("Non-empty stderr when running '{}'!", command);
            for (String errLine : stderr) {
                // First line seems to be generally empty, no need to print that line.
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
