package com.hartwig.hmftools.serve.transvar;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarRecord;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class TransvarProcessImpl implements TransvarProcess {

    private static final Logger LOGGER = LogManager.getLogger(TransvarProcessImpl.class);

    // DUPs give rise to errors in transvar but do lead to interpretable output so error can be ignored.
    private static final Set<String> ERRORS_TO_IGNORE = Sets.newHashSet("warning: unknown alternative: DUP, ignore alternative");

    private static final int TRANSVAR_TIMEOUT_SEC = 90;

    @NotNull
    private final RefGenomeVersion refGenomeVersion;
    @NotNull
    private final String refGenomeFastaFile;
    @NotNull
    private final TransvarCurator curator;

    TransvarProcessImpl(@NotNull RefGenomeVersion refGenomeVersion, @NotNull String refGenomeFastaFile) {
        this.refGenomeVersion = refGenomeVersion;
        this.refGenomeFastaFile = refGenomeFastaFile;
        this.curator = new TransvarCurator(refGenomeVersion);
    }

    @Override
    @NotNull
    public List<TransvarRecord> runTransvarPanno(@NotNull String gene, @NotNull String proteinAnnotation)
            throws InterruptedException, IOException {
        ProcessBuilder processBuilder = buildProcessBuilder(gene, proteinAnnotation);
        String command = command(processBuilder);

        LOGGER.debug("Running '{}'", command);
        Process process = processBuilder.start();
        if (!process.waitFor(TRANSVAR_TIMEOUT_SEC, TimeUnit.SECONDS)) {
            LOGGER.info("Timeout. '{}' took more than '{} {}' to execute", command, TRANSVAR_TIMEOUT_SEC, TimeUnit.SECONDS);
            // We still continue to wait for ever. Assume transvar will eventually succeed.
            process.waitFor();
        }

        if (process.exitValue() != 0) {
            throw new RuntimeException(String.format("'%s' failed with non-zero exit code '%s'", command, process.exitValue()));
        }

        List<String> stderr = filterErrorOutput(captureStderr(process));
        if (!stderr.isEmpty()) {
            LOGGER.warn("Non-empty stderr when running '{}'!", command);
            for (String errLine : stderr) {
                LOGGER.warn(" {}", errLine);
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
    private ProcessBuilder buildProcessBuilder(@NotNull String gene, @NotNull String proteinAnnotation) {
        String curatedGene = curator.curateGene(gene);
        String curatedProteinAnnotation = curator.curateProteinAnnotation(proteinAnnotation);

        ProcessBuilder processBuilder = new ProcessBuilder("transvar",
                "panno",
                "--reference",
                refGenomeFastaFile,
                "--refversion",
                toTransvarRefVersion(refGenomeVersion),
                "--noheader",
                "--ensembl",
                "-i",
                curatedGene + ":p." + curatedProteinAnnotation);

        // Below is required on environments where LC_CTYPE is not properly configured (usually on apple).
        processBuilder.environment().put("LC_CTYPE", "UTF-8");

        return processBuilder;
    }

    @NotNull
    private static String toTransvarRefVersion(@NotNull RefGenomeVersion refGenomeVersion) {
        switch (refGenomeVersion) {
            case V37:
                return "hg19";
            case V38:
                return "hg38";
            default:
                throw new IllegalStateException("Could not convert ref genome version to transvar ref version: " + refGenomeVersion);
        }
    }

    @NotNull
    private static List<String> filterErrorOutput(@NotNull List<String> errors) {
        List<String> filteredErrors = Lists.newArrayList();
        for (String error : errors) {
            if (!error.trim().isEmpty()) {
                boolean foundIgnoreString = false;
                for (String ignoreString : ERRORS_TO_IGNORE) {
                    if (error.contains(ignoreString)) {
                        LOGGER.debug("Suppressing error from transvar: '{}'", error);
                        foundIgnoreString = true;
                        break;
                    }
                }

                if (!foundIgnoreString) {
                    filteredErrors.add(error);
                }
            }
        }

        return filteredErrors;
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
