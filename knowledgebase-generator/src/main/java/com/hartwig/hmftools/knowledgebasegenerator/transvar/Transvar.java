package com.hartwig.hmftools.knowledgebasegenerator.transvar;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.knowledgebasegenerator.RefGenomeVersion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class Transvar {

    private static final Logger LOGGER = LogManager.getLogger(Transvar.class);

    private static final int TRANSVAR_TIMEOUT_SEC = 10;

    @NotNull
    private final RefGenomeVersion refGenomeVersion;
    @NotNull
    private final String refGenomeFastaFile;
    @NotNull
    private Map<String, HmfTranscriptRegion> transcriptPerGeneMap;

    public Transvar(@NotNull RefGenomeVersion refGenomeVersion, @NotNull String refGenomeFastaFile) {
        this.refGenomeVersion = refGenomeVersion;
        this.refGenomeFastaFile = refGenomeFastaFile;
        this.transcriptPerGeneMap = HmfGenePanelSupplier.allGenesMap37();
    }

    @NotNull
    public List<VariantHotspot> extractHotspotsFromProteinAnnotation(@NotNull String gene, @NotNull String proteinAnnotation)
            throws IOException, InterruptedException {
        HmfTranscriptRegion transcript = transcriptPerGeneMap.get(gene);
        if (transcript == null) {
            LOGGER.warn("Could not find gene '{}' in HMF gene panel. Skipping hotspot extraction for 'p.{}'", gene, proteinAnnotation);
            return Lists.newArrayList();
        }

        List<TransvarRecord> records = runTransvarPanno(gene, proteinAnnotation);
        if (records.isEmpty()) {
            LOGGER.warn("Could not resolve genomic coordinates. '{}:p.{}' likely does not exist", gene, proteinAnnotation);
            return Lists.newArrayList();
        }

        TransvarRecord best = pickBestRecord(records, transcript.transcriptID());
        if (best == null) {
            LOGGER.warn("Could not find acceptable record amongst {} records for '{}:p.{}'", records.size(), gene, proteinAnnotation);
            return Lists.newArrayList();
        }

        LOGGER.debug("Converting transvar record to hotspots: '{}'", best);
        // This is assuming every transcript on a gene lies on the same strand.
        List<VariantHotspot> hotspots = TransvarInterpreter.convertRecordToHotspots(best, transcript.strand());

        if (hotspots.isEmpty()) {
            LOGGER.warn("Could not derive any hotspots from record {} for  '{}:p.{}'", best, gene, proteinAnnotation);
        }

        return hotspots;
    }

    @Nullable
    private TransvarRecord pickBestRecord(@NotNull List<TransvarRecord> records, @NotNull String preferredTranscriptID) {
        TransvarRecord best = null;

        // We pick a random transcript if the preferred transcript is not found.
        // Ideally we follow "representative transcript from CiViC here
        // See also https://civic.readthedocs.io/en/latest/model/variants/coordinates.html#understanding-coordinates
        for (TransvarRecord record : records) {
            if (record.transcript().equals(preferredTranscriptID)) {
                return record;
            } else {
                best = record;
            }
        }

        return best;
    }

    @NotNull
    private List<TransvarRecord> runTransvarPanno(@NotNull String gene, @NotNull String proteinAnnotation)
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

        // Below is (somehow) necessary to run in intellij. Otherwise it can not find proper locale.
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
    private static List<String> captureStdout(@NotNull Process process) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
        List<String> stdout = toStringList(reader);

        reader.close();

        return stdout;
    }

    @NotNull
    private static List<String> captureStderr(@NotNull Process process) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getErrorStream()));

        List<String> stderr = toStringList(reader);

        reader.close();

        return stderr;
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
