package com.hartwig.hmftools.serve.transvar;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.RefGenomeVersion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class Transvar {

    private static final Logger LOGGER = LogManager.getLogger(Transvar.class);

    @NotNull
    private final TransvarProcess process;
    @NotNull
    private final TransvarInterpreter interpreter;
    @NotNull
    private final Map<String, HmfTranscriptRegion> transcriptPerGeneMap;

    @NotNull
    public static Transvar withRefGenome(@NotNull RefGenomeVersion refGenomeVersion, @NotNull String refGenomeFastaFile)
            throws FileNotFoundException {
        return new Transvar(new TransvarProcess(refGenomeVersion, refGenomeFastaFile),
                TransvarInterpreter.fromRefGenomeFastaFile(refGenomeFastaFile),
                HmfGenePanelSupplier.allGenesMap37());
    }

    private Transvar(@NotNull TransvarProcess process, @NotNull TransvarInterpreter interpreter,
            @NotNull Map<String, HmfTranscriptRegion> transcriptPerGeneMap) {
        this.process = process;
        this.interpreter = interpreter;
        this.transcriptPerGeneMap = transcriptPerGeneMap;
    }

    @NotNull
    public List<VariantHotspot> extractHotspotsFromProteinAnnotation(@NotNull String gene, @Nullable String transcriptID,
            @NotNull String proteinAnnotation) throws IOException, InterruptedException {
        List<TransvarRecord> records = process.runTransvarPanno(gene, proteinAnnotation);
        if (records.isEmpty()) {
            LOGGER.warn("Transvar could not resolve any genomic coordinates for '{}:p.{}'", gene, proteinAnnotation);
            return Lists.newArrayList();
        }

        HmfTranscriptRegion canonicalTranscript = transcriptPerGeneMap.get(gene);
        if (canonicalTranscript == null) {
            LOGGER.warn("Could not find canonical transcript for '{}' in HMF gene panel. Skipping hotspot extraction for 'p.{}'",
                    gene,
                    proteinAnnotation);
            return Lists.newArrayList();
        }

        TransvarRecord best = pickBestRecord(records, transcriptID, canonicalTranscript.transcriptID());

        LOGGER.debug("Interpreting transvar record: '{}'", best);
        // This is assuming every transcript on a gene lies on the same strand.
        List<VariantHotspot> hotspots = interpreter.convertRecordToHotspots(best, canonicalTranscript.strand());

        if (hotspots.isEmpty()) {
            LOGGER.warn("Could not derive any hotspots from record {} for '{}:p.{}'", best, gene, proteinAnnotation);
        }

        return hotspots;
    }

    @NotNull
    private TransvarRecord pickBestRecord(@NotNull List<TransvarRecord> records, @Nullable String preferredTranscriptId,
            @NotNull String canonicalTranscriptId) {
        assert !records.isEmpty();

        TransvarRecord preferredRecord = null;
        TransvarRecord canonicalRecord = null;
        TransvarRecord bestRecord = null;
        for (TransvarRecord record : records) {
            if (preferredTranscriptId != null && record.transcript().equals(preferredTranscriptId)) {
                preferredRecord = record;
            } else if (record.transcript().equals(canonicalTranscriptId)) {
                canonicalRecord = record;
            } else {
                bestRecord = record;
            }
        }

        if (preferredRecord != null) {
            return preferredRecord;
        }

        if (canonicalRecord != null) {
            return canonicalRecord;
        }

        assert bestRecord != null;
        return bestRecord;
    }
}
