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
    public List<VariantHotspot> extractHotspotsFromProteinAnnotation(@NotNull String gene, @NotNull String proteinAnnotation)
            throws IOException, InterruptedException {
        HmfTranscriptRegion transcript = transcriptPerGeneMap.get(gene);
        if (transcript == null) {
            LOGGER.warn("Could not find gene '{}' in HMF gene panel. Skipping hotspot extraction for 'p.{}'", gene, proteinAnnotation);
            return Lists.newArrayList();
        }

        List<TransvarRecord> records = process.runTransvarPanno(gene, proteinAnnotation);
        if (records.isEmpty()) {
            LOGGER.warn("Transvar could not resolve genomic coordinates for '{}:p.{}'", gene, proteinAnnotation);
            return Lists.newArrayList();
        }

        TransvarRecord best = pickBestRecord(records, transcript.transcriptID());
        if (best == null) {
            LOGGER.warn("Could not find acceptable record amongst {} records for '{}:p.{}'", records.size(), gene, proteinAnnotation);
            return Lists.newArrayList();
        }

        if (!best.transcript().equals(transcript.transcriptID())) {
            LOGGER.debug("Best transvar record not on canonical transcript for '{}:p.{}'", gene, proteinAnnotation);
        }

        LOGGER.debug("Interpreting transvar record: '{}'", best);
        // This is assuming every transcript on a gene lies on the same strand.
        List<VariantHotspot> hotspots = interpreter.convertRecordToHotspots(best, transcript.strand());

        if (hotspots.isEmpty()) {
            LOGGER.warn("Could not derive any hotspots from record {} for  '{}:p.{}'", best, gene, proteinAnnotation);
        }

        return hotspots;
    }

    @Nullable
    private TransvarRecord pickBestRecord(@NotNull List<TransvarRecord> records, @NotNull String preferredTranscriptID) {
        TransvarRecord best = null;

        // We pick a random transcript if the preferred (canonical) transcript is not found.
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
}
