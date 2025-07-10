package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_TRANS_NAME;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_GC_TARGET;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_GC_TOLERANCE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENERAL_GC_TARGET;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENERAL_GC_TOLERANCE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_EXON_FLANK_GAP;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_EXON_FLANK_REGION;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_LONG_INTRON_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_MAX_EXONS_TO_ADD_INTRON;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_MIN_INTRON_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_UPDOWNSTREAM_GAP;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_UPDOWNSTREAM_REGION;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_QUALITY_ACCEPT;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_QUALITY_REJECT;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

// Probes covering (regions of) selected genes.
// Methodology for gene regions:
//   - Coding: Cover the full coding region of each exon.
//   - UTR: 1 probe centered on each noncoding exon.
//   - Upstream/downstream: Select the best acceptable probe from a ~2kb region ~1kb upstream/downstream.
//   - Intronic: Only when there are not too many exons:
//     - Small introns: Select the best acceptable probe from a ~1kb region centered on the centre of the intron.
//     - Large introns: Select the best acceptable probe from each of ~1kb regions near the adjacent exons.
public class TargetGenes
{
    private static final ProbeSourceType PROBE_SOURCE = ProbeSourceType.GENE;

    private static final ProbeSelectCriteria EXON_PROBE_SELECT_CRITERIA = new ProbeSelectCriteria(
            new ProbeEvalCriteria(PROBE_QUALITY_REJECT, GENERAL_GC_TARGET, GENERAL_GC_TOLERANCE),
            ProbeSelectStrategy.MAX_QUALITY);
    private static final ProbeSelectCriteria CN_PROBE_SELECT_CRITERIA = new ProbeSelectCriteria(
            new ProbeEvalCriteria(PROBE_QUALITY_ACCEPT, CN_GC_TARGET, CN_GC_TOLERANCE),
            ProbeSelectStrategy.BEST_GC);

    private static final Logger LOGGER = LogManager.getLogger(TargetGenes.class);

    public static ProbeGenerationResult generateProbes(final String targetGeneFile, final EnsemblDataCache ensemblData,
            final ProbeGenerator probeGenerator)
    {
        LOGGER.info("Generating gene probes");

        List<GeneTranscriptId> geneIds = loadTargetGenesFile(targetGeneFile);
        List<GeneTranscript> genes = geneIds.stream()
                .map(gene -> loadGeneData(gene, ensemblData))
                .flatMap(Optional::stream)
                .toList();
        List<GeneRegion> geneRegions = genes.stream().flatMap(gene -> createGeneRegions(gene).stream()).toList();
        ProbeGenerationResult result = geneRegions.stream()
                .map(region -> generateProbes(region, probeGenerator))
                .reduce(new ProbeGenerationResult(), ProbeGenerationResult::add);

        LOGGER.info("Done generating gene probes");
        return result;
    }

    private record GeneTranscriptId(
            String geneName,
            String transcriptName
    )
    {
        public boolean canonical()
        {
            return transcriptName == null;
        }
    }

    private static List<GeneTranscriptId> loadTargetGenesFile(final String filePath)
    {
        try(DelimFileReader reader = new DelimFileReader(filePath))
        {
            int geneNameIdx = reader.getColumnIndex(FLD_GENE_NAME);
            int transNameIdx = reader.getColumnIndex(FLD_TRANS_NAME);

            List<GeneTranscriptId> genes = reader.stream().map(row ->
            {
                String geneName = row.get(geneNameIdx);
                String transcriptName = row.getOrNull(transNameIdx);
                return new GeneTranscriptId(geneName, transcriptName);
            }).toList();

            LOGGER.info("Loaded {} target genes from {}", genes.size(), filePath);
            return genes;
        }
    }

    public record GeneTranscript(
            GeneData gene,
            TranscriptData transcript
    )
    {
        @NotNull
        @Override
        public String toString()
        {
            return format("%s:%s", gene.GeneName, transcript.TransName);
        }
    }

    private static Optional<GeneTranscript> loadGeneData(final GeneTranscriptId geneId, final EnsemblDataCache ensemblData)
    {
        GeneData geneData = ensemblData.getGeneDataByName(geneId.geneName());
        if(geneData == null)
        {
            String error = format("Transcript gene for gene: %s not found", geneId.geneName());
            LOGGER.error(error);
            throw new RuntimeException(error);
        }

        TranscriptData transcriptData = geneId.canonical() ?
                ensemblData.getCanonicalTranscriptData(geneData.GeneId) :
                ensemblData.getTranscriptData(geneData.GeneId, geneId.transcriptName());
        if(transcriptData == null)
        {
            String error = format("gene(%s) transcript(%s) not found", geneId.geneName(), geneId.transcriptName());
            LOGGER.error(error);
            throw new RuntimeException(error);
        }

        if(transcriptData.nonCoding())
        {
            // User should add as a custom region instead.
            LOGGER.warn("gene({}) transcript({}) non-coding skipped", geneId.geneName(), geneId.transcriptName());
            return Optional.empty();
        }

        return Optional.of(new GeneTranscript(geneData, transcriptData));
    }

    private enum GeneRegionType
    {
        CODING,
        UTR,
        UP_STREAM,
        DOWN_STREAM,
        INTRONIC_LONG,
        INTRONIC_SHORT
    }

    private record GeneRegion(
            GeneTranscript gene,
            GeneRegionType type,
            ChrBaseRegion baseRegion
    )
    {
        public GeneRegion(GeneTranscript gene, GeneRegionType type, BaseRegion region)
        {
            this(gene, type, new ChrBaseRegion(gene.gene().Chromosome, region.start(), region.end()));
        }
    }

    private static List<GeneRegion> createGeneRegions(final GeneTranscript gene)
    {
        List<GeneRegion> regions = new ArrayList<>();

        GeneData geneData = gene.gene();
        TranscriptData transcriptData = gene.transcript();

        regions.add(new GeneRegion(
                gene,
                geneData.forwardStrand() ? GeneRegionType.UP_STREAM : GeneRegionType.DOWN_STREAM,
                new BaseRegion(
                        transcriptData.TransStart - GENE_UPDOWNSTREAM_GAP - GENE_UPDOWNSTREAM_REGION,
                        transcriptData.TransStart - GENE_UPDOWNSTREAM_GAP - 1)));

        int lastExonEnd = -1;

        boolean probeIntrons = transcriptData.exons().size() <= GENE_MAX_EXONS_TO_ADD_INTRON;

        for(ExonData exonData : transcriptData.exons())
        {
            if(probeIntrons && lastExonEnd != -1)
            {
                int intronLength = exonData.Start - lastExonEnd;

                if(intronLength >= GENE_MIN_INTRON_LENGTH)
                {
                    if(intronLength >= GENE_LONG_INTRON_LENGTH)
                    {
                        // GENE_LONG_INTRON_LENGTH should be large enough such that these two probes cannot overlap.
                        regions.add(new GeneRegion(
                                gene,
                                GeneRegionType.INTRONIC_LONG,
                                new BaseRegion(
                                        lastExonEnd + 1 + GENE_EXON_FLANK_GAP,
                                        lastExonEnd + GENE_EXON_FLANK_GAP + GENE_EXON_FLANK_REGION)));
                        regions.add(new GeneRegion(
                                gene,
                                GeneRegionType.INTRONIC_LONG,
                                new BaseRegion(
                                        exonData.Start - GENE_EXON_FLANK_GAP - GENE_EXON_FLANK_REGION,
                                        exonData.Start - GENE_EXON_FLANK_GAP - 1)));
                    }
                    else
                    {
                        int intronMid = (lastExonEnd + 1 + exonData.Start) / 2;
                        regions.add(new GeneRegion(
                                gene,
                                GeneRegionType.INTRONIC_SHORT,
                                new BaseRegion(
                                        intronMid - GENE_EXON_FLANK_REGION / 2,
                                        intronMid + GENE_EXON_FLANK_REGION / 2 - 1)));
                    }
                }
            }

            boolean isCoding = positionsOverlap(exonData.Start, exonData.End, transcriptData.CodingStart, transcriptData.CodingEnd);
            if(isCoding)
            {
                regions.add(new GeneRegion(
                        gene,
                        GeneRegionType.CODING,
                        new BaseRegion(
                                Math.max(exonData.Start, transcriptData.CodingStart),
                                Math.min(exonData.End, transcriptData.CodingEnd))));
            }
            else
            {
                int exonMid = (exonData.Start + exonData.End + 1) / 2;
                regions.add(new GeneRegion(
                        gene,
                        GeneRegionType.UTR,
                        new BaseRegion(exonMid - PROBE_LENGTH / 2, exonMid + PROBE_LENGTH / 2 - 1)));
            }

            lastExonEnd = exonData.End;
        }

        regions.add(new GeneRegion(
                gene,
                geneData.forwardStrand() ? GeneRegionType.DOWN_STREAM : GeneRegionType.UP_STREAM,
                new BaseRegion(
                        transcriptData.TransEnd + GENE_UPDOWNSTREAM_GAP,
                        transcriptData.TransEnd + GENE_UPDOWNSTREAM_GAP + GENE_UPDOWNSTREAM_REGION - 1)));

        regions.forEach(region -> LOGGER.trace("Gene region: {}", region));

        return regions;
    }

    private static ProbeGenerationResult generateProbes(final GeneRegion geneRegion, final ProbeGenerator probeGenerator)
    {
        LOGGER.trace("Generating probes for {}", geneRegion);

        ProbeSourceInfo source = createProbeSourceInfo(geneRegion);

        return switch(geneRegion.type())
        {
            case CODING, UTR -> probeGenerator.coverRegion(geneRegion.baseRegion(), source, EXON_PROBE_SELECT_CRITERIA);
            case UP_STREAM, DOWN_STREAM, INTRONIC_LONG, INTRONIC_SHORT ->
                    probeGenerator.coverOneSubregion(geneRegion.baseRegion(), source, CN_PROBE_SELECT_CRITERIA);
        };
    }

    private static ProbeSourceInfo createProbeSourceInfo(final GeneRegion geneRegion)
    {
        GeneData geneData = geneRegion.gene().gene();
        TranscriptData transcriptData = geneRegion.gene().transcript();
        String extraInfo = format("%s:%s:%s", geneData.GeneName, transcriptData.TransName, geneRegion.type().name());
        return new ProbeSourceInfo(PROBE_SOURCE, extraInfo);
    }
}
