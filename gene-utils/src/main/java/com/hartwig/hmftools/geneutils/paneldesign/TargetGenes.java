package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_TRANS_NAME;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_CANDIDATE_REGION_SIZE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_FLANKING_DISTANCE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_LONG_INTRON_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_MAX_CANDIDATE_PROBES;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_MAX_EXONS_TO_ADD_INTRON;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_MIN_INTRON_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;

import java.util.ArrayList;
import java.util.Collections;
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
//   - UTR: Probe centered on each noncoding exon.
//   - Upstream/downstream: Create several probes 1-2kb upstream/downstream, and choose the best probe.
//   - Intronic:
//     - For large introns, create several probes for 1-2kb on adjacent exons and choose the best probe for each.
//     - For small introns, create several probes centered on the intron and choose the best probe for each.
public class TargetGenes
{
    private static final ProbeSourceType PROBE_SOURCE = ProbeSourceType.GENE;

    private static final Logger LOGGER = LogManager.getLogger(TargetGenes.class);

    public static ProbeGenerationResult generateProbes(final String targetGeneFile, final EnsemblDataCache ensemblData,
            final ProbeEvaluator probeEvaluator)
    {
        LOGGER.info("Generating gene probes");

        List<GeneTranscriptId> geneIds = loadTargetGenesFile(targetGeneFile);
        List<GeneTranscript> genes = geneIds.stream()
                .map(gene -> loadGeneData(gene, ensemblData))
                .flatMap(Optional::stream)
                .toList();
        List<GeneRegion> geneRegions = genes.stream().flatMap(gene -> createGeneRegions(gene).stream()).toList();
        ProbeGenerationResult result = geneRegions.stream()
                .map(region -> generateProbe(region, probeEvaluator))
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
            BaseRegion region
    )
    {
        public ChrBaseRegion chrBaseRegion()
        {
            return new ChrBaseRegion(gene.gene().Chromosome, region.start(), region.end());
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
                        transcriptData.TransStart - GENE_FLANKING_DISTANCE - GENE_CANDIDATE_REGION_SIZE,
                        transcriptData.TransStart - GENE_FLANKING_DISTANCE - 1)));

        int lastExonEnd = -1;

        for(ExonData exonData : transcriptData.exons())
        {
            if(lastExonEnd != -1 && transcriptData.exons().size() <= GENE_MAX_EXONS_TO_ADD_INTRON)
            {
                int intronLength = exonData.Start - lastExonEnd;

                if(intronLength > GENE_MIN_INTRON_LENGTH)
                {
                    if(intronLength > GENE_LONG_INTRON_LENGTH)
                    {
                        regions.add(new GeneRegion(
                                gene,
                                GeneRegionType.INTRONIC_LONG,
                                new BaseRegion(
                                        lastExonEnd + 1 + GENE_FLANKING_DISTANCE,
                                        lastExonEnd + GENE_FLANKING_DISTANCE + GENE_CANDIDATE_REGION_SIZE)));

                        regions.add(new GeneRegion(
                                gene,
                                GeneRegionType.INTRONIC_LONG,
                                new BaseRegion(
                                        exonData.Start - GENE_FLANKING_DISTANCE - GENE_CANDIDATE_REGION_SIZE,
                                        exonData.Start - GENE_FLANKING_DISTANCE - 1)));
                    }
                    else
                    {
                        int intronMid = (lastExonEnd + 1 + exonData.Start) / 2;
                        regions.add(new GeneRegion(
                                gene,
                                GeneRegionType.INTRONIC_SHORT,
                                new BaseRegion(
                                        intronMid - GENE_CANDIDATE_REGION_SIZE / 2,
                                        intronMid + GENE_CANDIDATE_REGION_SIZE / 2 - 1)));
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
                        transcriptData.TransEnd + 1 + GENE_FLANKING_DISTANCE,
                        transcriptData.TransEnd + GENE_FLANKING_DISTANCE + GENE_CANDIDATE_REGION_SIZE)));

        regions.forEach(region -> LOGGER.trace("Gene region: {}", region));

        return regions;
    }

    private static ProbeGenerationResult generateProbe(final GeneRegion geneRegion, final ProbeEvaluator probeEvaluator)
    {
        return switch(geneRegion.type())
        {
            case CODING, UTR ->
                    RegionProbeTiling.fillRegionWithProbes(geneRegion.chrBaseRegion(), createProbeSourceInfo(geneRegion), probeEvaluator);
            case UP_STREAM, DOWN_STREAM, INTRONIC_LONG, INTRONIC_SHORT -> generateBestProbeInRegion(geneRegion, probeEvaluator);
        };
    }

    private static ProbeGenerationResult generateBestProbeInRegion(final GeneRegion geneRegion, final ProbeEvaluator probeEvaluator)
    {
        ProbeSourceInfo source = createProbeSourceInfo(geneRegion);

        // TODO: can probably look at all possible probes here instead, and put such a function in ProbeEvaluator
        List<CandidateProbe> candidates = RegionProbeTiling.tileBaseRegionsFrom(geneRegion.region().start())
                .limit(GENE_MAX_CANDIDATE_PROBES)
                .map(probeRegion ->
                {
                    ChrBaseRegion chrBaseRegion =
                            new ChrBaseRegion(geneRegion.gene().gene().Chromosome, probeRegion.start(), probeRegion.end());
                    return new CandidateProbe(source, chrBaseRegion, chrBaseRegion);
                })
                .toList();
        candidates.forEach(probe -> LOGGER.trace("Candidate probe: {}", probe));

        Optional<EvaluatedProbe> bestCandidate = probeEvaluator.selectBestProbe(candidates.stream());
        LOGGER.trace("{}: Best probe: {}", geneRegion, bestCandidate);

        ProbeGenerationResult result = bestCandidate
                .map(bestProbe -> new ProbeGenerationResult(List.of(bestProbe), Collections.emptyList()))
                .orElseGet(() ->
                        new ProbeGenerationResult(
                                Collections.emptyList(),
                                // TODO: rejection reason
                                List.of(new RejectedRegion(geneRegion.chrBaseRegion(), source, null))));
        return result;
    }

    private static ProbeSourceInfo createProbeSourceInfo(final GeneRegion geneRegion)
    {
        GeneData geneData = geneRegion.gene().gene();
        TranscriptData transcriptData = geneRegion.gene().transcript();
        String extraInfo = format("%s:%s:%s", geneData.GeneName, transcriptData.TransName, geneRegion.type().name());
        return new ProbeSourceInfo(PROBE_SOURCE, extraInfo);
    }
}
