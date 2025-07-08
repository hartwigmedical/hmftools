package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_TRANS_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_CANDIDATE_REGION_SIZE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_FLANKING_DISTANCE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_LONG_INTRON_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_MAX_CANDIDATE_PROBES;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_MAX_EXONS_TO_ADD_INTRON;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_MIN_INTRON_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;

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

public class TargetGenes
{
    private static final Logger LOGGER = LogManager.getLogger(TargetGenes.class);

    public List<ProbeCandidate> createProbeCandidates(final String targetGeneFile, final EnsemblDataCache ensemblData)
    {
        List<TargetGene> targetGenes = loadTargetGenesFile(targetGeneFile);
        List<GeneTranscript> genes = targetGenes.stream()
                .map(gene -> loadGeneData(gene, ensemblData))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .toList();
        List<GeneRegion> geneRegions = genes.stream().flatMap(gene -> createGeneRegions(gene).stream()).toList();
        List<ProbeCandidate> probeCandidates = geneRegions.stream().flatMap(region -> createProbeCandidates(region).stream()).toList();
        return probeCandidates;
    }

    private record TargetGene(
            String geneName,
            String transcriptName
    )
    {
        public boolean canonical()
        {
            return transcriptName == null;
        }
    }

    private static List<TargetGene> loadTargetGenesFile(final String filePath)
    {
        try(DelimFileReader reader = new DelimFileReader(filePath))
        {
            int geneNameIdx = reader.getColumnIndex(FLD_GENE_NAME);
            int transNameIdx = reader.getColumnIndex(FLD_TRANS_NAME);

            List<TargetGene> genes = reader.stream().map(row ->
            {
                String geneName = row.get(geneNameIdx);
                String transcriptName = row.getOrNull(transNameIdx);
                return new TargetGene(geneName, transcriptName);
            }).toList();

            LOGGER.info("Loaded {} target genes", genes.size());
            return genes;
        }
    }

    private static Optional<GeneTranscript> loadGeneData(final TargetGene gene, final EnsemblDataCache ensemblData)
    {
        GeneData geneData = ensemblData.getGeneDataByName(gene.geneName());
        if(geneData == null)
        {
            String error = format("Transcript gene for gene: %s not found", gene.geneName());
            GU_LOGGER.error(error);
            throw new RuntimeException(error);
        }

        TranscriptData transcriptData = gene.canonical() ?
                ensemblData.getCanonicalTranscriptData(geneData.GeneId) :
                ensemblData.getTranscriptData(geneData.GeneId, gene.transcriptName());
        if(transcriptData == null)
        {
            String error = format("gene(%s) transcript(%s) not found", gene.geneName(), gene.transcriptName());
            GU_LOGGER.error(error);
            throw new RuntimeException(error);
        }

        if(transcriptData.nonCoding())
        {
            // User should add as a custom region instead.
            GU_LOGGER.warn("gene({}) transcript({}) non-coding skipped", gene.geneName(), gene.transcriptName());
            return Optional.empty();
        }

        return Optional.of(new GeneTranscript(geneData, transcriptData));
    }

    private static List<GeneRegion> createGeneRegions(final GeneTranscript gene)
    {
        // Regions:
        //   - Coding: cover the full coding region of each exon
        //   - UTR: Add probe centered on each noncoding exon
        //   - Upstream/downstream: create 8 probes 1-2kb upstream/downstream, and choose the best probe.
        //   - Intronic:
        //     - For large introns, create 8 probes for 1-2kb on adjacent exons and choose the best probe for each.
        //     - For small introns, create 8 probes centered on the intron and choose the best probe for each.

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

        return regions;
    }

    private static List<ProbeCandidate> createProbeCandidates(final GeneRegion region) {
        // TODO
        switch(region.type())
        {
            case CODING:
            case UTR:
                // TODO: tile whole region
                break;
            case UP_STREAM:
            case DOWN_STREAM:
            case INTRONIC_LONG:
            case INTRONIC_SHORT:
                // TODO
                // create set of candidate probes within the region
                for(int i = 0; i < GENE_MAX_CANDIDATE_PROBES; ++i)
                {
                    int start = region.region().start() + i * PROBE_LENGTH;
                    int end = start + PROBE_LENGTH - 1; // end is inclusive

                    ProbeCandidate probe = createProbeCandidate(
                            new ChrBaseRegion(region.getGene().getGeneData().Chromosome, start, end), refGenomeInterface);

                    region.getProbeCandidates().add(probe);
                }
                break;
        }
    }
}
