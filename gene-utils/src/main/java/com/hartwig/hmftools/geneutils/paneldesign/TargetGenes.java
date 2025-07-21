package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.min;
import static java.lang.String.format;
import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_GC_OPTIMAL_TOLERANCE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_GC_TARGET;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CN_GC_TOLERANCE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENERAL_GC_TARGET;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENERAL_GC_TOLERANCE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_CN_QUALITY_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_CODING_REGION_EXPAND;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_EXON_FLANK_GAP;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_EXON_QUALITY_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_EXON_FLANK_REGION_MAX;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_EXON_FLANK_REGION_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_MAX_EXONS_TO_ADD_INTRON;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_UPDOWNSTREAM_GAP;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_UPDOWNSTREAM_REGION;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.regionCenteredAt;
import static com.hartwig.hmftools.geneutils.paneldesign.Utils.regionCentre;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

// TODO: unit test

// Probes covering (regions of) selected genes.
// Methodology for gene regions:
//   - Coding: Cover the full coding region of each exon, plus splice points.
//   - UTR: Select 1 probe within each noncoding exon.
//   - Upstream/downstream: Select the best acceptable probe from a 2kb region 1kb upstream/downstream.
//   - Exon flanks: Only when there are not too many exons:
//     - Small introns: Select the best acceptable probe from a 1kb region centered on the centre of the intron.
//     - Large introns: Select the best acceptable probe from each of 1-5kb regions near the adjacent exons.
public class TargetGenes
{
    private static final TargetMetadata.Type TARGET_REGION_TYPE = TargetMetadata.Type.GENE;

    private static final ProbeSelectCriteria EXON_PROBE_CRITERIA = new ProbeSelectCriteria(
            new ProbeEvaluator.Criteria(GENE_EXON_QUALITY_MIN, GENERAL_GC_TARGET, GENERAL_GC_TOLERANCE),
            new ProbeSelector.Strategy.MaxQuality());
    private static final ProbeSelectCriteria CN_PROBE_CRITERIA = new ProbeSelectCriteria(
            new ProbeEvaluator.Criteria(GENE_CN_QUALITY_MIN, CN_GC_TARGET, CN_GC_TOLERANCE),
            new ProbeSelector.Strategy.BestGc(CN_GC_OPTIMAL_TOLERANCE));

    private static final String FLD_INCLUDE_CODING = "IncludeCoding";
    private static final String FLD_INCLUDE_UTR = "IncludeUTR";
    private static final String FLD_INCLUDE_EXON_FLANK = "IncludeExonFlank";
    private static final String FLD_INCLUDE_UPSTREAM = "IncludeUpstream";
    private static final String FLD_INCLUDE_DOWNSTREAM = "IncludeDownstream";
    private static final String FLD_EXTRA_TRANSCRIPTS = "ExtraTransNames";

    private static final Logger LOGGER = LogManager.getLogger(TargetGenes.class);

    public static ProbeGenerationResult generateProbes(final String targetGeneFile, final EnsemblDataCache ensemblData,
            final ProbeGenerator probeGenerator)
    {
        LOGGER.info("Generating gene probes");

        List<GeneDefinition> geneDefs = loadTargetGenesFile(targetGeneFile);
        List<GeneTranscriptDefinition> geneTranscriptDefs = makeGeneTranscriptDefinitions(geneDefs);
        geneTranscriptDefs.forEach(gene -> LOGGER.debug("{}", gene));

        LOGGER.debug("Loading gene transcript data");
        List<GeneTranscriptData> geneTranscriptDatas = geneTranscriptDefs.stream()
                .map(gene -> loadGeneTranscriptData(gene, ensemblData))
                .flatMap(Optional::stream)
                .toList();

        LOGGER.debug("Creating gene target regions");
        List<GeneRegion> geneRegions = geneTranscriptDatas.stream().flatMap(gene -> createGeneRegions(gene).stream()).toList();

        // When generating probes, don't care about probe overlap. This is because:
        //   - Gene probes are generated first, so nothing is covered beforehand, and
        //   - Assume gene regions don't overlap, or if they do, it's small enough that the overlap is tolerable.

        LOGGER.debug("Generating probes");
        ProbeGenerationResult result = geneRegions.stream()
                .map(region -> generateProbes(region, probeGenerator))
                .reduce(new ProbeGenerationResult(), ProbeGenerationResult::add);

        LOGGER.info("Done generating gene probes");
        return result;
    }

    private record GeneOptions(
            boolean coding,
            boolean utr,
            boolean exonFlank,
            boolean upstream,
            boolean downstream
    )
    {
    }

    private record GeneDefinition(
            String geneName,
            GeneOptions options,
            List<String> extraTranscriptNames
    )
    {
        public List<String> transcriptNames()
        {
            // Always get the canonical transcript, plus any extra transcripts requested.
            List<String> names = new ArrayList<>();
            names.add(null);
            names.addAll(extraTranscriptNames);
            return names;
        }
    }

    private record GeneTranscriptDefinition(
            String geneName,
            // Null if canonical transcript.
            @Nullable String transcriptName,
            GeneOptions options
    )
    {
        public boolean canonical()
        {
            return transcriptName == null;
        }
    }

    private static List<GeneDefinition> loadTargetGenesFile(final String filePath)
    {
        LOGGER.debug("Loading target genes files: {}", filePath);

        try(DelimFileReader reader = new DelimFileReader(filePath))
        {
            List<GeneDefinition> genes = reader.stream().map(row ->
            {
                String geneName = row.get(FLD_GENE_NAME);
                boolean coding = row.getBoolean(FLD_INCLUDE_CODING);
                boolean utr = row.getBoolean(FLD_INCLUDE_UTR);
                boolean exonFlank = row.getBoolean(FLD_INCLUDE_EXON_FLANK);
                boolean upstream = row.getBoolean(FLD_INCLUDE_UPSTREAM);
                boolean downstream = row.getBoolean(FLD_INCLUDE_DOWNSTREAM);
                String extraTranscriptsStr = row.get(FLD_EXTRA_TRANSCRIPTS);
                List<String> extraTranscripts =
                        extraTranscriptsStr.isEmpty() ? emptyList() : Arrays.asList(extraTranscriptsStr.split(","));
                GeneOptions options = new GeneOptions(coding, utr, exonFlank, upstream, downstream);
                return new GeneDefinition(geneName, options, extraTranscripts);
            }).toList();

            LOGGER.info("Loaded {} target genes from {}", genes.size(), filePath);
            return genes;
        }
    }

    private static List<GeneTranscriptDefinition> makeGeneTranscriptDefinitions(final List<GeneDefinition> genes)
    {
        return genes.stream()
                .flatMap(gene -> gene.transcriptNames().stream()
                        .map(trans -> new GeneTranscriptDefinition(gene.geneName(), trans, gene.options()))
                )
                .toList();
    }

    private record GeneTranscriptData(
            GeneData gene,
            TranscriptData transcript,
            GeneOptions options
    )
    {
        @NotNull
        @Override
        public String toString()
        {
            return format("GeneTranscriptData[gene=%s, transcript=%s, options=%s]", gene.GeneName, transcript.TransName, options);
        }
    }

    private static Optional<GeneTranscriptData> loadGeneTranscriptData(final GeneTranscriptDefinition geneDef,
            final EnsemblDataCache ensemblData)
    {
        GeneData geneData = ensemblData.getGeneDataByName(geneDef.geneName());
        if(geneData == null)
        {
            throw new UserInputError(format("Gene not found: %s", geneDef.geneName()));
        }

        TranscriptData transcriptData = geneDef.canonical() ?
                ensemblData.getCanonicalTranscriptData(geneData.GeneId) :
                ensemblData.getTranscriptData(geneData.GeneId, geneDef.transcriptName());
        if(transcriptData == null)
        {
            throw new UserInputError(format("Gene not found: %s:%s", geneDef.geneName(), geneDef.transcriptName()));
        }

        if(transcriptData.nonCoding())
        {
            // User should add as a custom region instead.
            LOGGER.warn("Noncoding gene skipped: {}:{}", geneDef.geneName(), geneDef.transcriptName());
            return Optional.empty();
        }

        return Optional.of(new GeneTranscriptData(geneData, transcriptData, geneDef.options));
    }

    private enum GeneRegionType
    {
        CODING,
        UTR,
        UP_STREAM,
        DOWN_STREAM,
        EXON_FLANK
    }

    private record GeneRegion(
            GeneTranscriptData gene,
            GeneRegionType type,
            ChrBaseRegion region
    )
    {
        public GeneRegion(GeneTranscriptData gene, GeneRegionType type, BaseRegion region)
        {
            this(gene, type, ChrBaseRegion.from(gene.gene().Chromosome, region));
        }
    }

    private static List<GeneRegion> createGeneRegions(final GeneTranscriptData gene)
    {
        List<GeneRegion> regions = new ArrayList<>();

        GeneData geneData = gene.gene();
        TranscriptData transcriptData = gene.transcript();
        GeneOptions options = gene.options();

        if(geneData.forwardStrand() ? options.upstream() : options.downstream())
        {
            regions.add(new GeneRegion(
                    gene,
                    geneData.forwardStrand() ? GeneRegionType.UP_STREAM : GeneRegionType.DOWN_STREAM,
                    new BaseRegion(
                            transcriptData.TransStart - GENE_UPDOWNSTREAM_GAP - GENE_UPDOWNSTREAM_REGION,
                            transcriptData.TransStart - GENE_UPDOWNSTREAM_GAP - 1)));
        }

        int lastExonEnd = -1;
        boolean probeIntrons = options.exonFlank() && transcriptData.exons().size() <= GENE_MAX_EXONS_TO_ADD_INTRON;
        for(ExonData exonData : transcriptData.exons())
        {
            if(probeIntrons && lastExonEnd != -1)
            {
                int intronLength = exonData.Start - lastExonEnd;

                int minGap = 2 * GENE_EXON_FLANK_GAP;
                int available = intronLength - minGap;

                if(available >= GENE_EXON_FLANK_REGION_MIN)
                {
                    if(available >= GENE_EXON_FLANK_REGION_MIN * 3)
                    {
                        // Can fit 2 probes, with a gap in between.

                        // Remove space for the gap.
                        available -= GENE_EXON_FLANK_REGION_MIN;
                        // Allow the probe search regions to grow as the intron gets larger, up to a limit.
                        int regionSize = min(available / 2, GENE_EXON_FLANK_REGION_MAX);

                        regions.add(new GeneRegion(
                                gene,
                                GeneRegionType.EXON_FLANK,
                                new BaseRegion(
                                        lastExonEnd + 1 + GENE_EXON_FLANK_GAP,
                                        lastExonEnd + GENE_EXON_FLANK_GAP + regionSize)));
                        regions.add(new GeneRegion(
                                gene,
                                GeneRegionType.EXON_FLANK,
                                new BaseRegion(
                                        exonData.Start - GENE_EXON_FLANK_GAP - regionSize,
                                        exonData.Start - GENE_EXON_FLANK_GAP - 1)));
                    }
                    else
                    {
                        // Can fit 1 probe.
                        // TODO: doesn't expand region for larger introns - too restrictive?
                        int intronCentre = regionCentre(new BaseRegion(lastExonEnd, exonData.Start));
                        regions.add(new GeneRegion(
                                gene,
                                GeneRegionType.EXON_FLANK,
                                regionCenteredAt(intronCentre, GENE_EXON_FLANK_REGION_MIN)));
                    }
                }
            }

            boolean isCoding = positionsOverlap(exonData.Start, exonData.End, transcriptData.CodingStart, transcriptData.CodingEnd);
            if(isCoding)
            {
                if(options.coding())
                {
                    regions.add(new GeneRegion(
                            gene,
                            GeneRegionType.CODING,
                            new BaseRegion(
                                    Math.max(exonData.Start - GENE_CODING_REGION_EXPAND, transcriptData.CodingStart),
                                    min(exonData.End + GENE_CODING_REGION_EXPAND, transcriptData.CodingEnd))));
                }
            }
            else
            {
                if(options.utr())
                {
                    regions.add(new GeneRegion(
                            gene,
                            GeneRegionType.UTR,
                            new BaseRegion(exonData.Start, exonData.End)));
                }
            }

            lastExonEnd = exonData.End;
        }

        if(geneData.forwardStrand() ? options.downstream() : options.upstream())
        {
            regions.add(new GeneRegion(
                    gene,
                    geneData.forwardStrand() ? GeneRegionType.DOWN_STREAM : GeneRegionType.UP_STREAM,
                    new BaseRegion(
                            transcriptData.TransEnd + GENE_UPDOWNSTREAM_GAP,
                            transcriptData.TransEnd + GENE_UPDOWNSTREAM_GAP + GENE_UPDOWNSTREAM_REGION - 1)));
        }

        regions.forEach(region -> LOGGER.trace("Gene region: {}", region));

        return regions;
    }

    private static ProbeGenerationResult generateProbes(final GeneRegion geneRegion, final ProbeGenerator probeGenerator)
    {
        LOGGER.trace("Generating probes for {}", geneRegion);

        TargetMetadata metadata = createTargetMetadata(geneRegion);

        return switch(geneRegion.type())
        {
            case CODING ->
            {
                TargetRegion target = new TargetRegion(geneRegion.region(), metadata);
                CandidateProbeContext candidateContext = new CandidateProbeContext(target);
                yield probeGenerator.coverRegion(target.region(), candidateContext, EXON_PROBE_CRITERIA, null);
            }
            case UTR ->
            {
                BasePosition position = new BasePosition(geneRegion.region().chromosome(), regionCentre(geneRegion.region().baseRegion()));
                TargetRegion target = new TargetRegion(ChrBaseRegion.from(position), metadata);
                CandidateProbeContext candidateContext = new CandidateProbeContext(target);
                yield probeGenerator.coverPosition(position, candidateContext, EXON_PROBE_CRITERIA);
            }
            case UP_STREAM, DOWN_STREAM, EXON_FLANK ->
            {
                TargetRegion target = new TargetRegion(geneRegion.region(), metadata);
                CandidateProbeContext candidateContext = new CandidateProbeContext(target);
                yield probeGenerator.coverOneSubregion(target.region(), candidateContext, CN_PROBE_CRITERIA);
            }
        };
    }

    private static TargetMetadata createTargetMetadata(final GeneRegion geneRegion)
    {
        GeneData geneData = geneRegion.gene().gene();
        TranscriptData transcriptData = geneRegion.gene().transcript();
        String transcriptName = transcriptData.IsCanonical ? "canon" : transcriptData.TransName;
        String extraInfo = format("%s:%s:%s", geneData.GeneName, transcriptName, geneRegion.type().name());
        return new TargetMetadata(TARGET_REGION_TYPE, extraInfo);
    }
}
