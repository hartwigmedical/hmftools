package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;
import static java.lang.String.join;
import static java.util.Collections.emptyList;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CN_GC_OPTIMAL_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CN_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CN_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_CN_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_CODING_REGION_EXPAND;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_EXON_FLANK_GAP;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_EXON_FLANK_REGION_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_EXON_FLANK_REGION_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_GENERAL_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_GENERAL_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_GENERAL_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_PROMOTER_REGION;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_UPDOWNSTREAM_GAP;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_UPDOWNSTREAM_REGION;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionCenteredAt;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionCentre;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionEndingAt;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionOverlapsOrAdjacent;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionStartingAt;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

// Probes covering (regions of) selected genes.
// Methodology for gene regions:
//   - Coding: Cover the full coding region of each exon, plus splice points.
//   - UTR: Select one probe within each noncoding exon.
//   - Upstream/downstream: Select the best acceptable probe from a 2kb region 1kb upstream/downstream.
//   - Promoter: Cover the full region from the transcription start to 500b upstream.
//   - Exon flanks: Only when there are not too many exons:
//     - Small introns: Select the best acceptable probe from a 1kb region centered on the centre of the intron.
//     - Large introns: Select the best acceptable probe from each of 1-5kb regions near the adjacent exons.
public class Genes
{
    private static final TargetMetadata.Type TARGET_TYPE = TargetMetadata.Type.GENE;

    private static final ProbeEvaluator.Criteria GENERAL_PROBE_CRITERIA = new ProbeEvaluator.Criteria(
            GENE_GENERAL_QUALITY_MIN, GENE_GENERAL_GC_TARGET, GENE_GENERAL_GC_TOLERANCE);
    private static final ProbeSelector.Strategy GENERAL_PROBE_SELECT = new ProbeSelector.Strategy.MaxQuality();
    private static final ProbeEvaluator.Criteria CN_PROBE_CRITERIA = new ProbeEvaluator.Criteria(
            GENE_CN_QUALITY_MIN, CN_GC_TARGET, CN_GC_TOLERANCE);
    private static final ProbeSelector.Strategy CN_PROBE_SELECT = new ProbeSelector.Strategy.BestGc(CN_GC_OPTIMAL_TOLERANCE);

    private static final String FLD_INCLUDE_CODING = "IncludeCoding";
    private static final String FLD_INCLUDE_UTR = "IncludeUTR";
    private static final String FLD_INCLUDE_EXON_FLANK = "IncludeExonFlank";
    private static final String FLD_INCLUDE_UPSTREAM = "IncludeUpstream";
    private static final String FLD_INCLUDE_DOWNSTREAM = "IncludeDownstream";
    private static final String FLD_INCLUDE_PROMOTER = "IncludePromoter";
    private static final String FLD_EXTRA_TRANSCRIPTS = "ExtraTransNames";

    private static final Logger LOGGER = LogManager.getLogger(Genes.class);

    public record ExtraOutput(
            List<GeneStats> geneStats
    )
    {
    }

    public record GeneStats(
            String geneName,
            int probeCount
    )
    {
    }

    public static ExtraOutput generateProbes(final String targetGeneFile, final EnsemblDataCache ensemblData,
            final ProbeGenerator probeGenerator, PanelData panelData)
    {
        LOGGER.info("Generating gene probes");

        List<GeneDefinition> geneDefinitions = loadGenesFile(targetGeneFile);
        geneDefinitions.forEach(gene -> LOGGER.debug("{}", gene));

        LOGGER.debug("Loading gene transcript data");
        List<GeneTranscriptData> geneTranscriptDatas = loadGeneTranscriptDatas(geneDefinitions, ensemblData);

        // When generating probes, don't check probe overlap between genes. This is because:
        //   - Assume different genes don't overlap, or if they do, it's small enough that the overlap is tolerable, and
        //   - Within one gene, multiple transcripts are merged beforehand to avoid subregion overlap.

        LOGGER.debug("Generating probes");
        ProbeGenerationResult result = generateProbes(geneTranscriptDatas, probeGenerator, panelData);

        List<GeneStats> geneStats = computeGeneStats(result, geneTranscriptDatas);
        ExtraOutput extraOutput = new ExtraOutput(geneStats);

        panelData.addResult(result);

        LOGGER.info("Done generating gene probes");

        return extraOutput;
    }

    private record GeneOptions(
            boolean coding,
            boolean utr,
            boolean exonFlank,
            boolean upstream,
            boolean downstream,
            boolean promoter
    )
    {
    }

    private record GeneDefinition(
            String geneName,
            GeneOptions options,
            List<String> extraTranscriptNames
    )
    {
    }

    private static List<GeneDefinition> loadGenesFile(final String filePath)
    {
        LOGGER.debug("Loading genes file: {}", filePath);

        try(DelimFileReader reader = new DelimFileReader(filePath))
        {
            List<GeneDefinition> genes = reader.stream().map(row ->
            {
                String geneName = row.get(FLD_GENE_NAME);
                if(geneName.isBlank())
                {
                    throw new UserInputError("Gene name cannot be blank: " + geneName);
                }

                boolean coding = row.getBoolean(FLD_INCLUDE_CODING);
                boolean utr = row.getBoolean(FLD_INCLUDE_UTR);
                boolean exonFlank = row.getBoolean(FLD_INCLUDE_EXON_FLANK);
                boolean upstream = row.getBoolean(FLD_INCLUDE_UPSTREAM);
                boolean downstream = row.getBoolean(FLD_INCLUDE_DOWNSTREAM);
                boolean promoter = row.getBoolean(FLD_INCLUDE_PROMOTER);
                String extraTranscriptsStr = row.getOrNull(FLD_EXTRA_TRANSCRIPTS);
                List<String> extraTranscripts = parseGeneExtraTranscripts(extraTranscriptsStr);

                GeneOptions options = new GeneOptions(coding, utr, exonFlank, upstream, downstream, promoter);
                return new GeneDefinition(geneName, options, extraTranscripts);
            }).toList();

            LOGGER.debug("Loaded {} genes from {}", genes.size(), filePath);
            return genes;
        }
    }

    private static List<String> parseGeneExtraTranscripts(@Nullable final String field)
    {
        if(field == null || field.isEmpty())
        {
            return emptyList();
        }
        else
        {
            List<String> transcriptNames = Arrays.asList(field.strip().split(","));
            for(String transcriptName : transcriptNames)
            {
                if(transcriptName.isBlank())
                {
                    throw new UserInputError("Transcript name cannot be blank: " + transcriptName);
                }
            }
            return transcriptNames;
        }
    }

    private record GeneTranscriptData(
            GeneData gene,
            List<TranscriptData> transcripts,
            GeneOptions options
    )
    {
    }

    private static List<GeneTranscriptData> loadGeneTranscriptDatas(final List<GeneDefinition> geneDefinitions,
            final EnsemblDataCache ensemblData)
    {
        List<GeneTranscriptData> geneTranscriptDatas = new ArrayList<>();
        boolean error = false;
        for(GeneDefinition geneDef : geneDefinitions)
        {
            Optional<GeneTranscriptData> geneTranscriptData = loadGeneTranscriptDatas(geneDef, ensemblData);
            if(geneTranscriptData.isPresent())
            {
                geneTranscriptDatas.add(geneTranscriptData.get());
            }
            else
            {
                // Don't immediately fail, so we can log all the errors for the user. Only fail at the end of all the loading.
                error = true;
            }
        }
        if(error)
        {
            throw new UserInputError("Invalid genes (see error logs for details)");
        }
        else
        {
            return geneTranscriptDatas;
        }
    }

    private static Optional<GeneTranscriptData> loadGeneTranscriptDatas(final GeneDefinition geneDef,
            final EnsemblDataCache ensemblData)
    {
        GeneData geneData = ensemblData.getGeneDataByName(geneDef.geneName());
        if(geneData == null)
        {
            LOGGER.error("Gene not found: {}", geneDef.geneName());
            return Optional.empty();
        }

        List<TranscriptData> transcriptDatas = loadGeneTranscripts(geneData, geneDef.extraTranscriptNames(), ensemblData);
        if(transcriptDatas.isEmpty())
        {
            // Error occurred loading transcripts.
            return Optional.empty();
        }
        else
        {
            return Optional.of(new GeneTranscriptData(geneData, transcriptDatas, geneDef.options()));
        }
    }

    private static List<TranscriptData> loadGeneTranscripts(final GeneData geneData, final List<String> extraTranscriptNames,
            final EnsemblDataCache ensemblData)
    {
        List<String> transcriptNames = new ArrayList<>();
        // Null indicates the canonical transcript.
        transcriptNames.add(null);
        transcriptNames.addAll(extraTranscriptNames);

        List<TranscriptData> transcriptDatas = new ArrayList<>();
        boolean error = false;
        for(String transName : transcriptNames)
        {
            TranscriptData transcriptData = transName == null ?
                    ensemblData.getCanonicalTranscriptData(geneData.GeneId) :
                    ensemblData.getTranscriptData(geneData.GeneId, transName);
            if(transcriptData == null)
            {
                LOGGER.error("Gene transcript not found: {}:{}", geneData.GeneName, transName);
                error = true;
            }
            else if(transcriptData.nonCoding())
            {
                // The only support for noncoding regions we have is a single probe for UTR, which is probably not what the user wants.
                LOGGER.error("Noncoding gene transcript, add as custom region instead: {}:{}", geneData.GeneName, transcriptData.TransName);
                error = true;
            }
            else
            {
                transcriptDatas.add(transcriptData);
            }
        }

        if(error)
        {
            return emptyList();
        }
        else
        {
            return transcriptDatas;
        }
    }

    private enum GeneRegionType
    {
        CODING,
        UTR,
        UP_STREAM,
        DOWN_STREAM,
        EXON_FLANK,
        PROMOTER
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

    private static ProbeGenerationResult generateProbes(final List<GeneTranscriptData> genes, final ProbeGenerator probeGenerator,
            final PanelCoverage coverage)
    {
        Stream<ProbeGenerationSpec> probeGenerationSpecs = genes.stream()
                .flatMap(gene -> createGeneRegions(gene).stream().map(Genes::createProbeGenerationSpec));
        return probeGenerator.generateBatch(probeGenerationSpecs, coverage);
    }

    private static List<GeneRegion> createGeneRegions(final GeneTranscriptData gene)
    {
        // Methodology to handle multiple transcripts without producing overlapping probes:
        //   - Upstream/downstream: use the minimum/maximum start/end point of all transcripts.
        //   - Coding: superset of exons from all transcripts where any part of the exon is coding.
        //   - UTR: superset of exons from all transcripts where no part of the exon is coding.
        //   - Exon flanks: based on superset of exons from all transcripts.
        //   - Promoter: first/last exon in the superset of exons from all transcripts.

        List<GeneRegion> regions = new ArrayList<>();

        GeneData geneData = gene.gene();
        List<TranscriptData> transcripts = gene.transcripts();
        GeneOptions options = gene.options();

        if(geneData.forwardStrand() ? options.upstream() : options.downstream())
        {
            int minTransStart = transcripts.stream().mapToInt(trans -> trans.TransStart).min().orElseThrow();
            regions.add(new GeneRegion(
                    gene,
                    geneData.forwardStrand() ? GeneRegionType.UP_STREAM : GeneRegionType.DOWN_STREAM,
                    regionEndingAt(minTransStart - 1 - GENE_UPDOWNSTREAM_GAP, GENE_UPDOWNSTREAM_REGION)));
        }

        List<MergedExonRegion> mergedExons = mergeExons(transcripts);

        if(options.promoter() && !mergedExons.isEmpty())
        {
            int transStart = geneData.forwardStrand()
                    ? transcripts.stream().mapToInt(trans -> trans.TransStart).min().orElseThrow()
                    : transcripts.stream().mapToInt(trans -> trans.TransEnd).max().orElseThrow();
            BaseRegion region = geneData.forwardStrand()
                    ? regionEndingAt(transStart - 1, GENE_PROMOTER_REGION)
                    : regionStartingAt(transStart + 1, GENE_PROMOTER_REGION);
            regions.add(new GeneRegion(gene, GeneRegionType.PROMOTER, region));
        }

        int lastExonEnd = 0;
        for(MergedExonRegion mergedExon : mergedExons)
        {
            BaseRegion exonRegion = mergedExon.Region;
            if(options.exonFlank() && lastExonEnd > 0)
            {
                int intronLength = exonRegion.start() - lastExonEnd;

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
                                regionStartingAt(lastExonEnd + 1 + GENE_EXON_FLANK_GAP, regionSize)));
                        regions.add(new GeneRegion(
                                gene,
                                GeneRegionType.EXON_FLANK,
                                regionEndingAt(exonRegion.start() - 1 - GENE_EXON_FLANK_GAP, regionSize)));
                    }
                    else
                    {
                        // Can fit one probe.
                        int intronCentre = regionCentre(new BaseRegion(lastExonEnd, exonRegion.start()));
                        regions.add(new GeneRegion(
                                gene,
                                GeneRegionType.EXON_FLANK,
                                regionCenteredAt(intronCentre, GENE_EXON_FLANK_REGION_MIN)));
                    }
                }
            }

            if(mergedExon.IsCoding)
            {
                if(options.coding())
                {
                    regions.add(new GeneRegion(
                            gene,
                            GeneRegionType.CODING,
                            new BaseRegion(
                                    max(exonRegion.start() - GENE_CODING_REGION_EXPAND, mergedExon.CodingStart),
                                    min(exonRegion.end() + GENE_CODING_REGION_EXPAND, mergedExon.CodingEnd))));
                }
            }
            else
            {
                if(options.utr())
                {
                    // Just one probe approximately at the centre of the exon.
                    // Noncoding exons can be long and could produce a lot of probes which are not useful.
                    int centre = regionCentre(exonRegion);
                    regions.add(new GeneRegion(gene, GeneRegionType.UTR, new BaseRegion(centre, centre)));
                }
            }

            lastExonEnd = mergedExon.Region.end();
        }

        if(geneData.forwardStrand() ? options.downstream() : options.upstream())
        {
            int maxTransEnd = transcripts.stream().mapToInt(trans -> trans.TransEnd).max().orElseThrow();
            regions.add(new GeneRegion(
                    gene,
                    geneData.forwardStrand() ? GeneRegionType.DOWN_STREAM : GeneRegionType.UP_STREAM,
                    regionStartingAt(maxTransEnd + GENE_UPDOWNSTREAM_GAP, GENE_UPDOWNSTREAM_REGION)));
        }

        regions.forEach(region -> LOGGER.trace("Gene region: {}", region));

        return regions;
    }

    private static class MergedExonRegion
    {
        public BaseRegion Region;
        public boolean IsCoding = false;
        public int CodingStart = Integer.MAX_VALUE;
        public int CodingEnd = Integer.MIN_VALUE;
        public List<ExonData> Exons = new ArrayList<>();
    }

    private static List<MergedExonRegion> mergeExons(final List<TranscriptData> transcripts)
    {
        List<MergedExonRegion> mergedExons = new ArrayList<>();
        // Standard region merging algorithm...
        List<ImmutablePair<TranscriptData, ExonData>> sortedExons = transcripts.stream()
                .flatMap(trans -> trans.exons().stream().map(exon -> new ImmutablePair<>(trans, exon)))
                .sorted(Comparator.comparing(pair -> pair.getRight().Start))
                .toList();
        for(Pair<TranscriptData, ExonData> pair : sortedExons)
        {
            TranscriptData transcript = pair.getLeft();
            ExonData exon = pair.getRight();

            BaseRegion codingRegion = new BaseRegion(transcript.CodingStart, transcript.CodingEnd);
            BaseRegion exonRegion = new BaseRegion(exon.Start, exon.End);
            boolean isCoding = codingRegion.overlaps(exonRegion);

            mergedExons.stream()
                    .filter(region -> regionOverlapsOrAdjacent(region.Region, exonRegion))
                    .findFirst()
                    .ifPresentOrElse(merged ->
                            {
                                if(exon.End > merged.Region.end())
                                {
                                    // Don't mutate in place because we borrowed the object from the exon data.
                                    merged.Region = new BaseRegion(merged.Region.start(), exon.End);
                                }
                                if(isCoding)
                                {
                                    merged.IsCoding = true;
                                    merged.CodingStart = min(merged.CodingStart, codingRegion.start());
                                    merged.CodingEnd = max(merged.CodingEnd, codingRegion.end());
                                }
                                merged.Exons.add(exon);
                            },
                            () ->
                            {
                                MergedExonRegion merged = new MergedExonRegion();
                                merged.Region = exonRegion;
                                merged.IsCoding = isCoding;
                                if(isCoding)
                                {
                                    merged.CodingStart = codingRegion.start();
                                    merged.CodingEnd = codingRegion.end();
                                }
                                merged.Exons.add(exon);
                                mergedExons.add(merged);
                            });
        }
        return mergedExons;
    }

    private static ProbeGenerationSpec createProbeGenerationSpec(final GeneRegion geneRegion)
    {
        LOGGER.trace("Generating probes for {}", geneRegion);

        TargetMetadata metadata = createTargetMetadata(geneRegion);

        return switch(geneRegion.type())
        {
            case CODING, UTR, PROMOTER ->
                    new ProbeGenerationSpec.CoverRegion(geneRegion.region(), metadata, GENERAL_PROBE_CRITERIA, GENERAL_PROBE_SELECT);
            case UP_STREAM, DOWN_STREAM, EXON_FLANK ->
                    new ProbeGenerationSpec.CoverOneSubregion(geneRegion.region(), metadata, CN_PROBE_CRITERIA, CN_PROBE_SELECT);
        };
    }

    private static TargetMetadata createTargetMetadata(final GeneRegion geneRegion)
    {
        GeneData geneData = geneRegion.gene().gene();
        List<String> transcriptNames = geneRegion.gene().transcripts().stream().map(Genes::transcriptDataName).toList();
        // If there are multiple transcripts, merge their names, since the region is determined based on one or more transcripts.
        String transcripts = join("/", transcriptNames);
        String extraInfo = format("%s:%s:%s", geneData.GeneName, transcripts, geneRegion.type().name());
        // Store the gene region as extra data so it can be recovered later to generate summary statistics.
        return new TargetMetadata(TARGET_TYPE, extraInfo, geneRegion);
    }

    private static String transcriptDataName(final TranscriptData transcriptData)
    {
        if(transcriptData.IsCanonical)
        {
            return "canon";
        }
        else
        {
            return transcriptData.TransName;
        }
    }

    private static List<GeneStats> computeGeneStats(final ProbeGenerationResult result, List<GeneTranscriptData> genes)
    {
        Map<String, Integer> geneProbeCounts = new HashMap<>();
        for(GeneTranscriptData gene : genes)
        {
            geneProbeCounts.put(gene.gene.GeneName, 0);
        }
        for(Probe probe : result.probes())
        {
            GeneRegion geneRegion = (GeneRegion) requireNonNull(probe.metadata().extraData());
            geneProbeCounts.merge(geneRegion.gene().gene().GeneName, 1, Integer::sum);
        }
        return geneProbeCounts.entrySet().stream()
                // Sort by gene name to ensure deterministic output.
                .sorted(Map.Entry.comparingByKey())
                .map(entry -> new GeneStats(entry.getKey(), entry.getValue()))
                .toList();
    }
}
