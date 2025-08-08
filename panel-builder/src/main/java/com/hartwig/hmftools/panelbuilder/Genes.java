package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;
import static java.lang.String.join;
import static java.util.Collections.emptyList;

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

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Probes covering (regions of) selected genes.
// Methodology for gene regions:
//   - Coding: Cover the full coding region of each exon, plus splice points.
//   - UTR: Select 1 probe within each noncoding exon.
//   - Upstream/downstream: Select the best acceptable probe from a 2kb region 1kb upstream/downstream.
//   - Promoter: Cover the full region from the transcription start to 500b upstream.
//   - Exon flanks: Only when there are not too many exons:
//     - Small introns: Select the best acceptable probe from a 1kb region centered on the centre of the intron.
//     - Large introns: Select the best acceptable probe from each of 1-5kb regions near the adjacent exons.
public class Genes
{
    private static final TargetMetadata.Type TARGET_REGION_TYPE = TargetMetadata.Type.GENE;

    private static final ProbeSelectCriteria GENERAL_PROBE_CRITERIA = new ProbeSelectCriteria(
            new ProbeEvaluator.Criteria(GENE_GENERAL_QUALITY_MIN, GENE_GENERAL_GC_TARGET, GENE_GENERAL_GC_TOLERANCE),
            new ProbeSelector.Strategy.MaxQuality());
    private static final ProbeSelectCriteria CN_PROBE_CRITERIA = new ProbeSelectCriteria(
            new ProbeEvaluator.Criteria(GENE_CN_QUALITY_MIN, CN_GC_TARGET, CN_GC_TOLERANCE),
            new ProbeSelector.Strategy.BestGc(CN_GC_OPTIMAL_TOLERANCE));

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

        List<GeneDefinition> geneDefs = loadGenesFile(targetGeneFile);
        geneDefs.forEach(gene -> LOGGER.debug("{}", gene));

        LOGGER.debug("Loading gene transcript data");
        List<GeneTranscriptData> geneTranscriptDatas = geneDefs.stream()
                .map(gene -> loadGeneTranscriptData(gene, ensemblData))
                .flatMap(Optional::stream)
                .toList();

        // When generating probes, don't care about probe overlap. This is because:
        //   - Gene probes are generated first, so nothing is covered beforehand, and
        //   - Assume different genes don't overlap, or if they do, it's small enough that the overlap is tolerable, and
        //   - Within one gene, multiple transcripts are merged beforehand to avoid overlap.

        LOGGER.debug("Generating probes");
        ProbeGenerationResult result = new ProbeGenerationResult();
        ArrayList<GeneStats> geneStats = new ArrayList<>();
        for(GeneTranscriptData gene : geneTranscriptDatas)
        {
            List<GeneRegion> geneRegions = createGeneRegions(gene);

            ProbeGenerationResult geneResult = geneRegions.stream()
                    .map(geneRegion -> generateProbes(geneRegion, probeGenerator))
                    .reduce(new ProbeGenerationResult(), ProbeGenerationResult::add);
            result = result.add(geneResult);

            geneStats.add(new GeneStats(gene.gene().GeneName, geneResult.probes().size()));
        }
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
        public List<String> transcriptNames()
        {
            // Always get the canonical transcript, plus any extra transcripts requested.
            ArrayList<String> names = new ArrayList<>();
            names.add(null);
            names.addAll(extraTranscriptNames);
            return names;
        }
    }

    private static List<GeneDefinition> loadGenesFile(final String filePath)
    {
        LOGGER.debug("Loading genes file: {}", filePath);

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
                boolean promoter = row.getBoolean(FLD_INCLUDE_PROMOTER);
                String extraTranscriptsStr = row.get(FLD_EXTRA_TRANSCRIPTS);
                List<String> extraTranscripts =
                        extraTranscriptsStr.isEmpty() ? emptyList() : Arrays.asList(extraTranscriptsStr.split(","));
                GeneOptions options = new GeneOptions(coding, utr, exonFlank, upstream, downstream, promoter);
                return new GeneDefinition(geneName, options, extraTranscripts);
            }).toList();

            LOGGER.info("Loaded {} genes from {}", genes.size(), filePath);
            return genes;
        }
    }

    private record GeneTranscriptData(
            GeneData gene,
            List<TranscriptData> transcripts,
            GeneOptions options
    )
    {
    }

    private static Optional<GeneTranscriptData> loadGeneTranscriptData(final GeneDefinition geneDef,
            final EnsemblDataCache ensemblData)
    {
        GeneData geneData = ensemblData.getGeneDataByName(geneDef.geneName());
        if(geneData == null)
        {
            throw new UserInputError(format("Gene not found: %s", geneDef.geneName()));
        }

        List<TranscriptData> transcriptDatas = geneDef.transcriptNames().stream()
                .map(transName ->
                {
                    TranscriptData transcriptData = transName == null ?
                            ensemblData.getCanonicalTranscriptData(geneData.GeneId) :
                            ensemblData.getTranscriptData(geneData.GeneId, transName);
                    if(transcriptData == null)
                    {
                        throw new UserInputError(format("Gene not found: %s:%s", geneDef.geneName(), transName));
                    }
                    return transcriptData;
                })
                .filter(transcriptData ->
                {
                    if(transcriptData.nonCoding())
                    {
                        // User should add as a custom region instead.
                        LOGGER.warn("Noncoding gene skipped: {}:{}", geneDef.geneName(), transcriptData.TransName);
                        return false;
                    }
                    return true;
                })
                .toList();

        if(transcriptDatas.isEmpty())
        {
            return Optional.empty();
        }
        else
        {
            return Optional.of(new GeneTranscriptData(geneData, transcriptDatas, geneDef.options()));
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

    private static List<GeneRegion> createGeneRegions(final GeneTranscriptData gene)
    {
        // Methodology to handle multiple transcripts without producing overlapping probes:
        //   - Upstream/downstream: use the minimum/maximum start/end point of all transcripts.
        //   - Coding: superset of exons from all transcripts where any part of the exon is coding.
        //   - UTR: superset of exons from all transcripts where no part of the exon is coding.
        //   - Exon flanks: based on superset of exons from all transcripts.
        //   - Promoter: first/last exon in the superset of exons from all transcripts.

        ArrayList<GeneRegion> regions = new ArrayList<>();

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
                        // Can fit 1 probe.
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
                    regions.add(new GeneRegion(gene, GeneRegionType.UTR, exonRegion));
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
        public ArrayList<ExonData> Exons = new ArrayList<>();
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

    private static ProbeGenerationResult generateProbes(final GeneRegion geneRegion, final ProbeGenerator probeGenerator)
    {
        LOGGER.trace("Generating probes for {}", geneRegion);

        TargetMetadata metadata = createTargetMetadata(geneRegion);

        return switch(geneRegion.type())
        {
            case CODING, PROMOTER -> probeGenerator.coverRegion(geneRegion.region(), metadata, GENERAL_PROBE_CRITERIA, null);
            case UTR ->
            {
                BasePosition position = new BasePosition(geneRegion.region().chromosome(), regionCentre(geneRegion.region().baseRegion()));
                yield probeGenerator.coverPosition(position, metadata, GENERAL_PROBE_CRITERIA);
            }
            case UP_STREAM, DOWN_STREAM, EXON_FLANK -> probeGenerator.coverOneSubregion(geneRegion.region(), metadata, CN_PROBE_CRITERIA);
        };
    }

    private static TargetMetadata createTargetMetadata(final GeneRegion geneRegion)
    {
        GeneData geneData = geneRegion.gene().gene();
        List<String> transcriptNames = geneRegion.gene().transcripts().stream().map(Genes::transcriptDataName).toList();
        // If there are multiple transcripts, merge their names, since the region is determined based on 1 or more transcripts.
        String transcripts = join("/", transcriptNames);
        String extraInfo = format("%s:%s:%s", geneData.GeneName, transcripts, geneRegion.type().name());
        return new TargetMetadata(TARGET_REGION_TYPE, extraInfo);
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
}
