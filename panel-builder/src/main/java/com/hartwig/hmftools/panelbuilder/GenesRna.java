package com.hartwig.hmftools.panelbuilder;

import static java.lang.String.format;
import static java.lang.String.join;
import static java.util.Collections.emptyList;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.panelbuilder.GeneUtils.mergeExons;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_RNA_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_RNA_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_RNA_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.Utils.findDuplicates;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.OptionalInt;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

// Probes covering the transcribed (exonic) sequence of selected genes, for RNA panels.
// Coding sequence is always covered; the 5' and 3' UTRs are optional per gene. Probes are exon-aware: they only cover exonic bases, tiled
// within each exon and flush to splice junctions, with short exons padded across the junction into the adjacent exon. This is all driven off
// a per-gene RegionMapping (the merged exons as a contiguous probe-space) and RnaProbeGenerator.
public class GenesRna
{
    private static final TargetMetadata.Type TARGET_TYPE = TargetMetadata.Type.GENE_RNA;

    private static final ProbeEvaluator.Criteria PROBE_CRITERIA = new ProbeEvaluator.Criteria(
            GENE_RNA_QUALITY_MIN, GENE_RNA_GC_TARGET, GENE_RNA_GC_TOLERANCE);
    private static final ProbeSelector.Strategy PROBE_SELECT = new ProbeSelector.Strategy.MaxQuality();

    private static final Logger LOGGER = LogManager.getLogger(GenesRna.class);

    private enum Column
    {
        GeneName,
        Include5UTR,
        Include3UTR,
        TransNames
    }

    public record ExtraOutput(
            List<GeneStats> geneStats
    )
    {
    }

    public static ExtraOutput generateProbes(final String rnaGeneFile, final EnsemblDataCache ensemblData,
            final RnaProbeGenerator probeGenerator, PanelData panelData)
    {
        LOGGER.info("Generating gene RNA probes");

        List<GeneDefinition> geneDefinitions = loadGenesFile(rnaGeneFile);
        geneDefinitions.forEach(gene -> LOGGER.debug("{}", gene));

        LOGGER.debug("Loading gene transcript data");
        List<GeneTranscriptData> geneTranscriptDatas = loadGeneTranscriptDatas(geneDefinitions, ensemblData);
        checkNoDuplicateGenes(geneTranscriptDatas);

        LOGGER.debug("Generating probes");
        ProbeGenerationResult result = generateProbes(geneTranscriptDatas, probeGenerator, panelData);

        List<GeneStats> geneStats = computeGeneStats(result, geneTranscriptDatas);
        ExtraOutput extraOutput = new ExtraOutput(geneStats);

        LOGGER.info("Done generating gene RNA probes");

        return extraOutput;
    }

    record GeneOptions(
            boolean utr5,
            boolean utr3
    )
    {
    }

    private record GeneDefinition(
            String geneName,
            GeneOptions options,
            // Empty means all transcripts of the gene (merged); otherwise the specific transcripts to use (ENST or NM ids).
            List<String> transcriptNames
    )
    {
    }

    private static List<GeneDefinition> loadGenesFile(final String filePath)
    {
        LOGGER.debug("Loading genes RNA file: {}", filePath);

        try(DelimFileReader reader = new DelimFileReader(filePath))
        {
            List<GeneDefinition> genes = reader.stream().map(row ->
            {
                String geneName = row.get(Column.GeneName);
                if(geneName.isBlank())
                {
                    throw new UserInputError("Gene name cannot be blank");
                }
                boolean utr5 = row.getBoolean(Column.Include5UTR);
                boolean utr3 = row.getBoolean(Column.Include3UTR);
                List<String> transcripts = parseTranscripts(row.getStringOrNull(Column.TransNames));
                return new GeneDefinition(geneName, new GeneOptions(utr5, utr3), transcripts);
            }).toList();

            LOGGER.debug("Loaded {} genes from {}", genes.size(), filePath);
            return genes;
        }
    }

    private static List<String> parseTranscripts(@Nullable final String field)
    {
        if(field == null || field.isEmpty())
        {
            return emptyList();
        }
        List<String> transcriptNames = Arrays.asList(field.strip().split(","));
        for(String transcriptName : transcriptNames)
        {
            if(transcriptName.isBlank())
            {
                throw new UserInputError("Transcript name cannot be blank");
            }
        }
        return transcriptNames;
    }

    private record GeneTranscriptData(
            GeneData gene,
            List<TranscriptData> transcripts,
            GeneOptions options,
            // Whether the user requested a specific subset of transcripts (rather than the default of all transcripts).
            boolean specificTranscripts
    )
    {
    }

    private static void checkNoDuplicateGenes(final List<GeneTranscriptData> genes)
    {
        List<String> names = genes.stream().map(gene -> gene.gene().GeneName).toList();
        List<String> duplicates = findDuplicates(names);
        if(!duplicates.isEmpty())
        {
            duplicates.forEach(name -> LOGGER.error("Duplicate gene: {}", name));
            throw new UserInputError("Duplicate genes");
        }
    }

    private static List<GeneTranscriptData> loadGeneTranscriptDatas(final List<GeneDefinition> geneDefinitions,
            final EnsemblDataCache ensemblData)
    {
        List<GeneTranscriptData> geneTranscriptDatas = new ArrayList<>();
        boolean error = false;
        for(GeneDefinition geneDef : geneDefinitions)
        {
            // Don't fail immediately, so all errors are logged for the user; fail only after loading everything.
            Optional<GeneTranscriptData> geneTranscriptData = loadGeneTranscriptData(geneDef, ensemblData);
            if(geneTranscriptData.isPresent())
            {
                geneTranscriptDatas.add(geneTranscriptData.get());
            }
            else
            {
                error = true;
            }
        }
        if(error)
        {
            throw new UserInputError("Invalid genes (see error logs for details)");
        }
        return geneTranscriptDatas;
    }

    private static Optional<GeneTranscriptData> loadGeneTranscriptData(final GeneDefinition geneDef, final EnsemblDataCache ensemblData)
    {
        GeneData geneData = ensemblData.getGeneDataByName(geneDef.geneName());
        if(geneData == null)
        {
            LOGGER.error("Gene not found: {}", geneDef.geneName());
            return Optional.empty();
        }

        List<TranscriptData> transcripts = resolveTranscripts(geneData, geneDef.transcriptNames(), ensemblData);
        if(transcripts.isEmpty())
        {
            // No transcripts, or an error resolving them.
            LOGGER.error("No transcripts resolved for gene: {}", geneDef.geneName());
            return Optional.empty();
        }
        boolean specificTranscripts = !geneDef.transcriptNames().isEmpty();
        return Optional.of(new GeneTranscriptData(geneData, transcripts, geneDef.options(), specificTranscripts));
    }

    // Resolves the transcripts to use: all of the gene's transcripts if none were specified, otherwise the specified Ensembl transcripts.
    private static List<TranscriptData> resolveTranscripts(final GeneData geneData, final List<String> transcriptNames,
            final EnsemblDataCache ensemblData)
    {
        List<TranscriptData> allTranscripts = ensemblData.getTranscripts(geneData.GeneId);
        if(allTranscripts == null || allTranscripts.isEmpty())
        {
            return emptyList();
        }

        if(transcriptNames.isEmpty())
        {
            return allTranscripts;
        }

        List<TranscriptData> resolved = new ArrayList<>();
        boolean error = false;
        for(String transcriptName : transcriptNames)
        {
            TranscriptData transcript = resolveTranscript(geneData, allTranscripts, transcriptName);
            if(transcript == null)
            {
                error = true;
            }
            else
            {
                resolved.add(transcript);
            }
        }
        return error ? emptyList() : resolved;
    }

    // Resolves a single Ensembl transcript by its TransName. RefSeq (NM) resolution is not yet supported (needs additional validation of the
    // non-1:1 Ensembl<->RefSeq mapping), so any non-Ensembl id is reported as not found.
    // TODO: validate the Ensembl<->RefSeq mapping (RefSeqId may be null or multi-valued) and re-enable NM resolution via a RefSeqId fallback
    //  with clear not-found / ambiguous errors.
    @Nullable
    static TranscriptData resolveTranscript(final GeneData geneData, final List<TranscriptData> allTranscripts,
            final String transcriptName)
    {
        TranscriptData match = allTranscripts.stream()
                .filter(transcript -> transcript.TransName.equals(transcriptName))
                .findFirst()
                .orElse(null);
        if(match == null)
        {
            LOGGER.error("Gene transcript not found: {}:{}", geneData.GeneName, transcriptName);
        }
        return match;
    }

    enum RnaRegionType
    {
        CODING,
        UTR_5,
        UTR_3
    }

    // A single-exon target range in probe-space [spaceStart, spaceEnd).
    record RnaTarget(
            RnaRegionType type,
            int spaceStart,
            int spaceEnd
    )
    {
    }

    record GeneTargets(
            RegionMapping mapping,
            List<RnaTarget> targets
    )
    {
    }

    private static ProbeGenerationResult generateProbes(final List<GeneTranscriptData> genes, final RnaProbeGenerator probeGenerator,
            PanelData panelData)
    {
        ProbeGenerationResult total = new ProbeGenerationResult();
        for(GeneTranscriptData gene : genes)
        {
            GeneTargets geneTargets = createTargets(gene.gene(), gene.transcripts(), gene.options());
            if(geneTargets.targets().isEmpty())
            {
                LOGGER.warn("No features produced for gene {}", gene.gene().GeneName);
            }
            for(RnaTarget target : geneTargets.targets())
            {
                TargetMetadata metadata = createTargetMetadata(gene, target);
                ProbeGenerationResult result = probeGenerator.coverExonRange(
                        geneTargets.mapping(), target.spaceStart(), target.spaceEnd(), metadata, PROBE_CRITERIA, PROBE_SELECT);
                // coverExonRange does not add the candidate target region (it mirrors coverUncoveredRegion); add it here.
                ChrBaseRegion targetRegion = geneTargets.mapping().toGenomeRegions(target.spaceStart(), target.spaceEnd()).get(0);
                result = result.add(new ProbeGenerationResult(
                        emptyList(), List.of(new TargetRegion(targetRegion, metadata)), emptyList()));
                panelData.addResult(result);
                total = total.add(result);
            }
        }
        return total;
    }

    // Computes the exon-aware target ranges for a gene: the coding sequence (always) and, if requested, the 5' and 3' UTRs. All exons of the
    // selected transcripts are merged into one RegionMapping (so short-exon probes can pad across junctions). Each merged exon is one target,
    // classified whole: an exon with any coding bases is treated as coding, otherwise it is a fully noncoding (UTR) exon, 5' or 3' by its
    // position relative to the coding span (strand-aware). Both edges of each target are therefore true splice junctions.
    // TODO: reconsider classifying a part-coding exon as entirely coding. A very long exon with only a few coding bases would be tiled fully
    //  as coding, covering a lot of UTR sequence that cannot be excluded and is not attributed to the 5'/3' UTR features.
    static GeneTargets createTargets(final GeneData geneData, final List<TranscriptData> transcripts, final GeneOptions options)
    {
        List<GeneUtils.MergedExonRegion> mergedExons = mergeExons(transcripts);

        List<ChrBaseRegion> exonRegions = mergedExons.stream()
                .map(merged -> ChrBaseRegion.from(geneData.Chromosome, merged.Region))
                .toList();
        RegionMapping mapping = new RegionMapping(exonRegions);

        // Gene-wide coding start (across coding transcripts); used to classify fully noncoding exons as 5' or 3'.
        OptionalInt codingStart = transcripts.stream().filter(t -> !t.nonCoding()).mapToInt(t -> t.CodingStart).min();
        if(codingStart.isEmpty())
        {
            // Noncoding gene: no coding to anchor the 5'/3' UTR split, so nothing is produced.
            LOGGER.warn("Gene {} has no coding transcripts; no RNA probes produced", geneData.GeneName);
            return new GeneTargets(mapping, emptyList());
        }

        List<RnaTarget> targets = new ArrayList<>();
        for(GeneUtils.MergedExonRegion mergedExon : mergedExons)
        {
            BaseRegion exon = mergedExon.Region;
            RnaRegionType type;
            if(mergedExon.IsCoding)
            {
                type = RnaRegionType.CODING;
            }
            else
            {
                // Fully noncoding exon: below the coding start is 5' UTR on the forward strand, 3' UTR on the reverse strand (and vice versa).
                boolean belowCoding = exon.end() < codingStart.getAsInt();
                type = belowCoding == geneData.forwardStrand() ? RnaRegionType.UTR_5 : RnaRegionType.UTR_3;
            }
            addTarget(targets, mapping, geneData, type, exon.start(), exon.end());
        }

        List<RnaTarget> selected = targets.stream()
                .filter(target -> switch(target.type())
                {
                    case CODING -> true;
                    case UTR_5 -> options.utr5();
                    case UTR_3 -> options.utr3();
                })
                .sorted(Comparator.comparingInt(RnaTarget::spaceStart))
                .toList();
        return new GeneTargets(mapping, selected);
    }

    private static void addTarget(final List<RnaTarget> targets, final RegionMapping mapping, final GeneData geneData,
            final RnaRegionType type, int genomeStart, int genomeEnd)
    {
        int spaceStart = mapping.toProbeSpacePosition(geneData.Chromosome, genomeStart).orElseThrow();
        int spaceEnd = mapping.toProbeSpacePosition(geneData.Chromosome, genomeEnd).orElseThrow() + 1;
        targets.add(new RnaTarget(type, spaceStart, spaceEnd));
    }

    private static TargetMetadata createTargetMetadata(final GeneTranscriptData gene, final RnaTarget target)
    {
        // List the transcripts only when the user requested a specific subset; for the default (all transcripts) the extra info is just the
        // gene name and target type.
        String geneName = gene.gene().GeneName;
        String extraInfo;
        if(gene.specificTranscripts())
        {
            List<String> transcriptNames = gene.transcripts().stream().map(GenesRna::formatTranscriptName).toList();
            extraInfo = format("%s:%s:%s", geneName, join("/", transcriptNames), target.type().name());
        }
        else
        {
            extraInfo = format("%s:%s", geneName, target.type().name());
        }
        // Store the gene name so per-gene statistics can be recovered later.
        return new TargetMetadata(TARGET_TYPE, extraInfo, geneName);
    }

    private static String formatTranscriptName(final TranscriptData transcript)
    {
        return transcript.TransName;
    }

    private static List<GeneStats> computeGeneStats(final ProbeGenerationResult result, final List<GeneTranscriptData> genes)
    {
        Map<String, Integer> geneProbeCounts = new HashMap<>();
        for(GeneTranscriptData gene : genes)
        {
            geneProbeCounts.put(gene.gene().GeneName, 0);
        }
        for(Probe probe : result.probes())
        {
            String geneName = (String) requireNonNull(probe.metadata().extraData());
            geneProbeCounts.merge(geneName, 1, Integer::sum);
        }
        return geneProbeCounts.entrySet().stream()
                .sorted(Map.Entry.comparingByKey())
                .map(entry -> new GeneStats(entry.getKey(), entry.getValue()))
                .toList();
    }
}
