package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;
import static java.lang.String.join;
import static java.util.Collections.emptyList;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.panelbuilder.GeneUtils.mergeExons;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_RNA_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_RNA_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_RNA_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
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
// a per-gene RegionMapping (the merged exons as a contiguous probe-space) and ProbeGenerator.coverExonRange.
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
            final ProbeGenerator probeGenerator, PanelData panelData)
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

    private static ProbeGenerationResult generateProbes(final List<GeneTranscriptData> genes, final ProbeGenerator probeGenerator,
            PanelData panelData)
    {
        // Build one CoverExonRange spec per exon target and submit them as a single batch, so the generator can batch the (expensive)
        // junction-crossing quality-score alignments across all exons rather than per exon.
        List<ProbeGenerationSpec> specs = new ArrayList<>();
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
                specs.add(new ProbeGenerationSpec.CoverExonRange(
                        geneTargets.mapping(), target.spaceStart(), target.spaceEnd(), metadata, PROBE_CRITERIA, PROBE_SELECT));
            }
        }
        return probeGenerator.generateBatch(specs.stream(), panelData);
    }

    // Computes the exon-aware target ranges for a gene: the coding sequence (always) and, if requested, the 5' and 3' UTRs. All exons of the
    // selected transcripts are merged into one RegionMapping (so short-exon probes can pad across junctions), then each merged exon is
    // classified by how many of its bases are coding vs noncoding:
    //   - coding exon (some coding, less than a probe of noncoding): one whole-exon coding target;
    //   - noncoding exon (some noncoding, less than a probe of coding): one whole-exon UTR target (5' or 3' by position, strand-aware);
    //   - partially coding exon (at least a probe of each): split into a coding target plus the flanking UTR target(s), so a long exon's UTR
    //     is not covered as coding.
    // A whole exon or exon part is therefore probed once, with a single feature type. UTR targets are kept only if the corresponding UTR is
    // enabled.
    static GeneTargets createTargets(final GeneData geneData, final List<TranscriptData> transcripts, final GeneOptions options)
    {
        List<GeneUtils.MergedExonRegion> mergedExons = mergeExons(transcripts);

        List<ChrBaseRegion> exonRegions = mergedExons.stream()
                .map(merged -> ChrBaseRegion.from(geneData.Chromosome, merged.Region))
                .toList();
        RegionMapping mapping = new RegionMapping(exonRegions);

        // Gene-wide coding start (across coding transcripts); used to classify a fully noncoding exon as 5' or 3'.
        OptionalInt geneCodingStart = transcripts.stream().filter(t -> !t.nonCoding()).mapToInt(t -> t.CodingStart).min();
        if(geneCodingStart.isEmpty())
        {
            // Noncoding gene: no coding span, so there is no 5'/3' distinction. Cover each exon as UTR if either UTR is requested.
            return noncodingGeneTargets(geneData, mergedExons, mapping, options);
        }
        int codingStart = geneCodingStart.getAsInt();
        boolean forwardStrand = geneData.forwardStrand();

        List<RnaTarget> targets = new ArrayList<>();
        for(GeneUtils.MergedExonRegion mergedExon : mergedExons)
        {
            BaseRegion exon = mergedExon.Region;
            // Coding bases are the exon's overlap with the coding span; the rest are noncoding (UTR).
            int codingStartInExon = mergedExon.IsCoding ? max(exon.start(), mergedExon.CodingStart) : 0;
            int codingEndInExon = mergedExon.IsCoding ? min(exon.end(), mergedExon.CodingEnd) : -1;
            int codingBases = mergedExon.IsCoding ? max(0, codingEndInExon - codingStartInExon + 1) : 0;
            int nonCodingBases = exon.baseLength() - codingBases;

            boolean isCodingExon = codingBases > 0 && nonCodingBases < PROBE_LENGTH;
            boolean isNonCodingExon = nonCodingBases > 0 && codingBases < PROBE_LENGTH;

            // TODO: small coding part of a boundary exon. When the coding part is shorter than a probe it is folded into a whole-exon target
            //  (coding here, or UTR if the coding part is the smaller side). To give it a coding-specific probe, split it out and pad the short
            //  coding probe into the adjacent same-exon UTR (contiguous, single-region) rather than across the splice junction into the next
            //  exon (spliced). See RNA_DESIGN_NOTES "Planned - small coding part padding" (follow-up #4).
            if(isCodingExon)
            {
                // Any noncoding portion is shorter than a probe, so cover the whole exon as coding.
                addTarget(targets, mapping, geneData, RnaRegionType.CODING, exon.start(), exon.end());
            }
            else if(isNonCodingExon)
            {
                // Any coding portion is shorter than a probe, so cover the whole exon as one UTR feature.
                addTarget(targets, mapping, geneData, utrType(exon.end() < codingStart, forwardStrand), exon.start(), exon.end());
            }
            else
            {
                // Partially coding: coding and noncoding parts are each at least a probe long, so cover them separately.
                addTarget(targets, mapping, geneData, RnaRegionType.CODING, codingStartInExon, codingEndInExon);
                if(exon.start() < codingStartInExon)
                {
                    addTarget(targets, mapping, geneData, utrType(true, forwardStrand), exon.start(), codingStartInExon - 1);
                }
                if(exon.end() > codingEndInExon)
                {
                    addTarget(targets, mapping, geneData, utrType(false, forwardStrand), codingEndInExon + 1, exon.end());
                }
            }
        }

        // Each exon (or exon part) produced exactly one target with one feature type, so enabling both UTRs cannot probe a region twice.
        // Keep only the enabled UTR features (coding is always kept).
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

    // Targets for a noncoding gene (no coding transcripts): the transcript is all noncoding exonic sequence, so it is covered as UTR when
    // either UTR is requested. With no coding span there is no real 5' vs 3' distinction, so each exon is one whole-exon UTR target labelled
    // with whichever UTR is enabled (5' if both). If neither UTR is requested, nothing is produced.
    private static GeneTargets noncodingGeneTargets(final GeneData geneData, final List<GeneUtils.MergedExonRegion> mergedExons,
            final RegionMapping mapping, final GeneOptions options)
    {
        if(!options.utr5() && !options.utr3())
        {
            LOGGER.debug("Gene {} is noncoding and no UTR requested; no RNA probes produced", geneData.GeneName);
            return new GeneTargets(mapping, emptyList());
        }
        RnaRegionType utrType = options.utr5() ? RnaRegionType.UTR_5 : RnaRegionType.UTR_3;
        List<RnaTarget> targets = new ArrayList<>();
        for(GeneUtils.MergedExonRegion mergedExon : mergedExons)
        {
            addTarget(targets, mapping, geneData, utrType, mergedExon.Region.start(), mergedExon.Region.end());
        }
        targets.sort(Comparator.comparingInt(RnaTarget::spaceStart));
        return new GeneTargets(mapping, targets);
    }

    // 5' UTR below the coding span on the forward strand (3' on the reverse strand), and vice versa above it.
    private static RnaRegionType utrType(boolean belowCoding, boolean forwardStrand)
    {
        return belowCoding == forwardStrand ? RnaRegionType.UTR_5 : RnaRegionType.UTR_3;
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
