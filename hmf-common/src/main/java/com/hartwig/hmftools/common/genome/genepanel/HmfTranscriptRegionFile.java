package com.hartwig.hmftools.common.genome.genepanel;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.ImmutableHmfExonRegion;
import com.hartwig.hmftools.common.genome.region.ImmutableHmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class HmfTranscriptRegionFile
{
    private static final Logger LOGGER = LogManager.getLogger(HmfTranscriptRegionFile.class);

    public static final String DEFAULT_DELIM = "\t";
    public static final String ITEM_DELIM = ";";
    public static final String EXON_DATA_DELIM = ":";

    public static String header() { return header(DEFAULT_DELIM); }

    public static String header(final String delim)
    {
        StringJoiner header = new StringJoiner(delim);
        header.add("GeneId");
        header.add("GeneName");
        header.add("Chromosome");
        header.add("GeneStart");
        header.add("GeneEnd");
        header.add("Strand");
        header.add("EntrezIds");
        header.add("ChrBand");
        header.add("TranscriptId");
        header.add("IsCanonical");
        header.add("TransStart");
        header.add("TransEnd");
        header.add("CodingStart");
        header.add("CodingEnd");
        header.add("ExonData");
        return header.toString();
    }

    public String toString(final HmfTranscriptRegion transRegion) { return toString(transRegion, DEFAULT_DELIM); }

    public String toString(final HmfTranscriptRegion transRegion, final String delim)
    {
        StringJoiner sj = new StringJoiner(delim);
        sj.add(transRegion.geneId());
        sj.add(transRegion.geneName());
        sj.add(transRegion.chromosome());
        sj.add(String.valueOf(transRegion.geneStart()));
        sj.add(String.valueOf(transRegion.geneEnd()));
        sj.add(String.valueOf(transRegion.strand()));

        StringJoiner entrezIds = new StringJoiner(ITEM_DELIM);
        transRegion.entrezId().forEach(x -> entrezIds.add(String.valueOf(x)));
        sj.add(entrezIds.toString());
        sj.add(transRegion.chromosomeBand());

        sj.add(transRegion.transName());
        sj.add(String.valueOf(transRegion.isCanonical()));
        sj.add(String.valueOf(transRegion.start()));
        sj.add(String.valueOf(transRegion.end()));

        sj.add(String.valueOf(transRegion.codingStart()));
        sj.add(String.valueOf(transRegion.codingEnd()));

        StringJoiner exonData = new StringJoiner(ITEM_DELIM);
        transRegion.exons().forEach(x -> exonData.add(String.format("%d%s%s", x.start(), EXON_DATA_DELIM, x.end())));
        sj.add(exonData.toString());

        return sj.toString();
    }

    @NotNull
    public static List<HmfTranscriptRegion> fromInputStream(@NotNull InputStream genomeInputStream)
    {
        return fromLines(new BufferedReader(new InputStreamReader(genomeInputStream)).lines().collect(Collectors.toList()));
    }

    public static List<HmfTranscriptRegion> fromLines(final List<String> lines)
    {
        List<HmfTranscriptRegion> transcriptRegions = Lists.newArrayList();

        Set<String> uniqueGenes = Sets.newHashSet();

        if (lines.isEmpty())
        {
            LOGGER.error("empty transcript region data file");
            return transcriptRegions;
        }

        final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DEFAULT_DELIM);
        lines.remove(0);

        int geneIdIndex = fieldsIndexMap.get("GeneId");
        int geneNameIndex = fieldsIndexMap.get("GeneName");
        int chromosomeIndex = fieldsIndexMap.get("Chromosome");
        int strandIndex = fieldsIndexMap.get("Strand");
        int geneStartIndex = fieldsIndexMap.get("GeneStart");
        int geneEndIndex = fieldsIndexMap.get("GeneEnd");
        int entrezIdsIndex = fieldsIndexMap.get("EntrezIds");
        int chrBandIndex = fieldsIndexMap.get("ChrBand");
        int transIdIndex = fieldsIndexMap.get("TranscriptId");
        Integer canonicalIndex = fieldsIndexMap.get("IsCanonical");
        int transStartIndex = fieldsIndexMap.get("TransStart");
        int transEndIndex = fieldsIndexMap.get("TransEnd");
        int codingStartIndex = fieldsIndexMap.get("CodingStart");
        int codingEndIndex = fieldsIndexMap.get("CodingEnd");
        int exonDataIndex = fieldsIndexMap.get("ExonData");

        for(final String line : lines)
        {
            String[] values = line.split(DEFAULT_DELIM);

            String chromosome = values[chromosomeIndex].trim();
            if(!HumanChromosome.contains(chromosome))
            {
                LOGGER.warn("invalid chromosome from line: {}", line);
                continue;
            }

            String geneId = values[geneIdIndex];
            String geneName = values[geneNameIndex];

            if(uniqueGenes.contains(geneName))
                continue;

            uniqueGenes.add(geneName);

            String entrezIdsStr = values[entrezIdsIndex];

            List<Integer> entrezIds = entrezIdsStr.isEmpty()
                    ? Lists.newArrayList()
                    : Arrays.stream(entrezIdsStr.split(ITEM_DELIM, -1)).map(Integer::parseInt).collect(Collectors.toList());

            String[] exonsData = values[exonDataIndex].split(ITEM_DELIM, -1);

            List<HmfExonRegion> exons = Lists.newArrayList();

            for(String exonData : exonsData)
            {
                String[] exonItems = exonData.split(EXON_DATA_DELIM);
                exons.add(ImmutableHmfExonRegion.builder()
                        .chromosome(chromosome)
                        .exonRank(0) // unsed and will be deprecated
                        .start(Integer.parseInt(exonItems[0]))
                        .end(Integer.parseInt(exonItems[1]))
                        .build());
            }

            long codingStart = values[codingStartIndex].isEmpty() || values[codingStartIndex].equals("NA")
                    ? 0 : Integer.parseInt(values[codingStartIndex]);

            long codingEnd = values[codingEndIndex].isEmpty() || values[codingEndIndex].equals("NA")
                    ? 0 : Integer.parseInt(values[codingEndIndex]);

            transcriptRegions.add(ImmutableHmfTranscriptRegion.builder()
                    .geneId(geneId)
                    .geneName(geneName)
                    .chromosome(chromosome)
                    .strand(Strand.valueOf(Integer.parseInt(values[strandIndex])))
                    .geneStart(Integer.parseInt(values[geneStartIndex]))
                    .geneEnd(Integer.parseInt(values[geneEndIndex]))
                    .entrezId(entrezIds)
                    .chromosomeBand(values[chrBandIndex])
                    .transName(values[transIdIndex])
                    .isCanonical(canonicalIndex != null ? Boolean.parseBoolean(values[canonicalIndex]) : true)
                    .start(Integer.parseInt(values[transStartIndex]))
                    .end(Integer.parseInt(values[transEndIndex]))
                    .codingStart(codingStart)
                    .codingEnd(codingEnd)
                    .exons(exons)
                    .build());
        }

        return transcriptRegions;
    }
}