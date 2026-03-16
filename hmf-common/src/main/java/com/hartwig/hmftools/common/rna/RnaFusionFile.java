package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.RNA_LOGGER;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

public final class RnaFusionFile
{
    public static final String PASS_FUSION_FILE_ID = "pass_fusions.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ISF_FILE_ID + PASS_FUSION_FILE_ID;
    }

    private enum Columns
    {
        Name,
        KnownType,
        ChromosomeUp,
        ChromosomeDown,
        PositionUp,
        PositionDown,
        OrientationUp,
        OrientationDown,
        JunctionTypeUp,
        JunctionTypeDown,
        TranscriptUp,
        TranscriptDown,
        ExonUp,
        ExonDown,
        SvType,
        SplitFrags,
        RealignedFrags,
        DiscordantFrags,
        DepthUp,
        DepthDown,
        MaxAnchorLengthUp,
        MaxAnchorLengthDown,
        CohortFrequency;
    }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        Arrays.stream(Columns.values()).forEach(x -> sj.add(x.toString()));
        return sj.toString();
    }

    public static final String FUSION_GENE_DELIM = "_";

    public static String formFusionName(final String geneUp, final String geneDown)
    {
        return geneUp + FUSION_GENE_DELIM + geneDown;
    }

    public static String[] geneNames(final RnaFusion fusion) { return fusion.name().split(FUSION_GENE_DELIM, 2); }

    public static String geneUp(final RnaFusion fusion)
    {
        String[] genes = fusion.name().split(FUSION_GENE_DELIM, 2);
        return genes.length == 2 ? genes[0] : null;
    }

    public static String geneDown(final RnaFusion fusion)
    {
        String[] genes = fusion.name().split(FUSION_GENE_DELIM, 2);
        return genes.length == 2 ? genes[1] : null;
    }

    public static String write(final RnaFusion fusion)
    {
        return new StringJoiner(TSV_DELIM)
                .add(fusion.name())
                .add(String.valueOf(fusion.knownType()))
                .add(fusion.chromosomeUp())
                .add(fusion.chromosomeDown())
                .add(String.valueOf(fusion.positionUp()))
                .add(String.valueOf(fusion.positionDown()))
                .add(String.valueOf(fusion.orientationUp()))
                .add(String.valueOf(fusion.orientationDown()))
                .add(fusion.junctionTypeUp())
                .add(fusion.junctionTypeDown())
                .add(fusion.transcriptUp())
                .add(fusion.transcriptDown())
                .add(String.valueOf(fusion.exonUp()))
                .add(String.valueOf(fusion.exonDown()))
                .add(String.valueOf(fusion.svType()))
                .add(String.valueOf(fusion.splitFragments()))
                .add(String.valueOf(fusion.realignedFrags()))
                .add(String.valueOf(fusion.discordantFrags()))
                .add(String.valueOf(fusion.depthUp()))
                .add(String.valueOf(fusion.depthDown()))
                .add(String.valueOf(fusion.maxAnchorLengthUp()))
                .add(String.valueOf(fusion.maxAnchorLengthDown()))
                .add(String.valueOf(fusion.cohortFrequency()))
                .toString();
    }

    public static List<RnaFusion> read(final String filename)
    {
        try
        {
            List<RnaFusion> fusions = Lists.newArrayList();

            List<String> lines = Files.readAllLines(Paths.get(filename));
            String fileDelim = inferFileDelimiter(filename);
            Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), fileDelim);

            // pre v2.0 columns
            Integer geneUpIndex = fieldsIndexMap.get("GeneNameUp");
            Integer geneDownIndex = fieldsIndexMap.get("GeneNameDown");

            int chrUpIndex = fieldsIndexMap.getOrDefault(Columns.ChromosomeUp.toString(), fieldsIndexMap.get("ChrUp"));
            int chrDownIndex = fieldsIndexMap.getOrDefault(Columns.ChromosomeDown.toString(), fieldsIndexMap.get("ChrDown"));
            int posUpIndex = fieldsIndexMap.getOrDefault(Columns.PositionUp.toString(), fieldsIndexMap.get("PosUp"));
            int posDownIndex = fieldsIndexMap.getOrDefault(Columns.PositionDown.toString(), fieldsIndexMap.get("PosDown"));
            int orientUpIndex = fieldsIndexMap.getOrDefault(Columns.OrientationUp.toString(), fieldsIndexMap.get("OrientUp"));
            int orientDownIndex = fieldsIndexMap.getOrDefault(Columns.OrientationDown.toString(), fieldsIndexMap.get("OrientDown"));
            int juncTypeUpIndex = fieldsIndexMap.getOrDefault(Columns.JunctionTypeUp.toString(), fieldsIndexMap.get("JuncTypeUp"));
            int juncTypeDownIndex = fieldsIndexMap.getOrDefault(Columns.JunctionTypeDown.toString(), fieldsIndexMap.get("JuncTypeDown"));
            int svTypeIndex = fieldsIndexMap.getOrDefault(Columns.SvType.toString(), fieldsIndexMap.get("SVType"));
            int depthUpIndex = fieldsIndexMap.getOrDefault(Columns.DepthUp.toString(), fieldsIndexMap.get("CoverageUp"));
            int depthDownIndex = fieldsIndexMap.getOrDefault(Columns.DepthDown.toString(), fieldsIndexMap.get("CoverageDown"));
            int cohortIndex = fieldsIndexMap.getOrDefault(Columns.CohortFrequency.toString(), fieldsIndexMap.get("CohortCount"));

            // v2.0 columns
            Integer transcriptUpIndex = fieldsIndexMap.get(Columns.TranscriptUp);
            Integer transcriptDownIndex = fieldsIndexMap.get(Columns.TranscriptDown);
            Integer exonUpIndex = fieldsIndexMap.get(Columns.ExonUp);
            Integer exonDownIndex = fieldsIndexMap.get(Columns.ExonDown);
            Integer knownTypeIndex = fieldsIndexMap.get(Columns.KnownType);

            for(String line : lines.subList(1, lines.size()))
            {
                String[] values = line.split(fileDelim, -1);

                String transcriptUp = transcriptUpIndex != null ? values[transcriptUpIndex] : "";
                String transcriptDown = transcriptDownIndex != null ? values[transcriptDownIndex] : "";
                int exonUp = exonUpIndex != null ? Integer.parseInt(values[exonUpIndex]) : -1;
                int exonDown = exonDownIndex != null ? Integer.parseInt(values[exonDownIndex]) : -1;

                String fusionName;

                if(geneUpIndex != null && geneDownIndex != null)
                {
                    fusionName = formFusionName(values[geneUpIndex], values[geneDownIndex]);
                }
                else
                {
                    fusionName = values[fieldsIndexMap.get(Columns.Name.toString())];
                }

                KnownFusionType knownFusionType = knownTypeIndex != null ? parseKnownType(values[knownTypeIndex]) : KnownFusionType.NONE;

                fusions.add(ImmutableRnaFusion.builder()
                        .name(fusionName)
                        .knownType(knownFusionType)
                        .chromosomeUp(values[chrUpIndex])
                        .chromosomeDown(values[chrDownIndex])
                        .positionUp(Integer.parseInt(values[posUpIndex]))
                        .positionDown(Integer.parseInt(values[posDownIndex]))
                        .orientationUp(Byte.parseByte(values[orientUpIndex]))
                        .orientationDown(Byte.parseByte(values[orientDownIndex]))
                        .junctionTypeUp(values[juncTypeUpIndex])
                        .junctionTypeDown(values[juncTypeDownIndex])
                        .transcriptUp(transcriptUp)
                        .transcriptDown(transcriptDown)
                        .exonUp(exonUp)
                        .exonDown(exonDown)
                        .svType(StructuralVariantType.valueOf(values[svTypeIndex]))
                        .splitFragments(Integer.parseInt(values[fieldsIndexMap.get(Columns.SplitFrags.toString())]))
                        .realignedFrags(Integer.parseInt(values[fieldsIndexMap.get(Columns.RealignedFrags.toString())]))
                        .discordantFrags(Integer.parseInt(values[fieldsIndexMap.get(Columns.DiscordantFrags.toString())]))
                        .depthUp(Integer.parseInt(values[depthUpIndex]))
                        .depthDown(Integer.parseInt(values[depthDownIndex]))
                        .maxAnchorLengthUp(Integer.parseInt(values[fieldsIndexMap.get(Columns.MaxAnchorLengthUp.toString())]))
                        .maxAnchorLengthDown(Integer.parseInt(values[fieldsIndexMap.get(Columns.MaxAnchorLengthDown.toString())]))
                        .cohortFrequency(Integer.parseInt(values[cohortIndex]))
                        .build());
            }

            return fusions;
        }
        catch(IOException e)
        {
            RNA_LOGGER.error("failed to load Isofox fusion file({}): {}", filename, e.toString());
            return null;
        }
    }

    private static KnownFusionType parseKnownType(final String knownTypeStr)
    {
        try
        {
            return KnownFusionType.valueOf(knownTypeStr);
        }
        catch(Exception e)
        {
            if(knownTypeStr.contains("PROM5"))
                return KnownFusionType.PROMISCUOUS_5;

            if(knownTypeStr.contains("PROM3"))
                return KnownFusionType.PROMISCUOUS_3;

            return KnownFusionType.NONE;
        }
    }
}
