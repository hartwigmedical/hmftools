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
        SplitFragments,
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

            for(String line : lines.subList(1, lines.size()))
            {
                String[] values = line.split(fileDelim, -1);

                fusions.add(ImmutableRnaFusion.builder()
                        .name(values[fieldsIndexMap.get(Columns.Name.toString())])
                        .chromosomeUp(values[fieldsIndexMap.get(Columns.ChromosomeUp.toString())])
                        .chromosomeDown(values[fieldsIndexMap.get(Columns.ChromosomeDown.toString())])
                        .positionUp(Integer.parseInt(values[fieldsIndexMap.get(Columns.PositionUp)]))
                        .positionDown(Integer.parseInt(values[fieldsIndexMap.get(Columns.PositionDown)]))
                        .orientationUp(Byte.parseByte(values[fieldsIndexMap.get(Columns.OrientationUp)]))
                        .orientationDown(Byte.parseByte(values[fieldsIndexMap.get(Columns.OrientationDown)]))
                        .junctionTypeUp(values[fieldsIndexMap.get(Columns.JunctionTypeUp)])
                        .junctionTypeDown(values[fieldsIndexMap.get(Columns.JunctionTypeDown)])
                        .knownType(KnownFusionType.valueOf(values[fieldsIndexMap.get(Columns.KnownType)]))
                        .svType(StructuralVariantType.valueOf(values[fieldsIndexMap.get(Columns.SvType)]))
                        .splitFragments(Integer.parseInt(values[fieldsIndexMap.get(Columns.SplitFragments)]))
                        .realignedFrags(Integer.parseInt(values[fieldsIndexMap.get(Columns.RealignedFrags)]))
                        .discordantFrags(Integer.parseInt(values[fieldsIndexMap.get(Columns.DiscordantFrags)]))
                        .depthUp(Integer.parseInt(values[fieldsIndexMap.get(Columns.DepthUp)]))
                        .depthDown(Integer.parseInt(values[fieldsIndexMap.get(Columns.DepthDown)]))
                        .maxAnchorLengthUp(Integer.parseInt(values[fieldsIndexMap.get(Columns.MaxAnchorLengthUp)]))
                        .maxAnchorLengthDown(Integer.parseInt(values[fieldsIndexMap.get(Columns.MaxAnchorLengthDown)]))
                        .cohortFrequency(Integer.parseInt(values[fieldsIndexMap.get(Columns.CohortFrequency)]))
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
}
