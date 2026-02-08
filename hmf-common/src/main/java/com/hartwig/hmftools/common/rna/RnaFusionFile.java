package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.RNA_LOGGER;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

public final class RnaFusionFile
{
    public static final String UNFILTERED_FUSION_FILE_ID = "fusions.tsv";
    public static final String PASS_FUSION_FILE_ID = "pass_fusions.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + ISF_FILE_ID + PASS_FUSION_FILE_ID;
    }

    public static final String FLD_CHR = "Chr";
    public static final String FLD_POS = "Pos";
    public static final String FLD_ORIENT = "Orient";
    public static final String FLD_JUNC_TYPE = "JuncType";
    public static final String FLD_SV_TYPE = "SVType";
    public static final String FLD_COVERAGE = "Coverage";
    public static final String FLD_KNOWN_TYPE = "KnownFusionType";
    public static final String FLD_SPLIT_FRAGS = "SplitFrags";
    public static final String FLD_REALIGN_FLAGS = "RealignedFrags";
    public static final String FLD_DISCORD_FRAGS = "DiscordantFrags";
    public static final String FLD_COHORT_COUNT = "CohortCount";

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
                String[] items = line.split(fileDelim, -1);

                String fusionName = String.format("%s_%s",
                        items[fieldsIndexMap.get(formStreamField(FLD_GENE_NAME, FS_UP))],
                        items[fieldsIndexMap.get(formStreamField(FLD_GENE_NAME, FS_DOWN))]);

                fusions.add(ImmutableRnaFusion.builder()
                        .name(fusionName)
                        .chromosomeUp(items[fieldsIndexMap.get(formStreamField(FLD_CHR, FS_UP))])
                        .chromosomeDown(items[fieldsIndexMap.get(formStreamField(FLD_CHR, FS_DOWN))])
                        .positionUp(Integer.parseInt(items[fieldsIndexMap.get(formStreamField(FLD_POS, FS_UP))]))
                        .positionDown(Integer.parseInt(items[fieldsIndexMap.get(formStreamField(FLD_POS, FS_DOWN))]))
                        .orientationUp(Byte.parseByte(items[fieldsIndexMap.get(formStreamField(FLD_ORIENT, FS_UP))]))
                        .orientationDown(Byte.parseByte(items[fieldsIndexMap.get(formStreamField(FLD_ORIENT, FS_DOWN))]))
                        .junctionTypeUp(items[fieldsIndexMap.get(formStreamField(FLD_JUNC_TYPE, FS_UP))])
                        .junctionTypeDown(items[fieldsIndexMap.get(formStreamField(FLD_JUNC_TYPE, FS_DOWN))])
                        .knownType(KnownFusionType.valueOf(items[fieldsIndexMap.get(FLD_KNOWN_TYPE)]))
                        .svType(StructuralVariantType.valueOf(items[fieldsIndexMap.get(FLD_SV_TYPE)]))
                        .splitFragments(Integer.parseInt(items[fieldsIndexMap.get(FLD_SPLIT_FRAGS)]))
                        .realignedFrags(Integer.parseInt(items[fieldsIndexMap.get(FLD_REALIGN_FLAGS)]))
                        .discordantFrags(Integer.parseInt(items[fieldsIndexMap.get(FLD_DISCORD_FRAGS)]))
                        .depthUp(Integer.parseInt(items[fieldsIndexMap.get(formStreamField(FLD_COVERAGE, FS_UP))]))
                        .depthDown(Integer.parseInt(items[fieldsIndexMap.get(formStreamField(FLD_COVERAGE, FS_DOWN))]))
                        .maxAnchorLengthUp(Integer.parseInt(items[fieldsIndexMap.get(formStreamField(FLD_COVERAGE, FS_UP))]))
                        .maxAnchorLengthDown(Integer.parseInt(items[fieldsIndexMap.get(formStreamField(FLD_COVERAGE, FS_DOWN))]))
                        .cohortFrequency(Integer.parseInt(items[fieldsIndexMap.get(FLD_COHORT_COUNT)]))
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

    private static String formStreamField(final String field, int stream)
    {
        return field + (stream == FS_UP ? "Up" : "Down");
    }
}
