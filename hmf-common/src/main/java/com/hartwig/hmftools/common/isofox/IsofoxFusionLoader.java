package com.hartwig.hmftools.common.isofox;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.rna.ImmutableRnaFusion;
import com.hartwig.hmftools.common.rna.RnaCommon;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.jetbrains.annotations.NotNull;

public final class IsofoxFusionLoader
{
    private static final String FLD_CHR = "Chr";
    private static final String FLD_POS = "Pos";
    private static final String FLD_ORIENT = "Orient";
    private static final String FLD_JUNC_TYPE = "JuncType";
    private static final String FLD_SV_TYPE = "SVType";
    private static final String FLD_COVERAGE = "Coverage";
    private static final String FLD_SPLIT_FRAGS = "SplitFrags";
    private static final String FLD_REALIGN_FLAGS = "RealignedFrags";
    private static final String FLD_DISCORD_FRAGS = "DiscordantFrags";
    private static final String FLD_COHORT_COUNT = "CohortCount";

    @NotNull
    public static List<RnaFusion> load(@NotNull String isofoxFusionCsv) throws IOException
    {
        List<RnaFusion> fusions = Lists.newArrayList();

        List<String> lines = Files.readAllLines(Paths.get(isofoxFusionCsv));
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), RnaCommon.DELIMITER);

        for(String line : lines.subList(1, lines.size()))
        {
            String[] items = line.split(RnaCommon.DELIMITER, -1);

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

    @NotNull
    private static String formStreamField(@NotNull String field, int stream)
    {
        return field + (stream == FS_UP ? "Up" : "Down");
    }
}
