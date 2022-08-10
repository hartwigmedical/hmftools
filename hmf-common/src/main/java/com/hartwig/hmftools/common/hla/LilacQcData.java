package com.hartwig.hmftools.common.hla;

import static com.hartwig.hmftools.common.hla.LilacAllele.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class LilacQcData
{
    // full field list
    // Status,ScoreMargin,NextSolutionAlleles,MedianBaseQuality,HlaYAllele,DiscardedIndels,DiscardedIndelMaxFrags,DiscardedAlignmentFragments,
    // A_LowCoverageBases,B_LowCoverageBases,C_LowCoverageBases,ATypes,BTypes,CTypes,TotalFragments,FittedFragments,UnmatchedFragments,
    // UninformativeFragments,HlaYFragments,PercentUnique,PercentShared,PercentWildcard,UnusedAminoAcids,UnusedAminoAcidMaxFrags,UnusedHaplotypes,
    // UnusedHaplotypeMaxFrags,SomaticVariantsMatched,SomaticVariantsUnmatched

    public abstract String status();

    private static final String FILE_EXTENSION = ".lilac.qc.csv";
    public static final String FLD_QC_STATUS = "Status";

    public static String generateFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static LilacQcData read(final String filePath) throws IOException
    {
        List<String> lines = Files.readAllLines(Paths.get(filePath));

        final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);

        String[] values = lines.get(1).split(DELIMITER);

        return ImmutableLilacQcData.builder()
                    .status(values[fieldsIndexMap.get(FLD_QC_STATUS)])
                    .build();
    }
}
