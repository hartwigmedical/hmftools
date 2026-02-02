package com.hartwig.hmftools.common.hla;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.checkFileExtensionRename;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class LilacQcData
{
    public abstract String genes();
    public abstract String status();
    public abstract int totalFragments();
    public abstract int fittedFragments();
    public abstract int discardedIndels();
    public abstract int discardedAlignmentFragments();
    public abstract String hlaYAllele();

    private static final String FILE_EXTENSION = ".lilac.qc.tsv";

    public static final String FLD_GENES = "Genes";
    public static final String FLD_QC_STATUS = "Status";
    public static final String FLD_HLA_Y = "HlaYAllele";
    public static final String FLD_TOTAL_FRAGS = "TotalFragments";
    public static final String FLD_FIT_FRAGS = "FittedFragments";
    public static final String FLD_DISC_INDELS = "DiscardedIndels";
    public static final String FLD_DISC_ALIGN_FRAGS = "DiscardedAlignmentFragments";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + FILE_EXTENSION;
    }

    public static List<LilacQcData> read(final String filePath) throws IOException
    {
        String filename = checkFileExtensionRename(filePath);
        String delim = inferFileDelimiter(filename);

        List<String> lines = Files.readAllLines(Paths.get(filename));

        final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), delim);

        List<LilacQcData> qcData = Lists.newArrayList();
        for(int i = 1; i < lines.size(); i++)
        {
            String[] values = lines.get(i).split(delim);
            ImmutableLilacQcData data = ImmutableLilacQcData.builder()
                    .genes(values[fieldsIndexMap.get(FLD_GENES)])
                    .status(values[fieldsIndexMap.get(FLD_QC_STATUS)])
                    .totalFragments(Integer.parseInt(values[fieldsIndexMap.get(FLD_TOTAL_FRAGS)]))
                    .fittedFragments(Integer.parseInt(values[fieldsIndexMap.get(FLD_FIT_FRAGS)]))
                    .discardedAlignmentFragments(Integer.parseInt(values[fieldsIndexMap.get(FLD_DISC_ALIGN_FRAGS)]))
                    .discardedIndels(Integer.parseInt(values[fieldsIndexMap.get(FLD_DISC_INDELS)]))
                    .hlaYAllele(values[fieldsIndexMap.get(FLD_HLA_Y)])
                    .build();

            qcData.add(data);
        }

        return qcData;
    }
}
