package com.hartwig.hmftools.common.hla;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class LilacAllele {

    public abstract String allele();
    public abstract int refFragments();
    public abstract int refUnique();
    public abstract int refShared();
    public abstract int refWild();
    public abstract int tumorFragments();
    public abstract int tumorUnique();
    public abstract int tumorShared();
    public abstract int tumorWild();
    public abstract int rnaFragments();
    public abstract int rnaUnique();
    public abstract int rnaShared();
    public abstract int rnaWild();

    public abstract double tumorCopyNumber();

    public abstract double somaticMissense();
    public abstract double somaticNonsenseOrFrameshift();
    public abstract double somaticSplice();
    public abstract double somaticSynonymous();
    public abstract double somaticInframeIndel();

    public double somaticVariantCount()
    {
        return somaticMissense() + somaticNonsenseOrFrameshift() + somaticSplice() + somaticSynonymous() + somaticInframeIndel();
    }

    protected static final String DELIMITER = ",";
    private static final String FILE_EXTENSION = ".lilac.csv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static List<LilacAllele> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(final String filename, List<LilacAllele> alleles) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(alleles));
    }

    static List<String> toLines(final List<LilacAllele> alleles)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        alleles.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    private static final String FLD_ALLELE = "Allele";
    private static final String FLD_REF_TOTAL = "RefTotal";
    private static final String FLD_REF_UNIQUE = "RefUnique";
    private static final String FLD_REF_SHARED = "RefShared";
    private static final String FLD_REF_WILD = "RefWild";
    private static final String FLD_TUMOR_TOTAL = "TumorTotal";
    private static final String FLD_TUMOR_UNIQUE = "TumorUnique";
    private static final String FLD_TUMOR_SHARED = "TumorShared";
    private static final String FLD_TUMOR_WILD = "TumorWild";
    private static final String FLD_TUMOR_CN = "TumorCopyNumber";
    private static final String FLD_RNA_TOTAL = "RnaTotal";
    private static final String FLD_RNA_UNIQUE = "RnaUnique";
    private static final String FLD_RNA_SHARED = "RnaShared";
    private static final String FLD_RNA_WILD = "RnaWild";
    private static final String FLD_MISSENSE = "SomaticMissense";
    private static final String FLD_NFS = "SomaticNonsenseOrFrameshift";
    private static final String FLD_SPLICE = "SomaticSplice";
    private static final String FLD_SYNON = "SomaticSynonymous";
    private static final String FLD_INDEL = "SomaticInframeIndel";

    static List<LilacAllele> fromLines(final List<String> lines)
    {
        final String header = lines.get(0);
        lines.remove(0);

        final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header,DELIMITER);

        List<LilacAllele> alleles = Lists.newArrayList();

        for(int i = 0; i < lines.size(); ++i)
        {
            String[] values = lines.get(i).split(DELIMITER);

            alleles.add(ImmutableLilacAllele.builder()
                    .allele(values[fieldsIndexMap.get(FLD_ALLELE)])
                    .refFragments(Integer.parseInt(values[fieldsIndexMap.get(FLD_REF_TOTAL)]))
                    .refUnique(Integer.parseInt(values[fieldsIndexMap.get(FLD_REF_UNIQUE)]))
                    .refShared(Integer.parseInt(values[fieldsIndexMap.get(FLD_REF_SHARED)]))
                    .refWild(Integer.parseInt(values[fieldsIndexMap.get(FLD_REF_WILD)]))
                    .tumorFragments(Integer.parseInt(values[fieldsIndexMap.get(FLD_TUMOR_TOTAL)]))
                    .tumorUnique(Integer.parseInt(values[fieldsIndexMap.get(FLD_TUMOR_UNIQUE)]))
                    .tumorShared(Integer.parseInt(values[fieldsIndexMap.get(FLD_TUMOR_SHARED)]))
                    .tumorWild(Integer.parseInt(values[fieldsIndexMap.get(FLD_TUMOR_WILD)]))
                    .rnaFragments(Integer.parseInt(values[fieldsIndexMap.get(FLD_RNA_TOTAL)]))
                    .rnaUnique(Integer.parseInt(values[fieldsIndexMap.get(FLD_RNA_UNIQUE)]))
                    .rnaShared(Integer.parseInt(values[fieldsIndexMap.get(FLD_RNA_SHARED)]))
                    .rnaWild(Integer.parseInt(values[fieldsIndexMap.get(FLD_RNA_WILD)]))
                    .tumorCopyNumber(Double.parseDouble(values[fieldsIndexMap.get(FLD_TUMOR_CN)]))
                    .somaticMissense(Double.parseDouble(values[fieldsIndexMap.get(FLD_MISSENSE)]))
                    .somaticNonsenseOrFrameshift(Double.parseDouble(values[fieldsIndexMap.get(FLD_NFS)]))
                    .somaticSplice(Double.parseDouble(values[fieldsIndexMap.get(FLD_SPLICE)]))
                    .somaticSynonymous(Double.parseDouble(values[fieldsIndexMap.get(FLD_SYNON)]))
                    .somaticInframeIndel(Double.parseDouble(values[fieldsIndexMap.get(FLD_INDEL)]))
                    .build());
        }

        return alleles;
    }

    public static String header()
    {
        return new StringJoiner(DELIMITER)
                .add(FLD_ALLELE)
                .add(FLD_REF_TOTAL)
                .add(FLD_REF_UNIQUE)
                .add(FLD_REF_SHARED)
                .add(FLD_REF_WILD)
                .add(FLD_TUMOR_TOTAL)
                .add(FLD_TUMOR_UNIQUE)
                .add(FLD_TUMOR_SHARED)
                .add(FLD_TUMOR_WILD)
                .add(FLD_RNA_TOTAL)
                .add(FLD_RNA_UNIQUE)
                .add(FLD_RNA_SHARED)
                .add(FLD_RNA_WILD)
                .add(FLD_TUMOR_CN)
                .add(FLD_MISSENSE)
                .add(FLD_NFS)
                .add(FLD_SPLICE)
                .add(FLD_SYNON)
                .add(FLD_INDEL)
                .toString();
    }

    private static String toString(final LilacAllele allele)
    {
        return new StringJoiner(DELIMITER)
                .add(allele.allele())
                .add(String.valueOf(allele.refFragments()))
                .add(String.valueOf(allele.refUnique()))
                .add(String.valueOf(allele.refShared()))
                .add(String.valueOf(allele.refWild()))
                .add(String.valueOf(allele.tumorFragments()))
                .add(String.valueOf(allele.tumorUnique()))
                .add(String.valueOf(allele.tumorShared()))
                .add(String.valueOf(allele.tumorWild()))
                .add(String.valueOf(allele.rnaFragments()))
                .add(String.valueOf(allele.rnaUnique()))
                .add(String.valueOf(allele.rnaShared()))
                .add(String.valueOf(allele.rnaWild()))
                .add(String.valueOf(allele.tumorCopyNumber()))
                .add(String.valueOf(allele.somaticMissense()))
                .add(String.valueOf(allele.somaticNonsenseOrFrameshift()))
                .add(String.valueOf(allele.somaticSplice()))
                .add(String.valueOf(allele.somaticSynonymous()))
                .add(String.valueOf(allele.somaticInframeIndel()))
                .toString();
    }
}
