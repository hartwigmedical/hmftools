package com.hartwig.hmftools.common.hla;

import static com.hartwig.hmftools.common.hla.HlaCommon.MHC_CLASS_I;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_GENES;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.checkFileExtensionRename;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.BufferedWriter;
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
public abstract class LilacAllele
{
    public abstract String genes();

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

    private static final String FILE_EXTENSION = ".lilac.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static List<LilacAllele> read(final String filePath) throws IOException
    {
        String filename = checkFileExtensionRename(filePath);
        String delim = inferFileDelimiter(filename);

        List<String> lines = Files.readAllLines(new File(filename).toPath());

        final String header = lines.get(0);
        lines.remove(0);

        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, delim);

        Integer genesIndex = fieldsIndexMap.get(FLD_GENES);

        List<LilacAllele> alleles = Lists.newArrayList();

        for(String line : lines)
        {
            String[] values = line.split(delim);

            String allele = values[fieldsIndexMap.get(FLD_ALLELE)];
            String genes = genesIndex != null ? values[genesIndex] : MHC_CLASS_I;

            alleles.add(ImmutableLilacAllele.builder()
                    .genes(genes)
                    .allele(allele)
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

    public static void write(final BufferedWriter writer, final List<LilacAllele> alleles) throws IOException
    {
        List<String> lines = Lists.newArrayList(toLines(alleles));
        for(String line : lines)
        {
            writer.write(line);
            writer.newLine();
        }
    }

    private static List<String> toLines(final List<LilacAllele> alleles)
    {
        final List<String> lines = Lists.newArrayList();
        alleles.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    public static final String FLD_ALLELE = "Allele";
    public static final String FLD_REF_TOTAL = "RefTotal";
    public static final String FLD_REF_UNIQUE = "RefUnique";
    public static final String FLD_REF_SHARED = "RefShared";
    public static final String FLD_REF_WILD = "RefWild";
    public static final String FLD_TUMOR_TOTAL = "TumorTotal";
    public static final String FLD_TUMOR_UNIQUE = "TumorUnique";
    public static final String FLD_TUMOR_SHARED = "TumorShared";
    public static final String FLD_TUMOR_WILD = "TumorWild";
    public static final String FLD_TUMOR_CN = "TumorCopyNumber";
    public static final String FLD_RNA_TOTAL = "RnaTotal";
    public static final String FLD_RNA_UNIQUE = "RnaUnique";
    public static final String FLD_RNA_SHARED = "RnaShared";
    public static final String FLD_RNA_WILD = "RnaWild";
    public static final String FLD_MISSENSE = "SomaticMissense";
    public static final String FLD_NFS = "SomaticNonsenseOrFrameshift";
    public static final String FLD_SPLICE = "SomaticSplice";
    public static final String FLD_SYNON = "SomaticSynonymous";
    public static final String FLD_INDEL = "SomaticInframeIndel";

    public static String header()
    {
        return new StringJoiner(TSV_DELIM)
                .add(FLD_GENES)
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
        return new StringJoiner(TSV_DELIM)
                .add(allele.genes())
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
                .add(String.format("%.2f", allele.tumorCopyNumber()))
                .add(String.format("%.1f", allele.somaticMissense()))
                .add(String.format("%.1f", allele.somaticNonsenseOrFrameshift()))
                .add(String.format("%.1f", allele.somaticSplice()))
                .add(String.format("%.1f", allele.somaticSynonymous()))
                .add(String.format("%.1f", allele.somaticInframeIndel()))
                .toString();
    }
}
