package com.hartwig.hmftools.common.peach;

import static com.hartwig.hmftools.common.peach.PeachUtil.convertZygosityToAlleleCount;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;

public final class PeachGenotypeFile
{
    private static final String PYTHON_GENE_COLUMN_NAME = "gene";
    private static final String PYTHON_HAPLOTYPE_AND_ZYGOSITY_COLUMN_NAME = "haplotype";
    private static final String PYTHON_FUNCTION_COLUMN_NAME = "function";
    private static final String PYTHON_DRUGS_COLUMN_NAME = "linked_drugs";
    private static final String PYTHON_PRESCRIPTION_URL_COLUMN_NAME = "url_prescription_info";
    private static final String PYTHON_PANEL_VERSION_COLUMN_NAME = "panel_version";
    private static final String PYTHON_REPO_VERSION_COLUMN_NAME = "repo_version";
    private static final String PYTHON_HAPLOTYPE_ONLY_COLUMN_NAME = "haplotype_only";
    private static final String PYTHON_ZYGOSITY_ONLY_COLUMN_NAME = "zygosity_only";

    private static final String OLD_FORMAT_HEADER = new StringJoiner(TSV_DELIM).add(PYTHON_GENE_COLUMN_NAME)
            .add(PYTHON_HAPLOTYPE_AND_ZYGOSITY_COLUMN_NAME)
            .add(PYTHON_FUNCTION_COLUMN_NAME)
            .add(PYTHON_DRUGS_COLUMN_NAME)
            .add(PYTHON_PRESCRIPTION_URL_COLUMN_NAME)
            .add(PYTHON_PANEL_VERSION_COLUMN_NAME)
            .add(PYTHON_REPO_VERSION_COLUMN_NAME)
            .add(PYTHON_HAPLOTYPE_ONLY_COLUMN_NAME)
            .add(PYTHON_ZYGOSITY_ONLY_COLUMN_NAME)
            .toString();

    private static final String FILE_EXTENSION = ".peach.haplotypes.best.tsv";

    @NotNull
    public static String generateFileName(@NotNull String outputDir, @NotNull String sampleId)
    {
        return checkAddDirSeparator(outputDir) + sampleId + FILE_EXTENSION;
    }

    @NotNull
    public static List<PeachGenotype> read(@NotNull final String filename) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull final List<PeachGenotype> genotypes) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(genotypes));
    }

    @NotNull
    private static List<PeachGenotype> fromLines(@NotNull final List<String> lines)
    {
        final String header = lines.get(0);
        if(header.equals(OLD_FORMAT_HEADER))
        {
            return fromLinesPythonPeach(lines);
        }
        else
        {
            return fromLinesJavaPeach(lines);
        }
    }

    @NotNull
    private static List<PeachGenotype> fromLinesJavaPeach(@NotNull final List<String> lines)
    {
        final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);

        int geneIndex = fieldsIndexMap.get("gene");
        int alleleIndex = fieldsIndexMap.get("haplotype");
        int alleleCountIndex = fieldsIndexMap.get("count");
        int functionIndex = fieldsIndexMap.get("function");
        int linkedDrugsIndex = fieldsIndexMap.get("linkedDrugs");
        int urlPrescriptionInfoIndex = fieldsIndexMap.get("prescriptionUrls");

        List<PeachGenotype> genotypes = new ArrayList<>();
        lines.remove(0);
        for(final String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            final PeachGenotype genotype = ImmutablePeachGenotype.builder()
                    .gene(values[geneIndex])
                    .allele(values[alleleIndex])
                    .alleleCount(Integer.parseInt(values[alleleCountIndex]))
                    .function(values[functionIndex])
                    .linkedDrugs(values[linkedDrugsIndex])
                    .urlPrescriptionInfo(values[urlPrescriptionInfoIndex])
                    .build();

            genotypes.add(genotype);
        }

        return genotypes;
    }

    @NotNull
    private static List<PeachGenotype> fromLinesPythonPeach(@NotNull final List<String> lines)
    {
        final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);

        int geneIndex = fieldsIndexMap.get(PYTHON_GENE_COLUMN_NAME);
        int functionIndex = fieldsIndexMap.get(PYTHON_FUNCTION_COLUMN_NAME);
        int linkedDrugsIndex = fieldsIndexMap.get(PYTHON_DRUGS_COLUMN_NAME);
        int urlPrescriptionInfoIndex = fieldsIndexMap.get(PYTHON_PRESCRIPTION_URL_COLUMN_NAME);
        int panelVersionIndex = fieldsIndexMap.get(PYTHON_PANEL_VERSION_COLUMN_NAME);
        int repoVersionIndex = fieldsIndexMap.get(PYTHON_REPO_VERSION_COLUMN_NAME);
        int alleleIndex = fieldsIndexMap.get(PYTHON_HAPLOTYPE_ONLY_COLUMN_NAME);
        int zygosityIndex = fieldsIndexMap.get(PYTHON_ZYGOSITY_ONLY_COLUMN_NAME);

        List<PeachGenotype> genotypes = new ArrayList<>();
        lines.remove(0);
        for(final String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            final PeachGenotype genotype = ImmutablePeachGenotype.builder()
                    .gene(values[geneIndex])
                    .allele(values[alleleIndex])
                    .alleleCount(convertZygosityToAlleleCount(values[zygosityIndex]))
                    .function(values[functionIndex])
                    .linkedDrugs(values[linkedDrugsIndex])
                    .urlPrescriptionInfo(values[urlPrescriptionInfoIndex])
                    .panelVersion(values[panelVersionIndex])
                    .repoVersion(values[repoVersionIndex])
                    .build();

            genotypes.add(genotype);
        }

        return genotypes;
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull final List<PeachGenotype> genotypes)
    {
        List<String> lines = new ArrayList<>();
        lines.add(headerJavaPeach());
        genotypes.stream().map(PeachGenotypeFile::toLine).forEachOrdered(lines::add);
        return lines;
    }

    @NotNull
    private static String headerJavaPeach()
    {
        return new StringJoiner(TSV_DELIM).add("gene")
                .add("haplotype")
                .add("count")
                .add("function")
                .add("linkedDrugs")
                .add("prescriptionUrls")
                .toString();
    }

    @NotNull
    private static String toLine(@NotNull final PeachGenotype genotype)
    {
        return new StringJoiner(TSV_DELIM).add(genotype.gene())
                .add(genotype.allele())
                .add(Integer.toString(genotype.alleleCount()))
                .add(genotype.function())
                .add(genotype.linkedDrugs())
                .add(genotype.urlPrescriptionInfo())
                .toString();
    }
}
