package com.hartwig.hmftools.common.peach;

import static com.hartwig.hmftools.common.peach.PeachUtil.convertZygosityToAlleleCount;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import org.jetbrains.annotations.NotNull;

public final class PeachGenotypeFile
{
    private static final String OLD_FORMAT_HEADER = new StringJoiner(TSV_DELIM).add("gene")
            .add("haplotype")
            .add("function")
            .add("linked_drugs")
            .add("url_prescription_info")
            .add("panel_version")
            .add("repo_version")
            .add("haplotype_only")
            .add("zygosity_only")
            .toString();

    @NotNull
    public static List<PeachGenotype> read(@NotNull final String filename) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
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
    private static List<PeachGenotype> fromLinesJavaPeach(final @NotNull List<String> lines)
    {
        final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);

        int geneIndex = fieldsIndexMap.get("gene");
        int alleleIndex = fieldsIndexMap.get("allele");
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
                    .panelVersion(null)
                    .repoVersion(null)
                    .build();

            genotypes.add(genotype);
        }

        return genotypes;
    }

    @NotNull
    private static List<PeachGenotype> fromLinesPythonPeach(final @NotNull List<String> lines)
    {
        final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);

        int geneIndex = fieldsIndexMap.get("gene");
        int functionIndex = fieldsIndexMap.get("function");
        int linkedDrugsIndex = fieldsIndexMap.get("linked_drugs");
        int urlPrescriptionInfoIndex = fieldsIndexMap.get("url_prescription_info");
        int panelVersionIndex = fieldsIndexMap.get("panel_version");
        int repoVersionIndex = fieldsIndexMap.get("repo_version");
        int alleleIndex = fieldsIndexMap.get("haplotype_only");
        int zygosityIndex = fieldsIndexMap.get("zygosity_only");

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
}
