package com.hartwig.hmftools.common.virus;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.virus.VirusLikelihoodType.UNKNOWN;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class AnnotatedVirusFile
{
    private static final String ANNOTATED_VIRUS_EXTENSION = ".virus.annotated.tsv";

    @NotNull
    public static String generateFileName(@NotNull String outputDir, @NotNull String sampleId)
    {
        return checkAddDirSeparator(outputDir) + sampleId + ANNOTATED_VIRUS_EXTENSION;
    }

    @NotNull
    public static List<AnnotatedVirus> read(@NotNull String annotatedVirusTsv) throws IOException
    {
        return fromLines(Files.readAllLines(new File(annotatedVirusTsv).toPath()));
    }

    public static void write(@NotNull String annotatedVirusTsv, @NotNull List<AnnotatedVirus> annotatedViruses) throws IOException
    {
        Files.write(new File(annotatedVirusTsv).toPath(), toLines(annotatedViruses));
    }

    @VisibleForTesting
    @NotNull
    static List<String> toLines(@NotNull List<AnnotatedVirus> annotatedViruses)
    {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        annotatedViruses.stream().map(AnnotatedVirusFile::toString).forEach(lines::add);
        return lines;
    }

    @VisibleForTesting
    @NotNull
    static List<AnnotatedVirus> fromLines(@NotNull List<String> lines)
    {
        List<AnnotatedVirus> virusList = Lists.newArrayList();

        String header = lines.get(0);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        lines.remove(0);

        // support for version 1 not having coverage columns
        Integer percentageCoveredIndex = fieldsIndexMap.get("percentageCovered");
        Integer meanCoverageIndex = fieldsIndexMap.get("meanCoverage");
        Integer expectedClonalCoverageIndex = fieldsIndexMap.get("expectedClonalCoverage");
        Integer driverLikelihoodIndex = fieldsIndexMap.get("driverLikelihood");

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            Double expectedClonalCoverage = null;
            if(expectedClonalCoverageIndex != null)
            {
                if(!values[expectedClonalCoverageIndex].equals("null"))
                {
                    expectedClonalCoverage = Double.valueOf(values[expectedClonalCoverageIndex]);
                }
            }

            virusList.add(ImmutableAnnotatedVirus.builder()
                    .taxid(Integer.parseInt(values[fieldsIndexMap.get("taxid")]))
                    .name(values[fieldsIndexMap.get("name")])
                    .qcStatus(VirusBreakendQCStatus.valueOf(values[fieldsIndexMap.get("qcStatus")]))
                    .integrations(Integer.parseInt(values[fieldsIndexMap.get("integrations")]))
                    .interpretation(values[fieldsIndexMap.get("interpretation")].equals("null")
                            ? null
                            : VirusType.fromVirusName(values[fieldsIndexMap.get("interpretation")]))
                    .percentageCovered(percentageCoveredIndex != null ? Double.parseDouble(values[percentageCoveredIndex]) : 0)
                    .meanCoverage(meanCoverageIndex != null ? Double.parseDouble(values[meanCoverageIndex]) : 0)
                    .expectedClonalCoverage(expectedClonalCoverage)
                    .reported(Boolean.parseBoolean(values[fieldsIndexMap.get("reported")]))
                    .blacklisted(Boolean.parseBoolean(values[fieldsIndexMap.get("blacklisted")]))
                    .virusDriverLikelihoodType(driverLikelihoodIndex != null ?
                            VirusLikelihoodType.valueOf(values[driverLikelihoodIndex]) : UNKNOWN)
                    .build());
        }
        return virusList;
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(TSV_DELIM).add("taxid")
                .add("name")
                .add("qcStatus")
                .add("integrations")
                .add("interpretation")
                .add("percentageCovered")
                .add("meanCoverage")
                .add("expectedClonalCoverage")
                .add("reported")
                .add("blacklisted")
                .add("driverLikelihood")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull AnnotatedVirus annotatedVirus)
    {
        return new StringJoiner(TSV_DELIM).add(String.valueOf(annotatedVirus.taxid()))
                .add(annotatedVirus.name())
                .add(annotatedVirus.qcStatus().toString())
                .add(String.valueOf(annotatedVirus.integrations()))
                .add(String.valueOf(annotatedVirus.interpretation()))
                .add(String.valueOf(annotatedVirus.percentageCovered()))
                .add(String.valueOf(annotatedVirus.meanCoverage()))
                .add(String.valueOf(annotatedVirus.expectedClonalCoverage()))
                .add(String.valueOf(annotatedVirus.reported()))
                .add(String.valueOf(annotatedVirus.blacklisted()))
                .add(String.valueOf(annotatedVirus.virusDriverLikelihoodType()))
                .toString();
    }
}