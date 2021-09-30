package com.hartwig.hmftools.common.drivercatalog.panel;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;

import org.jetbrains.annotations.NotNull;

public final class DriverGeneFile
{
    private static final String DELIMITER = "\t";

    private DriverGeneFile()
    {
    }

    public static void write(@NotNull final String filename, @NotNull final List<DriverGene> driverGenes) throws IOException
    {
        List<DriverGene> sorted = Lists.newArrayList(driverGenes);
        Files.write(new File(filename).toPath(), toLines(sorted));
    }

    @NotNull
    public static List<DriverGene> read(@NotNull final String filename) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    @NotNull
    public static List<DriverGene> fromLines(@NotNull final List<String> lines)
    {
        List<DriverGene> driverGenes = Lists.newArrayList();

        String header = lines.get(0);
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIMITER);
        lines.remove(0);

        int geneIndex = fieldsIndexMap.get("gene");
        int missenseIndex = fieldsIndexMap.get("reportMissense");
        int nonsenseIndex = fieldsIndexMap.get("reportNonsense");
        int spliceIndex = fieldsIndexMap.get("reportSplice");
        int deletionIndex = fieldsIndexMap.get("reportDeletion");
        int disruptionIndex = fieldsIndexMap.get("reportDisruption");
        int amplificationIndex = fieldsIndexMap.get("reportAmplification");

        int somaticHotspotIndex = fieldsIndexMap.containsKey("reportSomaticHotspot") ?
                fieldsIndexMap.get("reportSomaticHotspot") : fieldsIndexMap.get("reportHotspot"); // for older files

        int likelihoodTypeIndex = fieldsIndexMap.get("likelihoodType");
        Integer germlineVariantIndex = fieldsIndexMap.get("reportGermlineVariant");
        Integer germlineHotspotIndex = fieldsIndexMap.get("reportGermlineHotspot");
        Integer germlineDisruptionIndex = fieldsIndexMap.get("reportGermlineDisruption");
        Integer altTransIndex = fieldsIndexMap.get("additionalReportedTranscripts");

        ImmutableDriverGene.Builder builder = ImmutableDriverGene.builder();

        for(String line : lines)
        {
            String[] values = line.split(DELIMITER, -1);

            builder.gene(values[geneIndex])
                    .reportMissenseAndInframe(Boolean.parseBoolean(values[missenseIndex]))
                    .reportNonsenseAndFrameshift(Boolean.parseBoolean(values[nonsenseIndex]))
                    .reportSplice(Boolean.parseBoolean(values[spliceIndex]))
                    .reportDeletion(Boolean.parseBoolean(values[deletionIndex]))
                    .reportDisruption(Boolean.parseBoolean(values[disruptionIndex]))
                    .reportAmplification(Boolean.parseBoolean(values[amplificationIndex]))
                    .reportSomaticHotspot(Boolean.parseBoolean(values[somaticHotspotIndex]))
                    .likelihoodType(DriverCategory.valueOf(values[likelihoodTypeIndex]))
                    .reportGermlineVariant(germlineVariantIndex != null
                            ?
                            DriverGeneGermlineReporting.valueOf(values[germlineVariantIndex].toUpperCase())
                            : DriverGeneGermlineReporting.NONE)
                    .reportGermlineHotspot(germlineHotspotIndex != null
                            ?
                            DriverGeneGermlineReporting.valueOf(values[germlineHotspotIndex].toUpperCase())
                            : DriverGeneGermlineReporting.NONE)
                    .reportGermlineDisruption(germlineDisruptionIndex != null ?
                            Boolean.parseBoolean(values[germlineDisruptionIndex]) : Boolean.parseBoolean(values[disruptionIndex]))
                    .additionalReportedTranscripts(altTransIndex != null ? values[altTransIndex] : "");

            driverGenes.add(builder.build());
        }

        return driverGenes;
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER).add("gene")
                .add("reportMissense")
                .add("reportNonsense")
                .add("reportSplice")
                .add("reportDeletion")
                .add("reportDisruption")
                .add("reportAmplification")
                .add("reportSomaticHotspot")
                .add("likelihoodType")
                .add("reportGermlineVariant")
                .add("reportGermlineHotspot")
                .add("reportGermlineDisruption")
                .add("additionalReportedTranscripts")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final DriverGene gene)
    {
        return new StringJoiner(DELIMITER)
                .add(gene.gene())
                .add(String.valueOf(gene.reportMissenseAndInframe()))
                .add(String.valueOf(gene.reportNonsenseAndFrameshift()))
                .add(String.valueOf(gene.reportSplice()))
                .add(String.valueOf(gene.reportDeletion()))
                .add(String.valueOf(gene.reportDisruption()))
                .add(String.valueOf(gene.reportAmplification()))
                .add(String.valueOf(gene.reportSomaticHotspot()))
                .add(String.valueOf(gene.likelihoodType()))
                .add(String.valueOf(gene.reportGermlineVariant()))
                .add(String.valueOf(gene.reportGermlineHotspot()))
                .add(String.valueOf(gene.reportGermlineDisruption()))
                .add(gene.additionalReportedTranscripts())
                .toString();
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<DriverGene> driverGenes)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        driverGenes.stream().map(DriverGeneFile::toString).forEach(lines::add);
        return lines;
    }
}
