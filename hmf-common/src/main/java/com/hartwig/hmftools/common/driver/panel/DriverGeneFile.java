package com.hartwig.hmftools.common.driver.panel;

import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.ANY;
import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.NONE;
import static com.hartwig.hmftools.common.purple.PurpleCommon.DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO;
import static com.hartwig.hmftools.common.purple.PurpleCommon.DEFAULT_DRIVER_HET_DELETION_THRESHOLD;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCategory;

import org.jetbrains.annotations.NotNull;

public final class DriverGeneFile
{
    private static final String OTHER_TRANS_DELIM = ";";

    private DriverGeneFile() {}

    public static void write(final String filename, final List<DriverGene> driverGenes) throws IOException
    {
        List<DriverGene> sorted = Lists.newArrayList(driverGenes);
        Files.write(new File(filename).toPath(), toLines(sorted));
    }

    public static List<DriverGene> read(final String filename) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    public static List<DriverGene> fromLines(final List<String> lines)
    {
        List<DriverGene> driverGenes = Lists.newArrayList();

        String header = lines.get(0);
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        lines.remove(0);

        int geneIndex = fieldsIndexMap.get("gene");
        int missenseIndex = fieldsIndexMap.get("reportMissense");
        int nonsenseIndex = fieldsIndexMap.get("reportNonsense");
        int spliceIndex = fieldsIndexMap.get("reportSplice");
        int deletionIndex = fieldsIndexMap.get("reportDeletion");
        Integer hetDelIndex = fieldsIndexMap.get("reportHetDeletion");
        Integer hetDelThresholdIndex = fieldsIndexMap.get("hetDeletionThreshold");
        int disruptionIndex = fieldsIndexMap.get("reportDisruption");
        int amplificationIndex = fieldsIndexMap.get("reportAmplification");
        Integer ampRatioIndex = fieldsIndexMap.get("amplificationRatio");

        int somaticHotspotIndex = fieldsIndexMap.get("reportSomaticHotspot");

        int likelihoodTypeIndex = fieldsIndexMap.get("likelihoodType");
        int germlineVariantIndex = fieldsIndexMap.get("reportGermlineVariant");
        int germlineHotspotIndex = fieldsIndexMap.get("reportGermlineHotspot");
        int germlineDisruptionIndex = fieldsIndexMap.get("reportGermlineDisruption");
        int germlineDeletionIndex = fieldsIndexMap.get("reportGermlineDeletion");
        int altTransIndex = fieldsIndexMap.get("additionalReportedTranscripts");
        int reportPGXIndex = fieldsIndexMap.get("reportPGX");

        ImmutableDriverGene.Builder builder = ImmutableDriverGene.builder();

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            List<String> otherReportableTrans = Arrays.stream(values[altTransIndex].split(OTHER_TRANS_DELIM))
                    .filter(x -> !x.isEmpty())
                    .collect(Collectors.toList());

            // backwards compatibility prior to pipeline v5.32
            String reportGermlineDisruptionStr = values[germlineDisruptionIndex];
            DriverGeneGermlineReporting reportGermlineDisruption;

            if(reportGermlineDisruptionStr.toLowerCase().equals(Boolean.TRUE.toString()))
                reportGermlineDisruption = ANY;
            else if(reportGermlineDisruptionStr.toLowerCase().equals(Boolean.FALSE.toString()))
                reportGermlineDisruption = NONE;
            else
                reportGermlineDisruption = DriverGeneGermlineReporting.valueOf(reportGermlineDisruptionStr);

            DriverGeneGermlineReporting reportGermlineDeletion = DriverGeneGermlineReporting.valueOf(values[germlineDeletionIndex]);

            builder.gene(values[geneIndex])
                    .reportMissenseAndInframe(Boolean.parseBoolean(values[missenseIndex]))
                    .reportNonsenseAndFrameshift(Boolean.parseBoolean(values[nonsenseIndex]))
                    .reportSplice(Boolean.parseBoolean(values[spliceIndex]))
                    .reportDeletion(Boolean.parseBoolean(values[deletionIndex]))
                    .reportHetDeletion(hetDelIndex != null ? Boolean.parseBoolean(values[hetDelIndex]) : false)
                    .hetDeletionThreshold(hetDelThresholdIndex != null ?
                            Double.parseDouble(values[hetDelThresholdIndex]) : DEFAULT_DRIVER_HET_DELETION_THRESHOLD)
                    .reportDisruption(Boolean.parseBoolean(values[disruptionIndex]))
                    .reportAmplification(Boolean.parseBoolean(values[amplificationIndex]))
                    .amplificationRatio(ampRatioIndex != null ?
                            Double.parseDouble(values[ampRatioIndex]) : DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO)
                    .reportSomaticHotspot(Boolean.parseBoolean(values[somaticHotspotIndex]))
                    .likelihoodType(DriverCategory.valueOf(values[likelihoodTypeIndex]))
                    .reportGermlineVariant(DriverGeneGermlineReporting.valueOf(values[germlineVariantIndex].toUpperCase()))
                    .reportGermlineHotspot(DriverGeneGermlineReporting.valueOf(values[germlineHotspotIndex].toUpperCase()))
                    .reportGermlineDisruption(reportGermlineDisruption)
                    .reportGermlineDeletion(reportGermlineDeletion)
                    .additionalReportedTranscripts(otherReportableTrans)
                    .reportPGX(Boolean.parseBoolean(values[reportPGXIndex]));

            driverGenes.add(builder.build());
        }

        return driverGenes;
    }

    private static String otherReportableTransStr(final List<String> otherTrans)
    {
        if(otherTrans.isEmpty())
            return "";

        StringJoiner sj = new StringJoiner(OTHER_TRANS_DELIM);
        otherTrans.forEach(x -> sj.add(x));
        return sj.toString();
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIM).add("gene")
                .add("reportMissense")
                .add("reportNonsense")
                .add("reportSplice")
                .add("reportDeletion")
                .add("reportHetDeletion")
                .add("hetDeletionThreshold")
                .add("reportDisruption")
                .add("reportAmplification")
                .add("amplificationRatio")
                .add("reportSomaticHotspot")
                .add("likelihoodType")
                .add("reportGermlineVariant")
                .add("reportGermlineHotspot")
                .add("reportGermlineDisruption")
                .add("reportGermlineDeletion")
                .add("additionalReportedTranscripts")
                .add("reportPGX")
                .toString();
    }

    private static String toString(final DriverGene gene)
    {
        return new StringJoiner(TSV_DELIM)
                .add(gene.gene())
                .add(String.valueOf(gene.reportMissenseAndInframe()))
                .add(String.valueOf(gene.reportNonsenseAndFrameshift()))
                .add(String.valueOf(gene.reportSplice()))
                .add(String.valueOf(gene.reportDeletion()))
                .add(String.valueOf(gene.reportHetDeletion()))
                .add(String.valueOf(gene.hetDeletionThreshold()))
                .add(String.valueOf(gene.reportDisruption()))
                .add(String.valueOf(gene.reportAmplification()))
                .add(String.valueOf(gene.amplificationRatio()))
                .add(String.valueOf(gene.reportSomaticHotspot()))
                .add(String.valueOf(gene.likelihoodType()))
                .add(String.valueOf(gene.reportGermlineVariant()))
                .add(String.valueOf(gene.reportGermlineHotspot()))
                .add(String.valueOf(gene.reportGermlineDisruption()))
                .add(String.valueOf(gene.reportGermlineDeletion()))
                .add(otherReportableTransStr(gene.additionalReportedTranscripts()))
                .add(String.valueOf(gene.reportPGX()))
                .toString();
    }

    private static List<String> toLines(final List<DriverGene> driverGenes)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        driverGenes.stream().map(DriverGeneFile::toString).forEach(lines::add);
        return lines;
    }
}
