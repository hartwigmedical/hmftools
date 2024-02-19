package com.hartwig.hmftools.linx.visualiser.file;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getIntValue;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getValue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegionBuilder;

public class VisGeneExon implements GenomeRegion
{
    public final String SampleId;
    public final int ClusterId;
    public final String Gene;
    public final String Transcript;
    public final String Chromosome;
    public final VisGeneAnnotationType AnnotationType;
    public final int ExonRank;
    public final int ExonStart;
    public final int ExonEnd;

    public VisGeneExon(
            final String sampleId, int clusterId, final String gene, final String transcript, final String chromosome,
            final VisGeneAnnotationType type, int exonRank, int exonStart, int exonEnd)
    {
        SampleId = sampleId;
        ClusterId = clusterId;
        Gene = gene;
        Transcript = transcript;
        Chromosome = chromosome;
        AnnotationType = type;
        ExonRank = exonRank;
        ExonStart = exonStart;
        ExonEnd = exonEnd;
    }

    public abstract static class Builder implements GenomeRegionBuilder<VisGeneExon> { }

    @Override
    public String chromosome() { return Chromosome; }

    @Override
    public int start() { return ExonStart; }

    @Override
    public int end() { return ExonEnd; }

    private static final String FILE_EXTENSION = ".linx.vis_gene_exon.tsv";
    private static final String GERMLINE_FILE_EXTENSION = ".linx.germline.vis_gene_exon.tsv";

    public static String generateFilename(final String basePath, final String sample, boolean isGermline)
    {
        return basePath + File.separator + sample + (isGermline ? GERMLINE_FILE_EXTENSION : FILE_EXTENSION);
    }

    public static List<VisGeneExon> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(final String filename, List<VisGeneExon> cnDataList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(cnDataList));
    }

    private static List<String> toLines(final List<VisGeneExon> cnDataList)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        cnDataList.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    private static List<VisGeneExon> fromLines(final List<String> lines)
    {
        String header = lines.get(0);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        lines.remove(0);

        List<VisGeneExon> data = Lists.newArrayList();

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM);

            data.add(new VisGeneExon(
                    getValue(fieldsIndexMap, FLD_SAMPLE_ID, "", values),
                    getIntValue(fieldsIndexMap, "ClusterId", 0, values),
                    getValue(fieldsIndexMap, "Gene", "", values),
                    getValue(fieldsIndexMap, "Transcript", "", values),
                    getValue(fieldsIndexMap, "Chromosome", "", values),
                    VisGeneAnnotationType.valueOf(values[fieldsIndexMap.get("AnnotationType")]),
                    getIntValue(fieldsIndexMap, "ExonRank", 0, values),
                    getIntValue(fieldsIndexMap, "ExonStart", 0, values),
                    getIntValue(fieldsIndexMap, "ExonEnd", 0, values)));
        }

        return data;
    }

    public static String header()
    {
        return new StringJoiner(TSV_DELIM)
                .add("ClusterId")
                .add("Gene")
                .add("Transcript")
                .add("Chromosome")
                .add("AnnotationType")
                .add("ExonRank")
                .add("ExonStart")
                .add("ExonEnd")
                .toString();
    }

    public static String toString(final VisGeneExon geData)
    {
        return new StringJoiner(TSV_DELIM)
                .add(String.valueOf(geData.ClusterId))
                .add(String.valueOf(geData.Gene))
                .add(String.valueOf(geData.Transcript))
                .add(String.valueOf(geData.Chromosome))
                .add(String.valueOf(geData.AnnotationType))
                .add(String.valueOf(geData.ExonRank))
                .add(String.valueOf(geData.ExonStart))
                .add(String.valueOf(geData.ExonEnd))
                .toString();
    }
}
