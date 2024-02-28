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

public class VisFusion
{
    public final String SampleId;
    public final boolean Reportable;
    public final int ClusterId;

    public final String GeneNameUp;
    public final String TranscriptUp;
    public final String ChrUp;
    public final int PosUp;
    public final int StrandUp;
    public final String RegionTypeUp;
    public final int FusedExonUp;

    public final String GeneNameDown;
    public final String TranscriptDown;
    public final String ChrDown;
    public final int PosDown;
    public final int StrandDown;
    public final String RegionTypeDown;
    public final int FusedExonDown;

    public VisFusion(final String sampleId, int clusterId, boolean reportable,
            final String geneNameUp, final String transcriptUp,
            final String chrUp, int posUp, int strandUp, final String regionTypeUp, int fusedExonUp,
            final String geneNameDown, final String transcriptDown,
            final String chrDown, int posDown, int strandDown, final String regionTypeDown, int fusedExonDown)
    {
        SampleId = sampleId;
        ClusterId = clusterId;
        Reportable = reportable;

        GeneNameUp = geneNameUp;
        TranscriptUp = transcriptUp;
        ChrUp = chrUp;
        PosUp = posUp;
        StrandUp = strandUp;
        RegionTypeUp = regionTypeUp;
        FusedExonUp = fusedExonUp;

        GeneNameDown = geneNameDown;
        TranscriptDown = transcriptDown;
        ChrDown = chrDown;
        PosDown = posDown;
        StrandDown = strandDown;
        RegionTypeDown = regionTypeDown;
        FusedExonDown = fusedExonDown;
    }

    public String name()
    {
        return GeneNameUp + "_" + GeneNameDown;
    }

    private static final String FILE_EXTENSION = ".linx.vis_fusion.tsv";
    private static final String GERMLINE_FILE_EXTENSION = ".linx.germline.vis_fusion.tsv";

    public static String generateFilename(final String basePath, final String sample, boolean isGermline)
    {
        return basePath + File.separator + sample + (isGermline ? GERMLINE_FILE_EXTENSION : FILE_EXTENSION);
    }

    public static List<VisFusion> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(final String filename, List<VisFusion> dataList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(dataList));
    }

    private static List<String> toLines(final List<VisFusion> dataList)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        dataList.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    private static List<VisFusion> fromLines(List<String> lines)
    {
        String header = lines.get(0);
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        lines.remove(0);

        List<VisFusion> data = Lists.newArrayList();

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM);

            data.add(new VisFusion(
                    getValue(fieldsIndexMap, FLD_SAMPLE_ID, "", values),
                    getIntValue(fieldsIndexMap, "ClusterId", 0, values),
                    Boolean.parseBoolean(values[fieldsIndexMap.get("Reportable")]),
                    getValue(fieldsIndexMap, "GeneNameUp", "", values),
                    getValue(fieldsIndexMap, "TranscriptUp", "", values),
                    getValue(fieldsIndexMap, "ChrUp", "", values),
                    getIntValue(fieldsIndexMap, "PosUp", 0, values),
                    getIntValue(fieldsIndexMap, "StrandUp", 0, values),
                    getValue(fieldsIndexMap, "RegionTypeUp", "", values),
                    getIntValue(fieldsIndexMap, "FusedExonUp", 0, values),
                    getValue(fieldsIndexMap, "GeneNameDown", "", values),
                    getValue(fieldsIndexMap, "TranscriptDown", "", values),
                    getValue(fieldsIndexMap, "ChrDown", "", values),
                    getIntValue(fieldsIndexMap, "PosDown", 0, values),
                    getIntValue(fieldsIndexMap, "StrandDown", 0, values),
                    getValue(fieldsIndexMap, "RegionTypeDown", "", values),
                    getIntValue(fieldsIndexMap, "FusedExonDown", 0, values)));
        }

        return data;
    }

    public static String header()
    {
        return new StringJoiner(TSV_DELIM)
                .add("ClusterId")
                .add("Reportable")
                .add("GeneNameUp")
                .add("TranscriptUp")
                .add("ChrUp")
                .add("PosUp")
                .add("StrandUp")
                .add("RegionTypeUp")
                .add("FusedExonUp")
                .add("GeneNameDown")
                .add("TranscriptDown")
                .add("ChrDown")
                .add("PosDown")
                .add("StrandDown")
                .add("RegionTypeDown")
                .add("FusedExonDown")
                .toString();
    }

    public static String toString(final VisFusion data)
    {
        return new StringJoiner(TSV_DELIM)
                .add(String.valueOf(data.ClusterId))
                .add(String.valueOf(data.Reportable))
                .add(String.valueOf(data.GeneNameUp))
                .add(String.valueOf(data.TranscriptUp))
                .add(String.valueOf(data.ChrUp))
                .add(String.valueOf(data.PosUp))
                .add(String.valueOf(data.StrandUp))
                .add(String.valueOf(data.RegionTypeUp))
                .add(String.valueOf(data.FusedExonUp))
                .add(String.valueOf(data.GeneNameDown))
                .add(String.valueOf(data.TranscriptDown))
                .add(String.valueOf(data.ChrDown))
                .add(String.valueOf(data.PosDown))
                .add(String.valueOf(data.StrandDown))
                .add(String.valueOf(data.RegionTypeDown))
                .add(String.valueOf(data.FusedExonDown))
                .toString();
    }
}
