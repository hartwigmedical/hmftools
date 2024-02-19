package com.hartwig.hmftools.linx.visualiser.file;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getDoubleValue;
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

public class VisCopyNumber implements GenomeRegion
{
    public final String SampleId;
    public final String Chromosome;
    public int Start;
    public int End;
    public final double CopyNumber;
    public final double BAF;
    public boolean Truncated;

    public VisCopyNumber(final String sampleId, final String chromosome, int start, int end, double copyNumber, double baf)
    {
        SampleId = sampleId;
        Chromosome = chromosome;
        Start = start;
        End = end;
        CopyNumber = copyNumber;
        BAF = baf;
        Truncated = false;
    }

    public static VisCopyNumber from(final VisCopyNumber other)
    {
        VisCopyNumber newCn = new VisCopyNumber(
                other.SampleId, other.Chromosome, other.Start, other.End, other.CopyNumber, other.BAF);
        newCn.Truncated = other.Truncated;
        return newCn;
    }

    @Override
    public String chromosome() { return Chromosome; }

    @Override
    public int start() { return Start; }

    @Override
    public int end() { return End; }

    public double minorAlleleCopyNumber()
    {
        return Math.max(0, (1 - BAF) * CopyNumber);
    }

    private static final String FILE_EXTENSION = ".linx.vis_copy_number.tsv";
    private static final String GERMLINE_FILE_EXTENSION = ".linx.germline.vis_copy_number.tsv";

    public static String generateFilename(final String basePath, final String sample, boolean isGermline)
    {
        return basePath + File.separator + sample + (isGermline ? GERMLINE_FILE_EXTENSION : FILE_EXTENSION);
    }

    public static List<VisCopyNumber> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(final String filename, List<VisCopyNumber> cnDataList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(cnDataList));
    }

    private static List<String> toLines(final List<VisCopyNumber> cnDataList)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        cnDataList.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    private static List<VisCopyNumber> fromLines(final List<String> lines)
    {
        String header = lines.get(0);
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        lines.remove(0);

        List<VisCopyNumber> data = Lists.newArrayList();

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM);

            data.add(new VisCopyNumber(
                    getValue(fieldsIndexMap, FLD_SAMPLE_ID, "", values),
                    getValue(fieldsIndexMap, FLD_CHROMOSOME, "", values),
                    getIntValue(fieldsIndexMap, "Start", 0, values),
                    getIntValue(fieldsIndexMap, "End", 0, values),
                    getDoubleValue(fieldsIndexMap, "CopyNumber", 0, values),
                    getDoubleValue(fieldsIndexMap, "BAF", 0, values)));
        }

        return data;
    }

    public static String header()
    {
        return new StringJoiner(TSV_DELIM)
                .add(FLD_CHROMOSOME)
                .add("Start")
                .add("End")
                .add("CopyNumber")
                .add("BAF")
                .toString();
    }

    public static String toString(final VisCopyNumber cnData)
    {
        return new StringJoiner(TSV_DELIM)
                .add(String.valueOf(cnData.Chromosome))
                .add(String.valueOf(cnData.Start))
                .add(String.valueOf(cnData.End))
                .add(VisDataWriter.FORMAT.format(cnData.CopyNumber))
                .add(VisDataWriter.FORMAT.format(cnData.BAF))
                .toString();
    }
}
