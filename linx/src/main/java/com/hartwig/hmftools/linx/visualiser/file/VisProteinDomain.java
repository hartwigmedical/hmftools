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

public class VisProteinDomain implements GenomeRegion
{
    public final String SampleId;
    public final int ClusterId;
    public final String Transcript;
    public String Chromosome;
    public int Start;
    public int End;
    public final String Info;

    public static final String PD_FIVE_PRIME_UTR = "5-Prime UTR";
    public static final String PD_THREE_PRIME_UTR = "3-Prime UTR";
    public static final String PD_NON_CODING = "Non Coding";

    public static final String PROTEIN_DOMAIN_UTR = "UTR/Non-coding";

    public VisProteinDomain(final String sampleId, int clusterId, final String transcript,
            final String chromosome, int start, int end, final String info)
    {
        SampleId = sampleId;
        ClusterId = clusterId;
        Transcript = transcript;
        Chromosome = chromosome;
        Info = info;
        Start = start;
        End = end;
    }

    public static VisProteinDomain from(final VisProteinDomain other)
    {
        return new VisProteinDomain(other.SampleId, other.ClusterId, other.Transcript, other.Chromosome, other.Start, other.End, other.Info);
    }

    @Override
    public String chromosome() { return Chromosome; }

    @Override
    public int start() { return Start; }

    @Override
    public int end() { return End; }

    public String name()
    {
        return Info.equals(PD_FIVE_PRIME_UTR) || Info.equals(PD_FIVE_PRIME_UTR) || Info.equals(PD_NON_CODING) ? PROTEIN_DOMAIN_UTR : Info;
    }

    private static final String FILE_EXTENSION = ".linx.vis_protein_domain.tsv";
    private static final String GERMLINE_FILE_EXTENSION = ".linx.germline.vis_protein_domain.tsv";

    public static String generateFilename(final String basePath, final String sample, boolean isGermline)
    {
        return basePath + File.separator + sample + (isGermline ? GERMLINE_FILE_EXTENSION : FILE_EXTENSION);
    }

    public static List<VisProteinDomain> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(final String filename, List<VisProteinDomain> cnDataList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(cnDataList));
    }

    private static List<String> toLines(final List<VisProteinDomain> cnDataList)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        cnDataList.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    private static List<VisProteinDomain> fromLines(List<String> lines)
    {
        String header = lines.get(0);
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        lines.remove(0);

        List<VisProteinDomain> data = Lists.newArrayList();

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM);

            data.add(new VisProteinDomain(
                    getValue(fieldsIndexMap, FLD_SAMPLE_ID, "", values),
                    getIntValue(fieldsIndexMap, "ClusterId", 0, values),
                    getValue(fieldsIndexMap, "Transcript", "", values),
                    getValue(fieldsIndexMap, "Chromosome", "", values),
                    getIntValue(fieldsIndexMap, "Start", 0, values),
                    getIntValue(fieldsIndexMap, "End", 0, values),
                    getValue(fieldsIndexMap, "Info", "", values)));
        }

        return data;
    }

    public static String header()
    {
        return new StringJoiner(TSV_DELIM)
                .add("ClusterId")
                .add("Transcript")
                .add("Chromosome")
                .add("Start")
                .add("End")
                .add("Info")
                .toString();
    }

    public static String toString(final VisProteinDomain proteinData)
    {
        return new StringJoiner(TSV_DELIM)
                .add(String.valueOf(proteinData.ClusterId))
                .add(String.valueOf(proteinData.Transcript))
                .add(String.valueOf(proteinData.Chromosome))
                .add(String.valueOf(proteinData.Start))
                .add(String.valueOf(proteinData.End))
                .add(String.valueOf(proteinData.Info))
                .toString();
    }
}
