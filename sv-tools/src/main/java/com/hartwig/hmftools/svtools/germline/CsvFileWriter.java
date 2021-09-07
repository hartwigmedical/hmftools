package com.hartwig.hmftools.svtools.germline;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.svtools.germline.GermlineUtils.GM_LOGGER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class CsvFileWriter
{
    private final GermlineVcfConfig mConfig;
    private BufferedWriter mCsvWriter;

    private static final String FILE_EXTENSION = ".linx.germline_sv.tsv";
    private static final String DELIMITER = ",";
    public static final int FIELD_COUNT = 25;

    public CsvFileWriter(final GermlineVcfConfig config)
    {
        mConfig = config;
        mCsvWriter = null;
    }

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<SvData> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<SvData> svList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(svList));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<SvData> svList)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        svList.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static List<SvData> fromLines(@NotNull List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("SampleId")).map(x -> fromString(x, false)).collect(toList());
    }

    @NotNull
    public static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("SampleId")
                .add("Id")
                .add("GridssFilter")
                .add("QualScore")
                .add("Type")
                .add("ChrStart")
                .add("ChrEnd")
                .add("PosStart")
                .add("PosEnd")
                .add("OrientStart")
                .add("OrientEnd")
                .add("NormalREF")
                .add("NormalRP")
                .add("NormalRPQ")
                .add("NormalSR")
                .add("NormalSRQ")
                .add("NormalVF")
                .add("InsertSequence")
                .add("Homology")
                .add("AffectedGenes")
                .add("AsmbStart")
                .add("AsmbEnd")
                .toString();
    }

    @NotNull
    public static String toString(@NotNull final SvData sv)
    {
        /*
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(sv.SampleId))
                .add(String.valueOf(sv.Id))
                .add(String.valueOf(sv.GridssFilter))
                .add(String.format("%.1f", sv.QualScore))
                .add(String.valueOf(sv.Type))
                .add(String.valueOf(sv.ChrStart))
                .add(String.valueOf(sv.ChrEnd))
                .add(String.valueOf(sv.PosStart))
                .add(String.valueOf(sv.PosEnd))
                .add(String.valueOf(sv.OrientStart))
                .add(String.valueOf(sv.OrientEnd))
                .add(String.valueOf(sv.NormalREF))
                .add(String.valueOf(sv.NormalRP))
                .add(String.format("%.2f", sv.NormalRPQ))
                .add(String.valueOf(sv.NormalSR))
                .add(String.format("%.2f", sv.NormalSRQ))
                .add(String.valueOf(sv.NormalVF))
                .add(String.valueOf(sv.InsertSequence))
                .add(String.valueOf(sv.Homology))
                .add(String.valueOf(sv.assemblySvIds(true)))
                .add(String.valueOf(sv.assemblySvIds(false)))
                .toString();

         */

        return "";
    }

    @NotNull
    public static SvData fromString(@NotNull final String sv, boolean reprocess)
    {
        String[] values = sv.split(DELIMITER, -1);
        int index = 0;

        return null;

        /*
        return new GermlineSV(
                values[index++],
                values[index++],
                values[index++],
                Double.parseDouble(values[index++]),
                StructuralVariantType.valueOf(values[index++]),
                values[index++],
                values[index++],
                Long.parseLong(values[index++]),
                Long.parseLong(values[index++]),
                Integer.parseInt(values[index++]),
                Integer.parseInt(values[index++]),
                Integer.parseInt(values[index++]),
                Integer.parseInt(values[index++]),
                Double.parseDouble(values[index++]),
                Integer.parseInt(values[index++]),
                Double.parseDouble(values[index++]),
                Integer.parseInt(values[index++]),
                values[index++],
                values[index++],
                reprocess ? "" : values[index++],
                reprocess ? "" : values[index++]);
         */
    }

    private void writeCsv(final SvData svData)
    {
        try
        {
            if(mCsvWriter == null)
            {
                String outputFileName = !mConfig.VcfFile.isEmpty() ?
                        CsvFileWriter.generateFilename(mConfig.OutputDir, mConfig.SampleId)
                        : mConfig.OutputDir + "LNX_GERMLINE_SVS.csv";

                mCsvWriter = createBufferedWriter(outputFileName, false);
                mCsvWriter.write(CsvFileWriter.header());

                mCsvWriter.newLine();
            }

            mCsvWriter.write(CsvFileWriter.toString(svData));

            mCsvWriter.newLine();
        }
        catch (final IOException e)
        {
            GM_LOGGER.error("error writing CSV output file: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mCsvWriter);

    }


}
