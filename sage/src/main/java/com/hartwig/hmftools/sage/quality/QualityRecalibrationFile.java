package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.sage.SageCommon.DELIM;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

public final class QualityRecalibrationFile
{
    public static String generateBqrFilename(final String sample, final String outputDir)
    {
        return outputDir + sample + ".sage.bqr.tsv";
    }

    public static void write(final String filename, final List<QualityRecalibrationRecord> counts) throws IOException
    {
        Collections.sort(counts);
        Files.write(new File(filename).toPath(), toLines(counts));
    }

    private static List<String> toLines(final Collection<QualityRecalibrationRecord> bafs)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        bafs.stream().map(QualityRecalibrationFile::toString).forEach(lines::add);
        return lines;
    }

    private static String toString(final QualityRecalibrationRecord baf)
    {
        StringJoiner sj = new StringJoiner(DELIM);
        sj.add(String.valueOf((char)baf.Key.Alt));
        sj.add(String.valueOf((char)baf.Key.Ref));
        sj.add(new String(baf.Key.TrinucleotideContext));
        sj.add(String.valueOf(baf.Count));
        sj.add(String.valueOf(baf.Key.Quality));
        sj.add(String.format("%.2f", baf.RecalibratedQuality));
        return sj.toString();
    }

    private static String header()
    {
        return new StringJoiner(DELIM, "", "")
                .add("alt")
                .add("ref")
                .add("trinucleotideContext")
                .add("count")
                .add("originalQual")
                .add("recalibratedQual")
                .toString();
    }

    public static List<QualityRecalibrationRecord> read(final String filename)
    {
        List<QualityRecalibrationRecord> counts = Lists.newArrayList();

        try
        {
            List<String> lines = Files.readAllLines(new File(filename).toPath());

            String header = lines.get(0);
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, "\t");
            lines.remove(0);

            for(String line : lines)
            {
                String[] values = line.split("\t", -1);

                byte alt = values[fieldsIndexMap.get("alt")].getBytes()[0];
                byte ref = values[fieldsIndexMap.get("ref")].getBytes()[0];
                String triContext = values[fieldsIndexMap.get("trinucleotideContext")];
                int count = Integer.parseInt(values[fieldsIndexMap.get("count")]);
                byte origQuality = (byte)Integer.parseInt(values[fieldsIndexMap.get("originalQual")]);
                double recalibQuality = Double.parseDouble(values[fieldsIndexMap.get("recalibratedQual")]);

                BaseQualityKey key = new BaseQualityKey(ref, alt, triContext.getBytes(), origQuality);

                counts.add(new QualityRecalibrationRecord(key, count, recalibQuality));
            }
        }
        catch(Exception e)
        {
            SG_LOGGER.error("failed to read BQR file({}) record index({}): {}", filename, counts.size(), e.toString());
        }

        return counts;

    }


}
