package com.hartwig.hmftools.common.qual;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class BqrFile
{
    private static final Logger LOGGER = LogManager.getLogger(BqrFile.class);

    public static void write(final String filename, final List<BqrRecord> counts) throws IOException
    {
        Collections.sort(counts);
        Files.write(new File(filename).toPath(), toLines(counts));
    }

    private static List<String> toLines(final Collection<BqrRecord> bafs)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        bafs.stream().map(BqrFile::toString).forEach(lines::add);
        return lines;
    }

    private static String toString(final BqrRecord baf)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(String.valueOf((char)baf.Key.Alt));
        sj.add(String.valueOf((char)baf.Key.Ref));
        sj.add(new String(baf.Key.TrinucleotideContext));
        sj.add(baf.Key.ReadType.toString());
        sj.add(String.valueOf(baf.Count));
        sj.add(String.valueOf(baf.Key.Quality));
        sj.add(String.format("%.2f", baf.RecalibratedQuality));
        return sj.toString();
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIM, "", "")
                .add("alt")
                .add("ref")
                .add("trinucleotideContext")
                .add("readType")
                .add("count")
                .add("originalQual")
                .add("recalibratedQual")
                .toString();
    }

    public static List<BqrRecord> read(final String filename)
    {
        List<BqrRecord> counts = Lists.newArrayList();

        try
        {
            List<String> lines = Files.readAllLines(new File(filename).toPath());

            String header = lines.get(0);
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
            lines.remove(0);

            Integer readTypeIndex = fieldsIndexMap.get("readType");

            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIM, -1);

                byte alt = values[fieldsIndexMap.get("alt")].getBytes()[0];
                byte ref = values[fieldsIndexMap.get("ref")].getBytes()[0];
                String triContext = values[fieldsIndexMap.get("trinucleotideContext")];
                BqrReadType readType = readTypeIndex != null ? BqrReadType.valueOf(values[readTypeIndex]) : BqrReadType.NONE;

                int count = Integer.parseInt(values[fieldsIndexMap.get("count")]);
                byte origQuality = (byte)Integer.parseInt(values[fieldsIndexMap.get("originalQual")]);
                double recalibQuality = Double.parseDouble(values[fieldsIndexMap.get("recalibratedQual")]);

                BqrKey key = new BqrKey(ref, alt, triContext.getBytes(), origQuality, readType);

                counts.add(new BqrRecord(key, count, recalibQuality));
            }
        }
        catch(Exception e)
        {
            LOGGER.error("failed to read BQR file({}) record index({}): {}", filename, counts.size(), e.toString());
            return null;
        }

        return counts;
    }
}
