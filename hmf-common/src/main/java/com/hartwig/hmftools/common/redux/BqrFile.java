package com.hartwig.hmftools.common.redux;

import static com.hartwig.hmftools.common.redux.ReduxCommon.REDUX_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getColumnIndex;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.ConsensusType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class BqrFile
{
    private static final Logger LOGGER = LogManager.getLogger(BqrFile.class);

    private static String FILE_EXTENSION = REDUX_FILE_ID + ".bqr.tsv";

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

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
                .add(COL_ALT)
                .add(COL_REF)
                .add(COL_TNC)
                .add(COL_READ_TYPE)
                .add(COL_COUNT)
                .add(COL_ORIG_QUAL)
                .add(COL_RECAL_QUAL)
                .toString();
    }

    private static final String COL_REF = "Ref";
    private static final String COL_ALT = "Alt";
    private static final String COL_TNC = "TrinucleotideContext";
    private static final String COL_READ_TYPE = "ReadType";
    private static final String COL_COUNT = "Count";
    private static final String COL_ORIG_QUAL = "OriginalQual";
    private static final String COL_RECAL_QUAL = "RecalibratedQual";

    public static List<BqrRecord> read(final String filename)
    {
        List<BqrRecord> counts = Lists.newArrayList();

        try
        {
            List<String> lines = Files.readAllLines(new File(filename).toPath());

            String header = lines.get(0);
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
            lines.remove(0);

            int refIndex = getColumnIndex(fieldsIndexMap, COL_REF);
            int altIndex = getColumnIndex(fieldsIndexMap, COL_ALT);
            int tncIndex = getColumnIndex(fieldsIndexMap, COL_TNC);
            int countIndex = getColumnIndex(fieldsIndexMap, COL_COUNT);
            int origQualIndex = getColumnIndex(fieldsIndexMap, COL_ORIG_QUAL);
            int recalQualIndex = getColumnIndex(fieldsIndexMap, COL_RECAL_QUAL);
            Integer readTypeIndex = getColumnIndex(fieldsIndexMap, COL_READ_TYPE);

            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIM, -1);

                byte alt = values[altIndex].getBytes()[0];
                byte ref = values[refIndex].getBytes()[0];
                String triContext = values[tncIndex];
                ConsensusType readType = readTypeIndex != null ? ConsensusType.valueOf(values[readTypeIndex]) : ConsensusType.NONE;

                int count = Integer.parseInt(values[countIndex]);
                byte origQuality = (byte)Integer.parseInt(values[origQualIndex]);
                double recalibQuality = Double.parseDouble(values[recalQualIndex]);

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
