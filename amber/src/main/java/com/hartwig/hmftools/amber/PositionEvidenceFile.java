package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

public final class PositionEvidenceFile
{

    private static final String TUMOR_EXTENSION = ".amber.tumor.raw.tsv.gz";

    public static String generateTumorDataFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + TUMOR_EXTENSION;
    }

    private static final String CHROMOSOME = "Chromosome";
    private static final String POSITION = "Position";
    private static final String REF = "Ref";
    private static final String ALT = "Alt";
    private static final String READ_DEPTH = "ReadDepth";
    private static final String INDEL_COUNT = "IndelCount";
    private static final String REF_COUNT = "RefCount";
    private static final String ALT_COUNT = "AltCount";
    private static final String BASE_QUAL_FILTERED = "BaseQualFiltered";
    private static final String MAP_QUAL_FILTERED = "MapQualFiltered";

    public static void write(final String filename, final List<PositionEvidence> bafs) throws IOException
    {
        DelimFileWriter.write(filename, columns(), bafs, (baf, row) ->
        {
            row.set(CHROMOSOME, baf.Chromosome);
            row.set(POSITION, baf.Position);
            row.set(REF, baf.Ref.name());
            row.set(ALT, baf.Alt.name());
            row.set(READ_DEPTH, baf.ReadDepth);
            row.set(INDEL_COUNT, baf.IndelCount);
            row.set(REF_COUNT, baf.RefSupport);
            row.set(ALT_COUNT, baf.AltSupport);
            row.set(BASE_QUAL_FILTERED, baf.BaseQualFiltered);
            row.set(MAP_QUAL_FILTERED, baf.MapQualFiltered);
        });
    }

    public static List<PositionEvidence> read(String filename)
    {
        List<PositionEvidence> result = new ArrayList<>();
        try(DelimFileReader reader = new DelimFileReader(filename))
        {
            for(DelimFileReader.Row row : reader)
            {
                PositionEvidence pe = new PositionEvidence(row.get(0), row.getInt(1), row.get(2), row.get(3));
                pe.ReadDepth = row.getInt(4);
                pe.IndelCount = row.getInt(5);
                pe.RefSupport = row.getInt(6);
                pe.AltSupport = row.getInt(7);
                pe.BaseQualFiltered = row.getInt(8);
                pe.MapQualFiltered = row.getInt(9);
                result.add(pe);
            }
        }
        return result;
    }

    private static List<String> columns()
    {
        return List.of(CHROMOSOME, POSITION, REF, ALT, READ_DEPTH, INDEL_COUNT, REF_COUNT, ALT_COUNT, BASE_QUAL_FILTERED, MAP_QUAL_FILTERED);
    }
}
