package com.hartwig.hmftools.compar.cobalt;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

public class RawCobaltRatioFile
{
    static List<RawCobaltRatio> read(final String filename)
    {
        List<RawCobaltRatio> ratios = new ArrayList<>();
        try(DelimFileReader reader = new DelimFileReader(filename))
        {
            Integer chrIndex = reader.getColumnIndex(CobaltRatioFile.Column.chromosome);
            Integer posIndex = reader.getColumnIndex(CobaltRatioFile.Column.position);
            Integer refReadCountIndex = reader.getColumnIndex(CobaltRatioFile.Column.referenceReadDepth);
            Integer tumorReadCountIndex = reader.getColumnIndex(CobaltRatioFile.Column.tumorReadDepth);
            Integer refGcRatioIndex = reader.getColumnIndex(CobaltRatioFile.Column.referenceGCRatio);
            Integer tumorGcRatioIndex = reader.getColumnIndex(CobaltRatioFile.Column.tumorGCRatio);
            Integer refGcDiploidRatioIndex = reader.getColumnIndex(CobaltRatioFile.Column.referenceGCDiploidRatio);

            for(DelimFileReader.Row row : reader)
            {
                RawCobaltRatio ratio = new RawCobaltRatio(
                        row.get(chrIndex),
                        row.getInt(posIndex),
                        row.getDouble(refReadCountIndex),
                        row.getDouble(tumorReadCountIndex),
                        row.getDouble(refGcRatioIndex),
                        row.getDouble(tumorGcRatioIndex),
                        row.getDouble(refGcDiploidRatioIndex));

                ratios.add(ratio);
            }
        }
        return ratios;
    }
}
