package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;

import java.util.List;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

public class AmberSites
{
    private static final String FLD_GNOMAD_FREQ = "GnomadFreq";
    private static final String FLD_GC_RATIO = "GcRatio";

    public static List<AmberSite> loadAmberSitesFile(String path)
    {
        try(DelimFileReader reader = new DelimFileReader(path))
        {
            int chrIndex = reader.getColumnIndex(FLD_CHROMOSOME);
            int posIndex = reader.getColumnIndex(FLD_POSITION);
            int gnomadIndex = reader.getColumnIndex(FLD_GNOMAD_FREQ);
            int gcRatioIndex = reader.getColumnIndex(FLD_GC_RATIO);

            return reader.stream().map(row ->
            {
                String chromosome = row.get(chrIndex);
                int position = row.getInt(posIndex);
                double gnomadFreq = row.getDouble(gnomadIndex);
                double gcRatio = row.getDouble(gcRatioIndex);
                return new AmberSite(new BasePosition(chromosome, position), gnomadFreq, gcRatio);
            }).toList();
        }
    }
}
