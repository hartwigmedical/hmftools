package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;

import java.util.List;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Amber heterozygous sites file loading.
public class AmberSites
{
    private static final String FLD_GNOMAD_FREQ = "GnomadFreq";

    private static final Logger LOGGER = LogManager.getLogger(AmberSites.class);

    public static List<AmberSite> loadAmberSitesFile(final String path)
    {
        LOGGER.debug("Loading Amber sites file: {}", path);

        try(DelimFileReader reader = new DelimFileReader(path))
        {
            int chrIndex = reader.getColumnIndex(FLD_CHROMOSOME);
            int posIndex = reader.getColumnIndex(FLD_POSITION);
            int gnomadIndex = reader.getColumnIndex(FLD_GNOMAD_FREQ);

            List<AmberSite> sites = reader.stream().map(row ->
            {
                String chromosome = row.get(chrIndex);
                int position = row.getInt(posIndex);
                double gnomadFreq = row.getDouble(gnomadIndex);
                return new AmberSite(new BasePosition(chromosome, position), gnomadFreq);
            }).toList();

            LOGGER.info("Loaded {} Amber sites from {}", sites.size(), path);
            return sites;
        }
    }
}
