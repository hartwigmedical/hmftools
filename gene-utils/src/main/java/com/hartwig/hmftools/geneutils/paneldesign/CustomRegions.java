package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

public class CustomRegions
{
    private final PanelConfig mConfig;
    private final PanelCache mPanelCache;

    public CustomRegions(final PanelConfig config, final PanelCache panelCache)
    {
        mConfig = config;
        mPanelCache = panelCache;
    }

    public void run()
    {
        if(mConfig.CustomRegionFile != null)
            loadRegions(mConfig.CustomRegionFile);
    }

    // load file of the form: Chromosome, Position
    enum Column
    {
        Chromosome,
        PositionStart,
        PositionEnd,
        SourceInfo;
    }

    private void loadRegions(final String filename)
    {
        DelimFileReader reader = new DelimFileReader(filename);

        int regionCount = 0;

        for(DelimFileReader.Row row : reader)
        {
            String chromosome = row.get(Column.Chromosome);
            int positionStart = row.getInt(Column.PositionStart);
            int positionEnd = row.getInt(Column.PositionEnd);
            String sourceInfo = row.get(Column.SourceInfo);

            PanelRegion region = new PanelRegion(new ChrBaseRegion(chromosome, positionStart, positionEnd), RegionType.CUSTOM, sourceInfo);
            mPanelCache.addRegion(region);
            ++regionCount;
        }

        GU_LOGGER.info("loaded {} custom panel regions", regionCount);
    }
}
