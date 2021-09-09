package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.Category.PURITY;

import java.util.List;

import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityQCContext;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.compress.utils.Lists;

public class PurityComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public PurityComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    public void processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        final MatchLevel matchLevel = mConfig.Categories.get(PURITY);

        final List<List<ComparableItem>> sourcePurities = Lists.newArrayList();

        for(String sourceName : mConfig.DbSourceNames)
        {
            sourcePurities.add(getSamplePurity(sampleId, mConfig.DbConnections.get(sourceName)));
        }

        for(int i = 0; i < mConfig.DbSourceNames.size() - 1; ++i)
        {
            final String source1 = mConfig.DbSourceNames.get(i);

            for(int j = i + 1; j < mConfig.DbSourceNames.size(); ++j)
            {
                final String source2 = mConfig.DbSourceNames.get(j);

                CommonUtils.compareItems(sampleId, mismatches, matchLevel, source1, source2, sourcePurities.get(i), sourcePurities.get(j));
            }
        }
    }

    private List<ComparableItem> getSamplePurity(final String sampleId, final DatabaseAccess dbAccess)
    {
        final PurityQCContext purityQcContext = dbAccess.readPurityContext(sampleId);

        List<ComparableItem> items = Lists.newArrayList();
        items.add(new PurityData(purityQcContext));
        return items;
    }

}
